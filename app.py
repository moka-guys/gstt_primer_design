import os
import csv
import sys
import uuid
import time
from io import StringIO
from pathlib import Path
import pandas as pd
from flask import Flask, request, jsonify, g, render_template, redirect, url_for, Response, session, send_from_directory, after_this_request
from flask_session import Session
import logging
import traceback
import tempfile
from datetime import datetime, timedelta
from dotenv import load_dotenv
from werkzeug.security import check_password_hash, generate_password_hash
import psycopg2
from database.query_db import *
from database.insert_db import insert_DB
from primer_design.primer3 import *
from primer_design.helper_function import generate_bed, vcf_to_bed, get_postgres_connection, prepare_df
from werkzeug.middleware.proxy_fix import ProxyFix
from functools import wraps
import warnings
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
load_dotenv()
logging.basicConfig(level=logging.DEBUG)

app = Flask(__name__, static_folder='static')
#app.config["APPLICATION_ROOT"] = "/gstt_primer_design"
app.wsgi_app = ProxyFix(app.wsgi_app, x_proto=1, x_host=1, x_prefix=1)
DB_HOST = os.environ["DB_HOST"]
DB_NAME = os.environ["DB_NAME"]
DB_USER = os.environ["DB_USER"]
DB_PASSWORD = os.environ["DB_PASSWORD"]
batch_editable_columns = ['notes', 'passed_validation', 'mix', 'arrival_date',
                          'tray', 'freezer', 'grid_fw', 'grid_rv', 'archive', 'manufacturer']
primer_table = "primers"
batch_table = "primer_batches"
schema = "primer_tool"


def wait_for_db():
    max_retries = 30
    retry_interval = 2
    for i in range(max_retries):
        try:
            conn = get_postgres_connection(DB_NAME, DB_USER, DB_PASSWORD, DB_HOST)
            conn.close()
            print("Database connected!")
            return
        except psycopg2.OperationalError:
            print(f"Waiting for database... ({i+1}/{max_retries})")
            time.sleep(retry_interval)
    raise Exception("Could not connect to database")


wait_for_db()

app.config["SECRET_KEY"] = os.environ["SECRET_KEY"]
app.config["SESSION_TYPE"] = "filesystem"
app.config["SESSION_FILE_DIR"] = os.environ["SESSION_FILE_DIR"]
app.config["SESSION_PERMANENT"] = True
#app.config["SESSION_COOKIE_EXPIRES"] = None
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(minutes=30)  # log out after 30 min idle
app.config['DOWNLOAD_FOLDER'] = os.environ["DOWNLOAD_FOLDER"]
Session(app)


# logging
app.logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler("/app/logs/primer_app.log")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(formatter)
app.logger.addHandler(file_handler)
app.logger.info("xxxxxxxxKEYSxxxxxx")
app.logger.info(DB_NAME)
app.logger.info(DB_USER)
app.logger.info(DB_PASSWORD)
app.logger.info(DB_HOST)


# Print stdout/stderr to logger
class StreamToLogger:
    def __init__(self, logger, level=logging.INFO):
        self.logger = logger
        self.level = level

    def write(self, message):
        message = message.strip()
        if message:
            self.logger.log(self.level, message)

    def flush(self):
        pass  # for compatibility


sys.stdout = StreamToLogger(app.logger, logging.INFO)
sys.stderr = StreamToLogger(app.logger, logging.ERROR)


@app.route('/', methods=['GET', 'POST'])
def gstt_primer_design():
    """
    app log in with user name and password
    """
    if request.method == 'POST':
        # Get username and password from UI input
        user = request.form.get('user')
        password = request.form.get('password')
        # check if user name and password are correct matches the ones in DB
        conn = get_postgres_connection(DB_NAME, DB_USER, DB_PASSWORD, DB_HOST)
        cur = conn.cursor()
        cur.execute("SELECT password_hash FROM app_users WHERE username = %s", (user,))
        result = cur.fetchone()
        conn.close()
        if result and check_password_hash(result[0], password):
            session["username"] = user
            app.logger.info(f"User logged in: {user}")
            return redirect(url_for('index'))
        else:
            return render_template('gstt_primer_design.html', error="Invalid username or password!")
    return render_template('gstt_primer_design.html')


def login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if "username" not in session:
            return redirect(url_for("gstt_primer_design"))

        g.user = session["username"]

        return f(*args, **kwargs)
    return decorated_function


@app.route('/index')
@login_required
def index():
    """
    load index page if log in is successful
    """
    return render_template('index.html', username=g.user)


@app.route('/moka_index')
@login_required
def moka_index():
    """
    load index page for moka query
    """
    return render_template('moka_index.html')


@app.route('/change_password', methods=['GET', 'POST'])
@login_required
def change_password():

    if request.method == "POST":

        current_password = request.form.get("current_password")
        new_password = request.form.get("new_password")
        confirm_password = request.form.get("confirm_password")

        if new_password != confirm_password:
            return render_template(
                "change_password.html",
                error="New passwords do not match."
            )

        conn = get_postgres_connection(
            DB_NAME,
            DB_USER,
            DB_PASSWORD,
            DB_HOST
        )

        cur = conn.cursor()

        # Get current password hash
        cur.execute(
            "SELECT password_hash FROM app_users WHERE username = %s",
            (g.user,)
        )

        result = cur.fetchone()

        if result is None:
            conn.close()
            return render_template(
                "change_password.html",
                error="User not found."
            )

        # Verify current password
        if not check_password_hash(result[0], current_password):
            conn.close()
            return render_template(
                "change_password.html",
                error="Current password is incorrect."
            )

        # Hash the new password
        new_hash = generate_password_hash(new_password)
        print(new_hash)

        # Update password
        cur.execute(
            """
            UPDATE app_users
            SET password_hash = %s
            WHERE username = %s
            """,
            (new_hash, g.user)
        )

        conn.commit()
        conn.close()

        return render_template(
            "change_password.html",
            success="Password changed successfully."
        )

    return render_template("change_password.html")

@app.route('/igv_view/<genome>')
@login_required
def igv_view(genome):
    """
    Load igv_view to visualize primers and common snp
    Generated primers are saved in bed file
    Common SNP are filtered using primer bed file
    Both generated primers and common SNP are loaded onto igv_view
    Genome can either be GRCh 37 or 38
    """
    # get primer output from design_primer
    primer_output = session.get('primer_output')
    if primer_output:
        df = pd.DataFrame(primer_output)
        # put left and right primers into bed file format
        df_left = prepare_df(df, "Left")
        df_right = prepare_df(df, "Right")
        df_combined = pd.concat([df_left, df_right], axis=0)
        df_combined = df_combined.drop(columns=["Primer_Pair"])
        if int(genome[2:]) == 19:
            build = 37
            prefix = "chr"
        else:
            build = 38
            prefix = ""
        # check if the primers are for build37 or 38
        temp_df = df_combined[df_combined["GRCh"] == build]
        if temp_df["GRCh"].iloc[0] == 38:
            temp_df["chr"] = "chr" + temp_df["chr"].astype(str)
        temp_df = temp_df.drop(columns=["GRCh"])
        # put primers into bed format
        job_id = uuid.uuid4().hex
        bed_file = generate_bed(temp_df, genome[2:], job_id)
        session["primer_bed"] = bed_file
        session.setdefault("temp_files", [])
        session["temp_files"].append(bed_file)
        # filter common SNP using primer bed
        snp_bed_file, vcf_temp_file = vcf_to_bed(bed_file, genome[2:], job_id)
        session["snp_bed"] = snp_bed_file
        session["temp_files"].append(snp_bed_file)
        session["temp_files"].append(vcf_temp_file)
    # define genome and initial focus for igv_view
        initial_query = {
                        "genome": genome,
                        "locus": f"{prefix}{temp_df.iloc[0]['chr']}:{temp_df.iloc[0]['start']}"
                        }

    return render_template('igv_view.html', initial_query=initial_query,
                           primer_bed=session["primer_bed"],
                           snp_bed=session["snp_bed"])


@app.route('/design_primer', methods=['GET', 'POST'])
@login_required
def design_primer():
    """
    Design primers using GenerateOrder from primer_design.primer3.py
    Input params are obtained from UI
    User can key in params for multiple ROI
    Keyed in values on UI are taken as list and saved in temp csv file 
    Temp csv file is checked by pydantic to confirm validity of input values
    Valid input is used to generate primers and temp csv file is deleted
    All generated primers are printed out as table on UI
    """
    cleanup_old_files("/app/static/temp", 18000)
    cleanup_old_files("/app/output", 18000)
    if request.method == "POST":
        # request.form contains multiple rows as lists
        # Each input field has name="chr[]", "position[]", etc.
        chr = request.form.getlist("chr[]")
        primer_name = request.form.getlist("primer_name[]")
        pos_start = request.form.getlist("pos_start[]")
        pos_end = request.form.getlist("pos_end[]")
        build = request.form.getlist("build[]")
        tag = request.form.getlist("tag[]")
        opt_tm = request.form.getlist("opt_tm[]")
        min_tm = request.form.getlist("min_tm[]")
        max_tm = request.form.getlist("max_tm[]")
        opt_gc = request.form.getlist("opt_gc[]")
        min_gc = request.form.getlist("min_gc[]")
        max_gc = request.form.getlist("max_gc[]")
        min_product_size = request.form.getlist("min_product_size[]")
        max_product_size = request.form.getlist("max_product_size[]")
        opt_primer_size = request.form.getlist("opt_primer_size[]")
        min_primer_size = request.form.getlist("min_primer_size[]")
        max_primer_size = request.form.getlist("max_primer_size[]")

        # Combine into rows
        rows = []
        for i in range(len(chr)):
            row = {
                "chr": chr[i],
                "primer_name": primer_name[i],
                "pos_start": pos_start[i],
                "pos_end": pos_end[i],
                "build": build[i],
                "tag": tag[i],
                "opt_tm": opt_tm[i],
                "min_tm": min_tm[i],
                "max_tm": max_tm[i],
                "opt_gc": opt_gc[i],
                "min_gc": min_gc[i],
                "max_gc": max_gc[i],
                "min_product_size": min_product_size[i],
                "max_product_size": max_product_size[i],
                "opt_primer_size": opt_primer_size[i],
                "min_primer_size": min_primer_size[i],
                "max_primer_size": max_primer_size[i]
            }
            rows.append(row)

        # Get the fieldnames from the first dictionary
        fieldnames = rows[0].keys()
        # save into temp csv
        temp = tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".csv",
            prefix="primer_input_",
            delete=False
        )

        csv_file = temp.name

        with temp:
            writer = csv.DictWriter(
                temp,
                fieldnames=fieldnames
            )
            writer.writeheader()
            writer.writerows(rows)
        empty_chr_rows = [row for row in rows if not row["chr"]]

        if csv_file and csv_file.endswith(".csv") and not empty_chr_rows:
            app_datetimestr = datetime.now().strftime("%Y%m%d%H%M%S%f")
            random_uuid = uuid.uuid4()
            session['primer_input'] = rows
            order_primer = GeneratePrimer(app_datetimestr)
            output = order_primer.parse_input(csv_file)
            app.logger.info("Primer design done")
            session['primer_output'] = output.to_dict(orient='records')
            session["order_sheet_name_auto"] = f'primer_order_sheet_TEST_VERSION_{random_uuid}_{app_datetimestr}.csv'
            del_file([csv_file])
            app.logger.info(
                            f"User '{g.user}' designed primer for {rows}"
                        )
            return redirect(url_for('success_primer_design'))

        else:
            return render_template("missing_input.html")

    return render_template('design_primer.html', saved_rows=[])


@app.route("/modify_primer", methods=["GET", "POST"])
def modify_primer():
    if request.method == "POST":
        return design_primer()

    return render_template(
        "design_primer.html",
        saved_rows=session.get("primer_input", [])
    )


@app.route('/query', methods=['GET', 'POST'])
@login_required
def query_data():
    if request.method == 'POST':
        chr_val = request.form.get('chr')
        gene_val = request.form.get('gene')
        primer_name = request.form.get('primer_name')
        primer_id = request.form.get('primer_id')
        passed_validation = request.form.get('passed_validation')
        archive = request.form.get('archive')
        grch = request.form.get('grch')
        tray = request.form.get('tray')
        start = request.form.get('start')
        end = request.form.get('end')
        notes = request.form.get('notes')
        start_date = request.form.get('start_date')
        end_date = request.form.get('end_date')
        variant_pos = request.form.get('variant_pos')
        # ensure at least one filter exists
        if not any([chr_val, gene_val, primer_name, primer_id, passed_validation,
                    grch, tray, start, end, archive, notes, start_date, end_date,
                    variant_pos]):
            return render_template(
                "query_result.html",
                results=None,
                editable_columns=None,
                error="Please provide at least one filter."
            )
        try:
            result_list, msg = search_postgres(
                DB_NAME, DB_USER, DB_PASSWORD, DB_HOST,
                chr_val, gene_val, primer_name, primer_id, passed_validation, grch,
                tray, archive, notes, variant_pos, start, end, start_date, end_date
            )
            session["query_results"] = result_list
            return render_template(
                "query_result.html",
                results=result_list,
                msg=msg,
                editable_columns=batch_editable_columns,
                error=None
            )
        except Exception as e:
            traceback.print_exc()
            return render_template("query_result.html", error=str(e))

    return render_template('query.html')


@app.route('/query_moka_id', methods=['GET', 'POST'])
@login_required
def query_moka_id():
    # Display the search form
    if request.method == 'GET':
        return render_template('query_moka_id.html')

    # Get and validate the input
    primer_id = request.form.get('primer_id', '').strip()

    if not primer_id:
        return render_template(
            'query_moka_id_result.html',
            results=[],
            error='Please enter a Primer ID.'
        )

    try:
        # Query the database
        result_list = query_moka_by_id(
            DB_NAME,
            DB_USER,
            DB_PASSWORD,
            DB_HOST,
            primer_id
        )

        # No records found
        if not result_list:
            return render_template(
                'query_moka_id_result.html',
                results=[],
                msg='No matching primers found.'
            )

        # Display results
        return render_template(
            'query_moka_id_result.html',
            results=result_list
        )

    except Exception:
        app.logger.exception("Error while querying the MOKA database.")
        return render_template(
            'query_moka_id_result.html',
            results=[],
            error='An unexpected error occurred while querying the database. Please try again later.'
        )


@app.route('/success_primer_design')
@login_required
def success_primer_design():
    """
    Success page for primer design. Generated primers are shown as table if any.
    Option to visualize primers on igv_view is provided.
    Primer not found message is shown if no primer is designed
    """
    output_dict = session.get('primer_output')

    if output_dict:
        df = pd.DataFrame(output_dict)
        builds = list(df["GRCh"].unique())
        # reorder df col to appear on UI
        df = df[['Primer_Pair', "chr", "GRCh", "primer_name", "start_POS", "end_POS",
                 'Left_Sequence', 'Right_Sequence', 'Pair_Product_Size',
                 'Left_Start', 'Left_End', 'Right_Start', 'Right_End', 'Specificity',
                 'FW_primer_snp', 'RV_primer_snp', 'snp_validity', 'gene', 'exon_num',
                 'transcript']]
        df = df.reset_index(drop=True)
        table_data = df.to_dict(orient="records")
        return render_template(
            "success_primer_design.html",
            table_data=table_data,
            builds=builds
        )

    else:
        return render_template("no_primers.html")


@app.route('/save_selected', methods=['POST'])
@login_required
def save_selected():
    output_dict = session.get('primer_output')
    if not output_dict:
        return "No data found in session."

    df = pd.DataFrame(output_dict)
    columns = ["chr", "primer_name", "tag_name",
               "primer", "tagged_primer"]
    rows = []
    selected_rows = request.form.getlist('selected_rows')

    if selected_rows:
        selected_df = df.iloc[[int(i) for i in selected_rows]]

        for index, row in selected_df.iterrows():
            temp_df = selected_df.loc[[index]]
            temp_df = temp_df.reset_index()
            temp_df["designer"] = "A"
            (FW_primer, RV_primer, tagged_FW,
             tagged_RV, tag_name_FW, tag_name_R) = prepare_order_sheet(temp_df)

            # Append FW
            rows.append({
                "chr": row["chr"],
                "primer_name": row["primer_name"],
                "tag_name": tag_name_FW,
                "primer": FW_primer,
                "tagged_primer": tagged_FW
            })

            # Append RV
            rows.append({
                "chr": row["chr"],
                "primer_name": row["primer_name"],
                "tag_name": tag_name_R,
                "primer": RV_primer,
                "tagged_primer": tagged_RV
            })
            inserted_ids = insert_DB(temp_df, tagged_FW,
                                     tagged_RV,temp_df["order_tag"][0],
                                     DB_USER, DB_PASSWORD, DB_NAME, DB_HOST)
            app.logger.info(
                            f"User '{g.user}' selected designed primers to insert DB for upi {inserted_ids}"
                            )
        df_to_order = pd.DataFrame(rows, columns=columns)
        df_to_order["Scale"] = "25RR"
        df_to_order["Purification"] = "STD"
        df_to_order.to_csv(f"/app/output/{session['order_sheet_name_auto']}", index=False)
        session.setdefault("order_files", [])
        session["order_files"].append(session["order_sheet_name_auto"])

        return render_template("save_complete.html")

    return "No rows selected."


@app.route("/save_complete")
@login_required
def save_complete():
    return render_template("save_complete.html")


@app.route("/download/<source>")
@login_required
def download(source):
    csv_to_download = None
    if source == "manual":
        csv_to_download = session.get("order_sheet_name_manual")
    else:
        csv_to_download = session.get("order_sheet_name_auto")
    if csv_to_download is None:
        return "No order sheet available for download.", 400
    directory = app.config['DOWNLOAD_FOLDER']
    full_path = os.path.join(directory, csv_to_download)

    if not os.path.exists(full_path):
        app.logger.error("download file not found")
        return render_template("download_not_found.html")
    else:
        app.logger.info("attempt to download")

        response = send_from_directory(
            directory=directory,
            path=csv_to_download,
            as_attachment=True
        )

        # Cleanup after sending
        @after_this_request
        def cleanup(resp):
            try:
                if os.path.exists(full_path):
                    os.remove(full_path)
                    print("Deleted file:", full_path)
            except Exception as e:
                print("Error deleting file:", e)
            return resp
    return response


@app.route('/update_row', methods=['POST'])
@login_required
def update_row():
    data = request.get_json()
    row_id = data['id']
    updated_fields = data['data']

    table_type = data.get("table_type", primer_table)
    # remove non-editable keys from payload
    updated_fields.pop('upi', None)
    updated_fields.pop('primer_id', None)

    try:
        if table_type == primer_table:
            return jsonify({
                "status": "error",
                "message": "Primer table cannot be updated"
            })

        elif table_type == "batches":
            allowed_columns = batch_editable_columns
            pk = "primer_id"
            table_name = f"{schema}.{batch_table}"

        # keep only allowed columns
        updated_fields = {
            k: v for k, v in updated_fields.items()
            if k in allowed_columns
        }

        if not updated_fields:
            return jsonify({
                "status": "error",
                "message": "No editable columns selected"
            })

        conn = psycopg2.connect(
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD,
            host=DB_HOST
        )
        cur = conn.cursor()
        # get existing value
        cur.execute(f"SELECT * FROM {table_name} WHERE {pk} = %s", (row_id,))
        current_row = cur.fetchone()

        if not current_row:
            return jsonify({
                "status": "error",
                "message": "Row not found"
            })

        columns = [desc[0] for desc in cur.description]
        current_dict = dict(zip(columns, current_row))
        # get changes
        real_changes = {}

        for k, v in updated_fields.items():
            db_val = current_dict.get(k)

            if str(db_val) != str(v):
                real_changes[k] = v

        if not real_changes:
            return jsonify({
                "status": "error",
                "message": "No actual changes detected"
            })

        # update query
        set_clause = ", ".join([f"{col} = %s" for col in real_changes.keys()])
        values = list(real_changes.values())
        values.append(row_id)

        query = f"""
        UPDATE {table_name}
        SET {set_clause}
        WHERE {pk} = %s
        """

        cur.execute(query, values)
        conn.commit()
        if "query_results" in session:
            for row in session["query_results"]:
                if str(row["primer_id"]) == str(row_id):
                    row.update(real_changes)
                    break

            session.modified = True
        # log for the changes
        updated_str = ", ".join(
            f"{col}={repr(val)}" for col, val in real_changes.items()
        )

        app.logger.info(
            f"User '{g.user}' updated primer_id {row_id}: {updated_str}"
        )

        cur.close()
        conn.close()

        return jsonify({"status": "success"})

    except Exception as e:
        app.logger.error(f"Update error: {str(e)}")
        return jsonify({"status": "error", "message": str(e)})


@app.route("/export_csv")
def export_csv():
    rows = session["query_results"]
    output = StringIO()

    if rows:
        writer = csv.DictWriter(output, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return Response(
        output.getvalue(),
        mimetype="text/csv",
        headers={
            "Content-Disposition": "attachment; filename=exported_primers.csv"
        }
    )


@app.route("/manual_insert", methods=["GET", "POST"])
@login_required
def manual_insert():
    datetimestr = datetime.now().strftime("%Y%m%d%H%M%S%f")
    random_uuid = uuid.uuid4()

    if request.method == "POST":
        try:
            chr = request.form.get("chr")
            primer_name = request.form.get("primer_name")
            left_seq = request.form.get("Left_Sequence")
            right_seq = request.form.get("Right_Sequence")
            left_start = request.form.get("Left_Start")
            left_end = request.form.get("Left_End")
            right_start = request.form.get("Right_Start")
            right_end = request.form.get("Right_End")
            product_size = request.form.get("Pair_Product_Size")
            gene = request.form.get("gene")
            grch = request.form.get("GRCh")
            tag = request.form.get("tag")
            notes = request.form.get("Notes")

            df_insert = pd.DataFrame([{
                "chr": chr,
                "start_POS": 1,
                "end_POS": 1,
                "primer_name": primer_name or "primer_name",
                "Left_Sequence": left_seq or "NNNNN",
                "Right_Sequence": right_seq or "NNNNN",
                "Left_Start": int(left_start) if left_start else 1,
                "Left_End": int(left_end) if left_end else 1,
                "Right_Start": int(right_start) if right_start else 1,
                "Right_End": int(right_end) if right_end else 1,
                "Pair_Product_Size": int(product_size) if product_size else 1,
                "gene": gene or "gene",
                "GRCh": int(grch),
                "order_tag": tag,
                "designer": "M"
            }])

            inserted_ids = insert_DB(
                                df_insert, None, None, tag,
                                username=DB_USER,
                                password=DB_PASSWORD,
                                db_name=DB_NAME,
                                db_host=DB_HOST,
                                notes=notes
                            )
            app.logger.info(f"{g.user} inserted primer manually: upi {inserted_ids} ")
            # generate order sheet for manual insert primer
            (FW_primer, RV_primer, tagged_FW,
             tagged_RV, tag_name_FW, tag_name_RV) = prepare_order_sheet(df_insert)
            rows = []
            if FW_primer != "NNNNN":
                rows.append({
                    "primer": FW_primer,
                    "tagged_primer": tagged_FW,
                    "tag_name": tag_name_FW,
                })

            if RV_primer != "NNNNN":
                rows.append({
                    "primer": RV_primer,
                    "tagged_primer": tagged_RV,
                    "tag_name": tag_name_RV,
                })

            if rows:
                # repeat base dataframe rows to match number of primers
                df_new = df_insert[["chr", "primer_name"]].iloc[
                    df_insert.index.repeat(len(rows))
                ].reset_index(drop=True)

                # merge primer-specific data
                df_rows = pd.DataFrame(rows)

                df_new = pd.concat(
                    [df_new.reset_index(drop=True), df_rows.reset_index(drop=True)],
                    axis=1
                )

                df_new["Scale"] = "25RR"
                df_new["Purification"] = "STD"

                if "order_sheet_name_manual" not in session:
                    session["order_sheet_name_manual"] = (
                        f"primer_order_sheet_TEST_VERSION_{random_uuid}_{datetimestr}.csv"
                    )

                filepath = os.path.join(
                    app.config["DOWNLOAD_FOLDER"],
                    session["order_sheet_name_manual"]
                )

                if os.path.exists(filepath):
                    df_new.to_csv(
                        filepath,
                        mode="a",
                        header=False,
                        index=False
                    )
                else:
                    df_new.to_csv(
                        filepath,
                        index=False
                    )
            # Success
            return render_template("insert_success.html")

        except Exception as e:
            # Error
            return render_template("insert_error.html", error_message=str(e))

    return render_template("manual_insert.html")


@app.route("/logout")
def logout():
    """
    Function to log out app
    """
    username = session.get('username')

    delete_session_files(
        "/app/static/temp",
        "temp_files"
    )

    delete_session_files(
        "/app/output",
        "order_files"
    )
    session.clear()
    if username:
        app.logger.info(
            f"User logged out: {username}"
        )
    else:
        app.logger.error(
            "No user name. Check!!!"
        )

    return redirect(url_for("gstt_primer_design"))


def delete_session_files(directory, session_key):
    """
    Delete files stored in session under session_key
    """
    folder = Path(directory)
    files = session.get(session_key, [])

    for filename in files:
        file_path = folder / filename

        if file_path.exists():
            file_path.unlink()

            app.logger.info(
                f"Deleted temporary file: {file_path}"
            )


def cleanup_old_files(directory, max_age_seconds):

    folder = Path(directory)

    if not folder.exists():
        return

    now = time.time()

    for file in folder.iterdir():

        if file.is_file():

            age = now - file.stat().st_mtime

            if age > max_age_seconds:
                file.unlink()

                app.logger.info(
                    f"Deleted expired file: {file}"
                )


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
