import os
import csv
import sys
import time
import pandas as pd
from flask import Flask, request, jsonify, render_template, redirect, url_for, session, send_from_directory, after_this_request
from flask_session import Session
import logging
import traceback
import json
from datetime import datetime, timedelta
from dotenv import load_dotenv
from werkzeug.security import check_password_hash
import psycopg2
from database.query_db import search_postgres, filter_multiple_postgres
from primer_design.primer3 import *
from primer_design.helper_function import generate_bed, vcf_to_bed, get_postgres_connection, prepare_df
from werkzeug.middleware.proxy_fix import ProxyFix
load_dotenv()
logging.basicConfig(level=logging.DEBUG)

app = Flask(__name__, static_folder='static')
#app.config["APPLICATION_ROOT"] = "/gstt_primer_design"
app.wsgi_app = ProxyFix(app.wsgi_app, x_proto=1, x_host=1, x_prefix=1)
DB_HOST = os.environ["DB_HOST"]
DB_NAME = os.environ["DB_NAME"]
DB_USER = os.environ["DB_USER"]
DB_PASSWORD = os.environ["DB_PASSWORD"]
default_schema = "primer_tool"
default_table = "ordered_primers"
primer_editable_columns = ['notes', 'passedvalidation', 'mix', 'dilute_time',
                           'tray', 'freezer', 'grid_fw', 'grid_rv']


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
app.logger.info(app.config["SECRET_KEY"])


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
            session['username'] = user
            app.logger.info(f"User logged in: {user}")
            return redirect(url_for('index'))
        else:
            return render_template('gstt_primer_design.html', error="Invalid username or password!")
    return render_template('gstt_primer_design.html')


@app.route("/protected")
def protected():
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
    # stored on server only, not sent to browser
    user = session["username"]
    password = session["password"]

    return f"user: {user}, password: {password}"


@app.route('/index')
def index():
    """
    load index page if log in is successful
    """
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
    username = session.get('username')
    return render_template('index.html', username=username)


@app.route('/igv_view/<genome>')
def igv_view(genome):
    """
    Load igv_view to visualize primers and common snp
    Generated primers are saved in bed file
    Common SNP are filtered using primer bed file
    Both generated primers and common SNP are loaded onto igv_view
    Genome can either be GRCh 37 or 38
    """
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
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
        else:
            build = 38
        # check if the primers are for build37 or 38
        temp_df = df_combined[df_combined["GRCh"] == build]
        if temp_df["GRCh"].iloc[0] == 38:
            temp_df["chr"] = "chr" + temp_df["chr"].astype(str)
        temp_df = temp_df.drop(columns=["GRCh"])
        # put primers into bed format
        bed_file = generate_bed(temp_df, genome[2:])
        #generate_filtered_vcf(bed_file, genome[2:])
        # filter common SNP using primer bed
        vcf_to_bed(bed_file, genome[2:])
    # define genome and initial focus for igv_view
    initial_query = {
                    "genome": genome,
                    "locus": "chr" + str(df["chr"][0]) + ":" + str(df["start_POS"][0])
                    }

    return render_template('igv_view.html', initial_query=initial_query)


@app.route('/design_primer', methods=['GET', 'POST'])
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
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))

    if request.method == "POST":
        # request.form contains multiple rows as lists
        # Each input field has name="chr[]", "position[]", etc.
        chr = request.form.getlist("chr[]")
        variant = request.form.getlist("variant[]")
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
                "variant": variant[i],
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
        csv_file = "temp_input.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()  # write the header row
            writer.writerows(rows)  # write all rows
        empty_chr_rows = [row for row in rows if not row["chr"]]

        if csv_file and csv_file.endswith(".csv") and not empty_chr_rows:
            app_datetimestr = datetime.now().strftime("%Y%m%d%H%M%S%f")
            order_primer = GeneratePrimer(app_datetimestr)
            output = order_primer.parse_input(csv_file)
            app.logger.info("Primer design done")
            session['primer_output'] = output.to_dict(orient='records')
            session["app_datetimestr"] = app_datetimestr
            session["order_sheet_name"] = f'primer_order_sheet_TEST_VERSION_{session["app_datetimestr"]}.csv'
            del_file(["temp_input.csv"])
            return redirect(url_for('success_primer_design'))

        else:
            return render_template("missing_input.html")

    return render_template('design_primer.html')


@app.route('/query', methods=['GET', 'POST'])
def query_data():
    """
    Function to query postgres tables
    """
    default_query_grch = 37
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
    if request.method == 'POST':
        schema = request.form.get('schema') or default_schema
        table = request.form.get('table') or default_table
        chromosome = request.form.get('chr')
        value1 = request.form.get('value1')
        gene = request.form.get('gene')
        value2 = request.form.get('value2')
        variant = request.form.get('variant')
        value3 = request.form.get('value3')
        passedvalidation = request.form.get('passedvalidation')
        value4 = request.form.get('value4')
        value5 = request.form.get('value5')
        value6 = request.form.get('value6')
        grch = request.form.get('grch')
        value7 = request.form.get('value7')
        file = request.files.get('file')

        if not table or not schema:
            return "Schema and Table name are required!", 400
        if not any([value1, value2, value3, value4, value5, value6, value7, file]):
            return render_template(
                "query_result.html",
                results=None,
                editable_columns=None,
                error="Please provide at least one search filter."
            )
        try:
            if not file:
                result_list = search_postgres(
                    DB_NAME, DB_USER, DB_PASSWORD, DB_HOST,
                    schema, table, chromosome, gene, variant, passedvalidation,
                    value1, value2, value3, value4, grch, value7,
                    value5, value6
                )
            else:
                filters = json.load(file)
                result_list = filter_multiple_postgres(
                    DB_NAME, DB_USER, DB_PASSWORD,
                    schema, table, filters
                )
            editable_columns = primer_editable_columns
            # Successful query
            return render_template(
                "query_result.html",
                results=result_list,
                editable_columns=editable_columns,
                error=None
            )
        except Exception as e:
            # Failed query
            print("Error type:", type(e).__name__)
            print("Error message:", str(e))
            traceback.print_exc()

            return render_template(
                "query_result.html",
                results=None,
                editable_columns=None,
                error=str(e)
            )

    return render_template(
        'query.html',
        default_schema=default_schema,
        default_table=default_table,
        default_grch=default_query_grch
    )


@app.route('/success_primer_design')
def success_primer_design():
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
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
        df = df[['Primer_Pair', "chr", "GRCh", "variant", "start_POS", "end_POS",
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
def save_selected():

    output_dict = session.get('primer_output')
    if not output_dict:
        return "No data found in session."

    df = pd.DataFrame(output_dict)
    columns = ["chr", "start_POS", "end_POS", "variant", "GRCh", "tag_name",
               "primer", "tagged_primer"]
    rows = []
    selected_rows = request.form.getlist('selected_rows')

    if selected_rows:
        selected_df = df.iloc[[int(i) for i in selected_rows]]

        for index, row in selected_df.iterrows():
            temp_df = selected_df.loc[[index]]
            temp_df = temp_df.reset_index()
            (FW_primer, RV_primer, tagged_FW,
             tagged_RV, tag_name_FW, tag_name_R) = prepare_order_sheet(temp_df)

            # Append FW
            rows.append({
                "chr": row["chr"],
                "start_POS": row["start_POS"],
                "end_POS": row["end_POS"],
                "variant": row["variant"],
                "GRCh": row["GRCh"],   
                "tag_name": tag_name_FW,
                "primer": FW_primer,
                "tagged_primer": tagged_FW
            })

            # Append RV
            rows.append({
                "chr": row["chr"],
                "start_POS": row["start_POS"],
                "end_POS": row["end_POS"],
                "variant": row["variant"],
                "GRCh": row["GRCh"],
                "tag_name": tag_name_R,
                "primer": RV_primer,
                "tagged_primer": tagged_RV
            })
            insert_DB(temp_df, tagged_FW, tagged_RV,temp_df["order_tag"][0],
                      DB_USER, DB_PASSWORD, DB_NAME, DB_HOST)
        df_to_order = pd.DataFrame(rows, columns=columns)
        df_to_order.to_csv(f"/app/output/{session['order_sheet_name']}", index=False)
        return render_template("save_complete.html")

    return "No rows selected."


@app.route("/save_complete")
def save_complete():
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
    return render_template("save_complete.html")


@app.route("/download")
def download():
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
    csv_to_download = session["order_sheet_name"]
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


@app.route('/update-row', methods=['POST'])
def update_row():
    data = request.get_json()
    row_id = data['id']
    updated_fields = data['data']

    # Remove primary key from fields to update
    updated_fields.pop('primerid', None)

    # Only allow these columns to be updated
    allowed_columns = primer_editable_columns
    updated_fields = {k: v for k, v in updated_fields.items() if k in allowed_columns}

    if not updated_fields:
        return jsonify({"status": "error", "message": "No editable columns selected"})

    schema_name = default_schema
    table_name = default_table

    # Get user from session for logging
    user = session.get('username', 'unknown')

    # Log what is being updated
    app.logger.info(
        f"User '{user}' updated row with primerid={row_id}. "
        f"Columns changed: {list(updated_fields.keys())}. "
        f"New values: {updated_fields}"
    )

    try:
        set_clause = ", ".join([f"{col} = %s" for col in updated_fields.keys()])
        values = list(updated_fields.values())
        values.append(row_id)  # for WHERE primerid = %s

        query = f"UPDATE {schema_name}.{table_name} SET {set_clause} WHERE primerid = %s"

        conn = psycopg2.connect(
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD,
            host=DB_HOST
        )
        cur = conn.cursor()
        cur.execute(query, values)
        conn.commit()
        cur.close()
        conn.close()

        return jsonify({"status": "success"})

    except Exception as e:
        app.logger.error(f"Error updating row {row_id} by user '{user}': {str(e)}")
        return jsonify({"status": "error", "message": str(e)})


@app.route("/manual_insert", methods=["GET", "POST"])
def manual_insert():
    if request.method == "POST":
        try:
            chr = request.form.get("chr")
            start_pos = request.form.get("start_POS")
            end_pos = request.form.get("end_POS")
            variant = request.form.get("variant")
            left_seq = request.form.get("Left_Sequence")
            right_seq = request.form.get("Right_Sequence")
            left_start = request.form.get("Left_Start")
            left_end = request.form.get("Left_End")
            right_start = request.form.get("Right_Start")
            right_end = request.form.get("Right_End")
            product_size = request.form.get("Pair_Product_Size")
            gene = request.form.get("gene")
            grch = request.form.get("GRCh")
            tagged_FW = request.form.get("tagged_FW")
            tagged_RV = request.form.get("tagged_RV")
            tag = request.form.get("tag")
            notes = request.form.get("Notes")
            passval = request.form.get("PassedValidation")

            df_insert = pd.DataFrame([{
                "chr": chr,
                "start_POS": int(start_pos),
                "end_POS": int(end_pos),
                "variant": variant or "variant",
                "Left_Sequence": left_seq or "NNNNN",
                "Right_Sequence": right_seq or "NNNNN",
                "Left_Start": int(left_start) if left_start else 1,
                "Left_End": int(left_end) if left_end else 1,
                "Right_Start": int(right_start) if right_start else 1,
                "Right_End": int(right_end) if right_end else 1,
                "Pair_Product_Size": int(product_size) if product_size else 1,
                "gene": gene or "gene",
                "GRCh": int(grch)
            }])

            insert_DB(
                df_insert, tagged_FW, tagged_RV, tag,
                username=DB_USER,
                password=DB_PASSWORD,
                db_name=DB_NAME,
                db_host=DB_HOST,
                notes=notes,
                passval=passval
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

    session.clear()
    if username:
        app.logger.info(f"User logged out: {username}")
    else:
        app.logger.error("No user name. Check!!!")
    return redirect(url_for("gstt_primer_design"))


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
