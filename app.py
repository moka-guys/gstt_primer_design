import os
import csv
import sys
import time
import pandas as pd
from flask import Flask, request, render_template, redirect, url_for, session, send_from_directory, after_this_request
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

load_dotenv()
logging.basicConfig(level=logging.DEBUG)

app = Flask(__name__, static_folder='/app/static')
DB_HOST = os.environ["DB_HOST"]
DB_NAME = os.environ["DB_NAME"]
DB_USER = os.environ["DB_USER"]
DB_PASSWORD = os.environ["DB_PASSWORD"]


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
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(minutes=5)  # log out after 5 min idle
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


@app.route('/gstt_primer_design', methods=['GET', 'POST'])
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


@app.route('/igv/<genome>')
def igv(genome):
    """
    Load IGV to visualize primers and common snp
    Generated primers are saved in bed file
    Common SNP are filtered using primer bed file
    Both generated primers and common SNP are loaded onto IGV
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
        temp_df = temp_df.drop(columns=["GRCh"])
        # put primers into bed format
        bed_file = generate_bed(temp_df, genome[2:])
        #generate_filtered_vcf(bed_file, "primer_design/config.json", genome[2:])
        # filter common SNP using primer bed
        vcf_to_bed(bed_file, "primer_design/config.json", genome[2:])
    # define genome and initial focus for IGV
    initial_query = {
                    "genome": genome,
                    "locus": "chr" + str(df["chr"][0]) + ":" + str(df["POS"][0])
                    }

    return render_template('igv.html', initial_query=initial_query)


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
            order_primer = GenerateOrder(app_datetimestr)
            output = order_primer.order(csv_file, DB_USER, DB_PASSWORD, DB_NAME, DB_HOST)
            app.logger.info("Primer design done")
            session['primer_output'] = output.to_dict(orient='records')
            session["app_datetimestr"] = app_datetimestr
            session["order_sheet_name"] = f'primer_order_sheet_TEST_VERSION_{session["app_datetimestr"]}.csv'
            del_file(["temp_input.csv"])
            return redirect(url_for('success_primer_design'))

        else:
            return '''
                    Missing valid input(s) to design primers. Try again!!! 
                    <p>Go to <a href="/index">Home Page</a> or <a href="/design_primer">Design Primer</a></p>
                   '''

    return render_template('design_primer.html')


@app.route('/query', methods=['GET', 'POST'])
def query_data():
    """
    Function to query postgres tables
    """
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
    col1 = None
    col2 = None
    col3 = None
    value1 = None
    value2 = None
    value3 = None
    if request.method == 'POST':
        print("query submitted")
        schema = request.form.get('schema')
        table = request.form.get('table')
        col1 = request.form.get('col1')
        value1 = request.form.get('value1')
        col2 = request.form.get('col2')
        value2 = request.form.get('value2')
        col3 = request.form.get('col3')
        value3 = request.form.get('value3')
        file = request.files.get('file')

        if not table:
            return "File and table name are required!", 400
        try:
            if not file:
                result_list = search_postgres(DB_NAME, DB_USER, DB_PASSWORD, DB_HOST,
                                              schema, table, col1, col2, col3,
                                              value1, value2, value3)
            else:
                filters = json.load(file)
                result_list = filter_multiple_postgres(DB_NAME, DB_USER,
                                                       DB_PASSWORD, schema, table, filters)

            output = ""
            for row_dict in result_list:
                output += f"{row_dict}<br>"
            output += '<p>Go to <a href="/index">Home Page</a> or <a href="/query">Query Page</a></p>'
            return output
        except Exception as e:
            print("Error type:", type(e).__name__)
            print("Error message:", str(e))
            traceback.print_exc()
            return f'''
                    Error: {e}
                    <p>Your query is invalid. Try again!!<p>
                    <p>Go to <a href="/index">Home Page</a> or <a href="/query">Query Page</a></p>
                    '''

    return render_template('query.html')


@app.route('/success_primer_design')
def success_primer_design():
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
    """
    Success page for primer design. Generated primers are shown as table if any.
    Option to visualize primers on IGV is provided.
    Primer not found message is shown if no primer is designed
    """
    output_dict = session.get('primer_output')

    if output_dict:
        df = pd.DataFrame(output_dict)
        builds = list(df["GRCh"].unique())
        # reorder df col to appear on UI
        df = df[['Primer_Pair', "chr", "GRCh", "variant", "POS",
                 'Left_Sequence', 'Right_Sequence', 'Pair_Product_Size',
                 'Left_Start', 'Left_End', 'Right_Start', 'Right_End', 'Specificity',
                 'FW_primer_snp', 'RV_primer_snp', 'snp_validity', 'gene', 'exon_num',
                 'transcript', 'ref_genome_source', 'common_snp_source']]
        table_html = df.to_html(index=False)
        if 37 in builds and 38 in builds:
            return f'''
                    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css" rel="stylesheet">

                    <h2>Primer design completed!</h2>
                    <div class="table-responsive">
                        {table_html.replace('<table border="1"', '<table class="table table-striped table-bordered table-hover"')}
                    </div>
                    <p>Go to <a href="/index">Home Page</a> or <a href="/design_primer">Design Primer</a> or <a href="/igv/hg19">IGV build37</a> or <a href="/igv/hg38">IGV build38</a></p>
                    <p>
                        <a href="/download" class="btn btn-primary btn-lg">
                            Download Order Sheet
                        </a>
                    </p>
                    '''
        elif 37 in builds:
            return f'''
                    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css" rel="stylesheet">

                    <h2>Primer design completed!</h2>
                    <div class="table-responsive">
                        {table_html.replace('<table border="1"', '<table class="table table-striped table-bordered table-hover"')}
                    </div>
                    <p>Go to <a href="/index">Home Page</a> or <a href="/design_primer">Design Primer</a> or <a href="/igv/hg19">IGV build37</a></p>
                    <p>
                        <a href="/download" class="btn btn-primary btn-lg">
                            Download Order Sheet
                        </a>
                    </p>
                    '''
        elif 38 in builds:
            return f'''
                    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css" rel="stylesheet">

                    <h2>Primer design completed!</h2>
                    <div class="table-responsive">
                        {table_html.replace('<table border="1"', '<table class="table table-striped table-bordered table-hover"')}
                    </div>
                    <p>Go to <a href="/index">Home Page</a> or <a href="/design_primer">Design Primer</a> or <a href="/igv/hg38">IGV build38</a></p>
                    <p>
                        <a href="/download" class="btn btn-primary btn-lg">
                            Download Order Sheet
                        </a>
                    </p>
                    '''

    else:
        return '''
        <h2>No primers generated with max padding. Change config and try again.</h2>
        <p>Go to <a href="/index">Home Page</a> or <a href="/design_primer">Design Primer</a> or <a href="/logout">Log Out</a> </p>
        '''


@app.route("/download")
def download():
    if "username" not in session:
        return redirect(url_for("gstt_primer_design"))
    csv_to_download = session["order_sheet_name"]
    directory = app.config['DOWNLOAD_FOLDER']
    full_path = os.path.join(directory, csv_to_download)

    if not os.path.exists(full_path):
        app.logger.error("download file not found")
        return '''
        <h2>File to download is not found. It is either because the file is not generated or you've already downloaded it</h2>
        <p>Go to <a href="/index">Home Page</a> or <a href="/design_primer">Design Primer</a> or <a href="/logout">Log Out</a> </p>
        '''
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
