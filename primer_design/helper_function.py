import os
import logging
import json
import numpy as np
import pandas as pd
import math
import subprocess
from pydantic import BaseModel, ValidationError, field_validator, model_validator
import psycopg2
from pathlib import Path


class PrimerRecord(BaseModel):
    chr: str | int
    variant: str
    pos_start: int
    pos_end: int
    build: int
    tag: str
    opt_tm: int | float
    min_tm: int | float
    max_tm: int | float
    opt_gc: int | float
    min_gc: int | float
    max_gc: int | float
    min_product_size: int
    max_product_size: int
    opt_primer_size: int
    min_primer_size: int
    max_primer_size: int

    @field_validator("pos_start", "pos_end", "min_tm", "max_tm", "opt_tm",
                     "min_gc", "max_gc", "opt_gc", "min_product_size",
                     "max_product_size", "opt_primer_size",
                     "min_primer_size", "max_primer_size")
    def non_negative(cls, v, field):
        if v < 0:
            raise ValueError(f"{field.name} must not be negative")
        return v

    @field_validator("chr")
    def validate_chr(cls, v):
        allowed = {str(i) for i in range(1, 23)} | {"X", "Y"}
        v_str = str(v).strip()
        if v_str not in allowed:
            raise ValueError("chr must be 1-22, X, or Y")
        return v_str

    @field_validator("build")
    def validate_build(cls, v):
        if v not in (37, 38):
            raise ValueError("build must be 37 or 38")
        return v

    @field_validator("tag")
    def validate_tag(cls, v):
        if v not in ("M13", "T1", "T2", "T3", "T4", "sT1", "FAM",
                     "VIC", "NED", "PET", "no_tag"):
            raise ValueError("invalid tag")
        return v

    @model_validator(mode="after")
    def check_opt_tm(cls, model):
        if not (model.min_tm <= model.opt_tm <= model.max_tm):
            raise ValueError("opt_tm must be between min_tm and max_tm")
        if model.pos_end < model.pos_start:
            raise ValueError("pos_end must be greater than or equal to pos_start")
        if not (model.min_gc <= model.opt_gc <= model.max_gc):
            raise ValueError("opt_gc must be between min_gc and max_gc")
        if not (model.min_primer_size <= model.opt_primer_size <= model.max_primer_size):
            raise ValueError("opt_primer_size must be between min_primer_size and max_primer_size")
        if model.max_product_size < model.min_product_size:
            raise ValueError("max_product_size must be greater than min_product_size")
        return model


def get_log(file_dir, datetimestr) -> logging.Logger:
    """
    Setup for primer design log file and return the logger.
    """

    logger = logging.getLogger(f"primer_log_{datetimestr}")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False  # prevent logging to root/app logger
    if not logger.hasHandlers():
        # Console handler
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # File handler
        outdir = os.path.join(file_dir, "primer_log")
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        log_file = os.path.join(outdir, f"primer_design_{datetimestr}.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_tag(config, df, tag):
    """
    get primer tag for primer order
    """
    with open(config, "r") as file:
        config = json.load(file)
    FW_primer = df["Left_Sequence"][0]
    RV_primer = df["Right_Sequence"][0]
    FW_tag = config[tag]["F"]
    RV_tag = config[tag]["R"]
    tagged_FW = FW_tag + FW_primer
    tagged_RV = RV_tag + RV_primer
    tag_name_FW = config[tag]["F_name"]
    tag_name_RV = config[tag]["R_name"]

    return FW_primer, RV_primer, tagged_FW, tagged_RV, tag_name_FW, tag_name_RV


def chr_to_int64(x):
    try:
        return pd.Int64Dtype().type(int(x))
    except ValueError:
        return x


def validate_primer_csv(df: pd.DataFrame):
    """
    check input params are valid by Pydantic
    """
    valid_records = []
    errors = []

    for idx, row in df.iterrows():
        try:
            record = PrimerRecord(**row.to_dict())
            valid_records.append(record.model_dump())
        except ValidationError as e:
            errors.append({"row": idx, "errors": e.errors()})

    return errors


def parse_csv(file):
    """
    parse input csv
    """
    df = pd.read_csv(file)
    # drop empty rows
    cols = ["chr", "pos_start", "pos_end", "build", "tag"]
    df = df[df[cols].notna().all(axis=1)]
    for c in ["pos_start", "pos_end", "build"]:
        df[c] = df[c].astype("Int64")
    df["chr"] = df["chr"].apply(chr_to_int64)

    if df.empty:
        print("At least one valid region is required to design primers")
        raise SystemExit
    errors = validate_primer_csv(df)
    if not errors:
        chrom = df["chr"].to_list()
        variant = df["variant"].to_list()
        pos_start = df["pos_start"].to_list()
        pos_end = df["pos_end"].to_list()
        build = df["build"].to_list()
        tag = df["tag"].to_list()
        opt_tm = df["opt_tm"].to_list()
        min_tm = df["min_tm"].to_list()
        max_tm = df["max_tm"].to_list()
        opt_gc = df["opt_gc"].to_list()
        min_gc = df["min_gc"].to_list()
        max_gc = df["max_gc"].to_list()
        min_p_size = df["min_product_size"].to_list()
        max_p_size = df["max_product_size"].to_list()
        opt_primer_size = df["opt_primer_size"].to_list()
        min_primer_size = df["min_primer_size"].to_list()
        max_primer_size = df["max_primer_size"].to_list()
        return {
            "chrom": chrom,
            "variant": variant,
            "pos_start": pos_start,
            "pos_end": pos_end,
            "build": build,
            "tag": tag,
            "opt_tm": opt_tm,
            "min_tm": min_tm,
            "max_tm": max_tm,
            "opt_gc": opt_gc,
            "min_gc": min_gc,
            "max_gc": max_gc,
            "min_p_size": min_p_size,
            "max_p_size": max_p_size,
            "opt_primer_size": opt_primer_size,
            "min_primer_size": min_primer_size,
            "max_primer_size": max_primer_size
        }
    else:
        print(errors)
        raise SystemExit


def get_value(param_value, default):
    """
    Return param_value if it is valid; otherwise return default from config
    """
    try:
        if isinstance(param_value, float) and math.isnan(param_value):
            return default
        elif param_value is None:
            return default
        else:
            return param_value
    except Exception:
        return default


def make_list(x):
    """Convert string or None to empty list"""
    if x is None:
        return []
    if isinstance(x, list):
        return x
    if isinstance(x, str):
        return []
    return [x]


def insert_DB(order_primer, tagged_FW, tagged_RV, tag, username,
              password, db_name, db_host, config):
    with open(config, "r") as file:
        config = json.load(file)
    df_insert = order_primer[["chr", "POS", "variant", "Left_Sequence",
                              "Right_Sequence", "Left_Start", "Left_End",
                              "Right_Start", "Right_End", "Pair_Product_Size",
                              "gene", "GRCh"]]
    df_insert["tagged_FW"] = tagged_FW
    df_insert["tagged_RV"] = tagged_RV
    df_insert["tag"] = tag

    db_cols = ["chr", "position", "variant", "left_primer_seq",
               "right_primer_seq", "left_primer_start", "left_primer_end",
               "right_primer_start", "right_primer_end", "product_size",
               "gene", "grch", "tagged_left", "tagged_right", "tag"]
    df_insert = df_insert.applymap(
        lambda x: x.item() if isinstance(x, (np.integer, np.floating)) else x
    )
    df_insert = df_insert.where(pd.notnull(df_insert), None)
    connection = get_postgres_connection(db_name, username, password, db_host)
    cursor = connection.cursor()
    # Insert each row if not duplicated
    for idx, row in df_insert.iterrows():
        values = tuple(row.values)
        query = f"""
        INSERT INTO {config["db"]["schema_name"]}.{config["db"]["table_name"]} ({", ".join(db_cols)})
        VALUES ({", ".join(["%s"] * len(db_cols))})
        ON CONFLICT DO NOTHING;
        """
        cursor.execute(query, values)
    connection.commit()


def generate_bed(primers, build):
    """
    put generated primers into bed format to plot on IGV
    """
    temp_dir = Path("/app/static/temp")
    temp_dir.mkdir(parents=True, exist_ok=True)
    bed_file = temp_dir / f"primers_{build}.bed"
    primers.to_csv(bed_file, sep="\t", header=False, index=False)
    return str(bed_file)


def vcf_to_bed(bed_file, config, build):
    """
    get common SNP within primers to plot on IGV
    """
    with open(config, "r") as file:
        config = json.load(file)
    vcf_dir = Path("/app/static/temp")
    vcf_dir.mkdir(parents=True, exist_ok=True)

    if int(build) == 19:
        vcf_in = config["ref_b37"]["snp_ref"]
        bed_file = config["inter_file"]["bed_file_37"]
        intermediate_vcf = config["inter_file"]["inter_vcf_37"]
        vcf_bed = config["inter_file"]["vcf_bed_37"]
    else:
        vcf_in = config["ref_b38"]["snp_ref"]
        bed_file = config["inter_file"]["bed_file_38"]
        intermediate_vcf = config["inter_file"]["inter_vcf_38"]
        vcf_bed = config["inter_file"]["vcf_bed_38"]
    print("inter", build, intermediate_vcf, vcf_in, bed_file)

    # filter variant with bed file
    subprocess.run([
            "bcftools", "view",
            "-R", bed_file,
            "-o", intermediate_vcf,
            vcf_in
            ], check=True)
    with open(intermediate_vcf, "r") as vcf, open(vcf_bed, "w") as bed:
        for line in vcf:
            if line.startswith("#"):
                continue  # skip header lines
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            vid = fields[2] if fields[2] != "." else "NA"
            ref = fields[3]
            alt = fields[4]
            name = f"{vid}_{ref}_{alt}"
            start = pos - 1  # BED is 0-based
            end = start + len(ref)  # end position
            bed.write(f"{chrom}\t{start}\t{end}\t{name}\n")


def generate_filtered_vcf(bed_file, config, build):

    with open(config, "r") as file:
        config = json.load(file)
    vcf_dir = Path("/app/static/temp")
    vcf_dir.mkdir(parents=True, exist_ok=True)

    if int(build) == 19:
        vcf_in = config["ref_b37"]["snp_ref"]
        bed_file = config["inter_file"]["bed_file_37"]
        vcf_out = config["inter_file"]["vcf_file_37"]
        intermediate_vcf = config["inter_file"]["inter_vcf_37"]
    else:
        vcf_in = config["ref_b38"]["snp_ref"]
        bed_file = config["inter_file"]["bed_file_38"]
        vcf_out = config["inter_file"]["vcf_file_38"]
        intermediate_vcf = config["inter_file"]["inter_vcf_38"]
    # filter variant with bed file
    subprocess.run([
        "bcftools", "view",
        "-R", bed_file,
        "-o", intermediate_vcf,
        vcf_in
    ], check=True)
    # remove duplicated variants
    try:
        subprocess.run([
            "bcftools", "norm",
            "-d", "both",
            "-Oz",
            "-o", vcf_out,
            intermediate_vcf
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print("Exit code:", e.returncode)
        print("stdout:", e.stdout)
        print("stderr:", e.stderr)
    # generate index file
    subprocess.run([
        "bcftools", "index", "--tbi",
        vcf_out
    ], check=True)


def del_file(files_to_remove=None, base_dir="/app"):
    """
    del intermediate files
    """
    if files_to_remove is None:
        files_to_remove = ["designed_primer.fa", "designed_primer.fa.sam"]

    for root, dirs, files in os.walk(base_dir):
        for f in files_to_remove:
            path = os.path.join(root, f)
            if os.path.exists(path):
                os.remove(path)
                print("Deleted:", path)
            else:
                print("None to delete")


def get_postgres_connection(db_name, db_user, db_password, db_host):
    """connect to PostgreSQL"""
    return psycopg2.connect(
        dbname=db_name,
        user=db_user,
        password=db_password,
        host=db_host,
    )


def prepare_df(df, side):
    """
    prepare df for bed file for IGV
    side: "Left" or "Right"
    """
    start_col = f"{side}_Start"
    end_col = f"{side}_End"
    cols = ["chr", start_col, end_col, "Primer_Pair", "GRCh"]

    df_side = df[cols].copy()
    df_side["UID"] = (
        df_side['chr'].astype(str) + "_" +
        df_side[start_col].astype(str) + "_" +
        df_side[end_col].astype(str) + f"_{side.upper()}_" +
        df_side['Primer_Pair'].astype(str)
    )

    df_side = df_side.rename(columns={
        start_col: "start",
        end_col: "end"
    })
    return df_side
