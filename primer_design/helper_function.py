import os
from typing import Any
import logging
import json
import pandas as pd
import math
import subprocess
from pydantic import BaseModel, ValidationError, field_validator, model_validator
import psycopg2
from psycopg2.extensions import connection
from pathlib import Path
from pyliftover import LiftOver
import tempfile
import shutil
import gzip
from pathlib import Path

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(BASE_DIR, "config.json")
with open(config_path, "r") as file:
    config = json.load(file)


class PrimerRecord(BaseModel):
    """Validate primer design input parameters.
    This model performs input validation for primer design
    requests submitted through the application.
    """
    chr: str | int
    primer_name: str
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
    @classmethod
    def non_negative(cls, v, info) -> int | float:
        """Validate that numeric values are non-negative.
        Args:
            v: Value being validated.
            info: Pydantic validation metadata containing the field name.

        Returns:
            The validated value.

        Raises:
            ValueError: If the value is less than zero.
        """
        if v < 0:
            raise ValueError(f"{info.field_name} must not be negative")
        return v

    @field_validator("chr")
    def validate_chr(cls, v) -> str:
        """Validate chromosome value.
        Args:
            v: Chromosome identifier.

        Returns:
            Chromosome as a normalized string.

        Raises:
            ValueError: If chromosome is not 1-22, X, or Y.
        """
        allowed = {str(i) for i in range(1, 23)} | {"X", "Y"}
        v_str = str(v).strip()
        if v_str not in allowed:
            raise ValueError("chr must be 1-22, X, or Y")
        return v_str

    @field_validator("build")
    def validate_build(cls, v) -> int:
        """Validate genome build.
        Args:
            v: Genome build number.

        Returns:
            The validated genome build.

        Raises:
            ValueError: If build is not 37 or 38.
        """
        if v not in (37, 38):
            raise ValueError("build must be 37 or 38")
        return v

    @field_validator("tag")
    def validate_tag(cls, v) -> str:
        """Validate primer tag.
        Args:
            v: Primer tag name.

        Returns:
            The validated tag.

        Raises:
            ValueError: If an unsupported tag is provided.
        """
        if v not in ("M13", "T1", "T2", "T3", "T4", "sT1", "FAM",
                     "VIC", "NED", "PET", "ATTO550", "CY5",
                     "ABY", "HEX", "no_tag"):
            raise ValueError("invalid tag")
        return v

    @model_validator(mode="after")
    def check_opt_tm(cls, model) -> "PrimerRecord":
        """Validate relationships between primer design parameters.
        Checks that optimal values fall within configured ranges and
        verifies that position and product-size ranges are valid.

        Args:
            model: Primer design parameters after field validation.

        Returns:
            The validated model.

        Raises:
            ValueError: If any parameter relationship is invalid.
        """
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
    """Create and configure a logger for primer design jobs.
    Creates a file-based logger that writes primer design messages to
    a timestamped log file under the primer_log directory.
    Args:
        file_dir: Base directory for log output.
        datetimestr: Timestamp string used in the log file name and
        logger name.

    Returns:
        logging.Logger: Configured logger instance.
    """

    logger = logging.getLogger(f"primer_log_{datetimestr}")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    if not logger.hasHandlers():
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )

        outdir = os.path.join(file_dir, "primer_log")
        os.makedirs(outdir, exist_ok=True)

        log_file = os.path.join(
            outdir,
            f"primer_design_{datetimestr}.log"
        )

        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_tag(df, tag) -> tuple[str, str, str, str, str, str]:
    """Generate tagged primer sequences for ordering.
    Args:
        df: DataFrame containing primer sequences.
        tag: Tag configuration name.

    Returns:
        tuple[str, str, str, str, str, str]:
            Forward primer sequence,
            reverse primer sequence,
            tagged forward primer,
            tagged reverse primer,
            forward tag name,
            reverse tag name.
    """
    FW_primer = df["Left_Sequence"][0]
    RV_primer = df["Right_Sequence"][0]
    FW_tag = config[tag]["F"]
    RV_tag = config[tag]["R"]
    tagged_FW = FW_tag + FW_primer
    tagged_RV = RV_tag + RV_primer
    tag_name_FW = config[tag]["F_name"]
    tag_name_RV = config[tag]["R_name"]

    return FW_primer, RV_primer, tagged_FW, tagged_RV, tag_name_FW, tag_name_RV


def chr_to_int64(x) -> Any:
    """Convert chromosome value to pandas Int64 when possible.
    Args:
        x: Chromosome value.

    Returns:
        pandas.Int64Dtype-compatible integer if conversion succeeds,
        otherwise the original value.
    """
    try:
        return pd.Int64Dtype().type(int(x))
    except ValueError:
        return x


def validate_primer_csv(df: pd.DataFrame) -> list[dict[str, Any]]:
    """Validate primer design records from a DataFrame.
    Each row is validated against the PrimerRecord model.
    Args:
        df: Input DataFrame containing primer design parameters.

    Returns:
        list[dict]: Validation errors. Empty list indicates all records
        are valid.
    """
    errors = []

    for idx, row in df.iterrows():
        try:
            record = PrimerRecord(**row.to_dict())

        except ValidationError as e:
            errors.append({"row": idx, "errors": e.errors()})

    return errors


def parse_csv(file) -> dict[str, Any]:
    """Parse and validate a primer design input CSV file.
    Reads the CSV file, validates all records, and returns primer
    design parameters as lists suitable for batch processing.
    Args:
        file: Path to the CSV file.

    Returns:
        dict: Parsed parameter values on success, or an error
        dictionary containing validation failures.
    """
    df = pd.read_csv(file)
    # drop empty rows
    #cols = ["chr", "pos_start", "pos_end", "build", "tag"]
    #df = df[df[cols].notna().all(axis=1)]
    for c in ["pos_start", "pos_end", "build"]:
        df[c] = df[c].astype("Int64")
    df["chr"] = df["chr"].apply(chr_to_int64)

    if df.empty:
        print("At least one valid region is required to design primers")
        return {
                "error": "At least one valid region is required to design primers"
                }
    errors = validate_primer_csv(df)
    if not errors:
        chrom = df["chr"].to_list()
        primer_name = df["primer_name"].to_list()
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
            "primer_name": primer_name,
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
        return {
                "error": errors
                }


def get_value(param_value, default) -> Any:
    """Return a parameter value or a default value.
    Handles missing values from user input, including None and NaN,
    and returns a configured default value when appropriate.
    Args:
        param_value: User-supplied parameter value.
        default: Default value to use when the parameter is invalid.

    Returns:
        The supplied parameter value if valid; otherwise the default
        value.
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


def make_list(x) -> list:
    """Convert a value to a list.
    Args:
        x: Value to convert.

    Returns:
        list: Empty list for None or string values, the original list
        if already a list, otherwise a single-item list.
    """
    if x is None:
        return []
    if isinstance(x, list):
        return x
    if isinstance(x, str):
        return []
    return [x]


def generate_bed(primers, build, job_id) -> str:
    """Create a BED file from primer coordinates.
    Args:
        primers: DataFrame containing primer coordinates.
        build: Genome build identifier.
        job_id: Unique job identifier.

    Returns:
        str: Generated BED file name.
    """
    temp_dir = Path(config["directory"]["temp_folder"])
    temp_dir.mkdir(parents=True, exist_ok=True)
    bed_file = temp_dir / f"primers_{build}_{job_id}.bed"
    primers.to_csv(bed_file, sep="\t", header=False, index=False)
    return str(os.path.basename(bed_file))


def vcf_to_bed(bed_file, build, job_id) -> tuple[str, str]:
    """Extract common variants overlapping primer regions.
    Filters variants from both exome and genome reference VCFs using
    primer regions, combines the results, removes exact duplicate
    variants, and converts the matching records to BED format.

    Args:
        bed_file: BED file containing primer coordinates.
        build: Genome build identifier.
        job_id: Unique job identifier.

    Returns:
        tuple[str, str]:
        Generated BED-format variant file name and filtered VCF
        file name.
    """
    vcf_dir = Path(config["directory"]["temp_folder"])
    vcf_dir.mkdir(parents=True, exist_ok=True)

    bed_file = f"{vcf_dir}/{bed_file}"

    if int(build) == 19:
        exome_vcf = config["ref_b37"]["snp_ref_exome"]
        genome_vcf = config["ref_b37"]["snp_ref_genome"]

        intermediate_exome = (
            f"{vcf_dir}/tmp_filtered_exome_19_{job_id}.vcf.gz"
        )
        intermediate_genome = (
            f"{vcf_dir}/tmp_filtered_genome_19_{job_id}.vcf.gz"
        )
        intermediate_vcf = (
            f"{vcf_dir}/tmp_filtered_19_{job_id}.vcf.gz"
        )
        vcf_bed = f"{vcf_dir}/vcf_bed_19_{job_id}.bed"

    else:
        exome_vcf = config["ref_b38"]["snp_ref_exome"]
        genome_vcf = config["ref_b38"]["snp_ref_genome"]

        intermediate_exome = (
            f"{vcf_dir}/tmp_filtered_exome_38_{job_id}.vcf.gz"
        )
        intermediate_genome = (
            f"{vcf_dir}/tmp_filtered_genome_38_{job_id}.vcf.gz"
        )
        intermediate_vcf = (
            f"{vcf_dir}/tmp_filtered_38_{job_id}.vcf.gz"
        )
        vcf_bed = f"{vcf_dir}/vcf_bed_38_{job_id}.bed"

    # Filter exome variants
    subprocess.run([
        "bcftools", "view",
        "-R", bed_file,
        "-i", "INFO/AF > 0.01",
        "-O", "z",
        "-o", intermediate_exome,
        exome_vcf
    ], check=True)

    # Filter genome variants
    subprocess.run([
        "bcftools", "view",
        "-R", bed_file,
        "-i", "INFO/AF > 0.01",
        "-O", "z",
        "-o", intermediate_genome,
        genome_vcf
    ], check=True)

    # Combine exome + genome and remove exact duplicates
    seen = set()
    header_written = False

    with gzip.open(intermediate_vcf, "wt") as output_vcf:

        for input_vcf in [intermediate_genome, intermediate_exome]:

            with gzip.open(input_vcf, "rt") as vcf:

                for line in vcf:

                    if line.startswith("#"):
                        if not header_written:
                            output_vcf.write(line)
                        continue

                    fields = line.rstrip("\n").split("\t")

                    chrom = fields[0]
                    pos = fields[1]
                    ref = fields[3]
                    alt = fields[4]

                    # Unique variant based on CHROM + POS + REF + ALT
                    variant_key = (
                        chrom,
                        pos,
                        ref,
                        alt
                    )

                    if variant_key in seen:
                        continue

                    seen.add(variant_key)

                    output_vcf.write(line)

            header_written = True

    # Convert combined VCF to BED
    with gzip.open(intermediate_vcf, "rt") as vcf, \
            open(vcf_bed, "w") as bed:

        for line in vcf:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[0]
            pos = int(fields[1])
            vid = fields[2] if fields[2] != "." else "NA"
            ref = fields[3]
            alt = fields[4]
            info_field = fields[7]

            info_dict = {
                item.split("=", 1)[0]: item.split("=", 1)[1]
                for item in info_field.split(";")
                if "=" in item
            }

            af = info_dict.get("AF")

            if af is not None:
                af_value = float(af)
            else:
                af_value = "NA"

            name = f"{vid}_{ref}_{alt}_{af_value}"

            # BED is 0-based
            start = pos - 1
            end = start + len(ref)

            bed.write(
                f"{chrom}\t{start}\t{end}\t{name}\n"
            )

    return (
        str(os.path.basename(vcf_bed)),
        str(os.path.basename(intermediate_vcf))
    )


def del_file(files_to_remove=None, base_dir="/app") -> None:
    """Delete intermediate files from a directory tree.
    Recursively searches a directory and removes matching files.
    Args:
        files_to_remove: List of file names to delete. Defaults to
        common primer-design intermediate files.
        base_dir: Root directory to search.

    Returns:
        None.
    """
    if files_to_remove is None:
        files_to_remove = ["designed_primer.fa", "designed_primer.fa.sam"]

    for root, dirs, files in os.walk(base_dir):
        for f in files_to_remove:
            path = os.path.join(root, f)
            if os.path.exists(path):
                os.remove(path)
            #     print("Deleted:", path)
            # else:
            #     print("None to delete")


def get_postgres_connection(db_name, db_user,
                            db_password, db_host) -> connection:
    """Create a PostgreSQL database connection.
    Args:
        db_name: Database name.
        db_user: Database username.
        db_password: Database password.
        db_host: Database host.

    Returns:
        psycopg2.extensions.connection: PostgreSQL connection object.
    """
    return psycopg2.connect(
        dbname=db_name,
        user=db_user,
        password=db_password,
        host=db_host,
    )


def prepare_df(df, side) -> pd.DataFrame:
    """Prepare primer coordinates for IGV BED file generation.
    Extracts primer coordinates for a specified primer side and
    creates a unique identifier for each primer.
    Args:
        df: DataFrame containing primer design results.
        side: Primer side to process ("Left" or "Right").

    Returns:
        pd.DataFrame: Reformatted DataFrame containing chromosome,
        start position, end position, genome build, primer pair,
        and a unique identifier.
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


def prepare_order_sheet(df) -> tuple[str, str, str, str, str, str]:
    """Prepare tagged primer information for ordering.
    Args:
        df: DataFrame containing primer information and order tag.

    Returns:
        tuple[str, str, str, str, str, str]:
            Forward primer,
            reverse primer,
            tagged forward primer,
            tagged reverse primer,
            forward tag name,
            reverse tag name.
    """
    (FW_primer, RV_primer,
     tagged_FW, tagged_RV,
     tag_name_FW, tag_name_RV) = get_tag(df, df["order_tag"][0])

    return FW_primer, RV_primer, tagged_FW, tagged_RV, tag_name_FW, tag_name_RV


def liftover(chrom, pos, build) -> int | None:
    """Convert a genomic coordinate between genome builds.
    Uses pyliftover to map a coordinate between GRCh37 and GRCh38.

    Args:
        chrom: Chromosome name.
        pos: Genomic position.
        build: Source genome build (37 or 38).

    Returns:
        int | None: Converted position if mapping succeeds,
        otherwise None.
    """

    try:
        build = int(build)

        if build == 37:
            liftover_ref = config["ref_b37"]["liftover"]
        elif build == 38:
            liftover_ref = config["ref_b38"]["liftover"]
        else:
            return None

        lo = LiftOver(liftover_ref)
        result = lo.convert_coordinate(chrom, int(pos))

        if not result:
            return None
        new_pos = result[0][1]

        return int(new_pos)

    except Exception as e:
        print(f"liftover error: {e}")
        return None


def liftover_crossmap(chrom, start, end, grch) -> int | None:
    """Convert genomic coordinates using CrossMap.
    Creates a temporary BED file, performs coordinate conversion
    using CrossMap, and returns the converted end position.

    Args:
        chrom: Chromosome name.
        start: Start coordinate.
        end: End coordinate.
        grch: Source genome build.

    Returns:
        int | None: Converted coordinate if successful, otherwise None.
    """

    crossmap_path = shutil.which("CrossMap")
    if int(grch) == 37:
        chain_file = config["ref_b37"]["crossmap_ref"]
    elif int(grch) == 38:
        chain_file = config["ref_b38"]["crossmap_ref"]
    # create temporary input/output files

    with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=config["directory"]["output_folder"]) as in_bed:
        in_bed.write(f"{chrom}\t{start}\t{end}\n")
        input_path = in_bed.name

    output_path = input_path + ".out.bed"
    try:
        # run CrossMap
        cmd = [
            crossmap_path,
            "bed",
            chain_file,
            input_path,
            output_path,
        ]

        subprocess.run(cmd, check=True)

        # read result
        with open(output_path) as f:
            line = f.readline().strip()

        if not line:
            return None

        fields = line.split("\t")

        return int(fields[2])

    finally:
        # cleanup temp files
        for path in [
            input_path,
            output_path,
            output_path + ".unmap",
        ]:
            if os.path.exists(path):
                os.remove(path)
