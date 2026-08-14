import pandas as pd
from primer_design.helper_function import get_postgres_connection


def insert_DB(order_primer, tagged_FW, tagged_RV, tag, username,
              password, db_name, db_host, notes="NA", passval="Not_Done",
              archive="No"):

    df_insert = order_primer[[
        "chr", "start_POS", "end_POS", "primer_name", "Left_Sequence",
        "Right_Sequence", "Left_Start", "Left_End",
        "Right_Start", "Right_End", "Pair_Product_Size",
        "gene", "GRCh", "designer"
    ]].copy()

    df_insert["tagged_FW"] = tagged_FW
    df_insert["tagged_RV"] = tagged_RV
    df_insert["tag"] = tag
    df_insert["notes"] = notes
    df_insert["passval"] = passval
    df_insert["archive"] = archive

    # clean empty values
    for col in ["tagged_FW", "tagged_RV", "tag", "notes"]:
        df_insert[col] = df_insert[col].replace([None, ""], "NA")

    for col in ["archive"]:
        df_insert[col] = df_insert[col].replace([None, ""], "No")

    for col in ["passval"]:
        df_insert[col] = df_insert[col].replace([None, ""], "Not_Done")

    # rename to DB columns
    df_insert.rename(columns={
        "start_POS": "start_pos",
        "end_POS": "end_pos",
        "Left_Sequence": "left_primer_seq",
        "Right_Sequence": "right_primer_seq",
        "Left_Start": "left_primer_start",
        "Left_End": "left_primer_end",
        "Right_Start": "right_primer_start",
        "Right_End": "right_primer_end",
        "Pair_Product_Size": "product_size",
        "GRCh": "grch",
        "tagged_FW": "p_tagged_left",
        "tagged_RV": "p_tagged_right",
        "tag": "p_tag",
        "notes": "p_notes",
        "passval": "p_passed_validation",
        "archive": "p_archive",
        "designer": "p_designer"
    }, inplace=True)

    df_insert = df_insert.where(pd.notnull(df_insert), None)

    connection = get_postgres_connection(db_name, username, password, db_host)
    cursor = connection.cursor()
    query = """
    SELECT * FROM primer_tool.insert_primer_with_batch(
        %s::varchar,
        %s::int,
        %s::int,
        %s::int,
        %s::varchar,
        %s::varchar,
        %s::varchar,
        %s::int,
        %s::int,
        %s::int,
        %s::int,
        %s::int,
        %s::varchar,

        %s::varchar,
        %s::varchar,
        %s::varchar,
        %s::varchar,
        %s::text,
        %s::text,
        %s::varchar,
        %s::varchar,
        %s::varchar,
        %s::varchar,
        %s::varchar,
        %s::varchar,
        %s::varchar,
        %s::varchar
    );
    """
    inserted_ids = []
    for _, row in df_insert.iterrows():

        values = (
            str(row["chr"]),
            int(row["start_pos"]),
            int(row["end_pos"]),
            int(row["grch"]),
            str(row["primer_name"]) if row["primer_name"] else "NA",
            str(row["left_primer_seq"]),
            str(row["right_primer_seq"]),
            int(row["left_primer_start"]) if row["left_primer_start"] is not None else 0,
            int(row["left_primer_end"]) if row["left_primer_end"] is not None else 0,
            int(row["right_primer_start"]) if row["right_primer_start"] is not None else 0,
            int(row["right_primer_end"]) if row["right_primer_end"] is not None else 0,
            int(row["product_size"]) if row["product_size"] is not None else 0,
            str(row["gene"]),

            str(row["p_tagged_left"]),
            str(row["p_tagged_right"]),
            str(row["p_tag"]),
            str(row["p_notes"]),
            str(row["p_passed_validation"]),
            str(row["p_archive"]),
            None,  # p_mix
            None,  # p_dilution_date
            None,  # p_tray
            None,  # p_freezer
            None,  # p_grid_fw
            None,   # p_grid_rv
            None,  # p_manufacturer
            str(row["p_designer"])
        )

        cursor.execute(query, values)
        upi, primer_id = cursor.fetchone()
        inserted_ids.append([upi, primer_id])

    connection.commit()
    cursor.close()
    connection.close()
    return inserted_ids