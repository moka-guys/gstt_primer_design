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
            int(row["start_POS"]),
            int(row["end_POS"]),
            int(row["GRCh"]),
            str(row["primer_name"]) if row["primer_name"] else "NA",
            str(row["Left_Sequence"]),
            str(row["Right_Sequence"]),
            int(row["Left_Start"]) if row["Left_Start"] is not None else 0,
            int(row["Left_End"]) if row["Left_End"] is not None else 0,
            int(row["Right_Start"]) if row["Right_Start"] is not None else 0,
            int(row["Right_End"]) if row["Right_End"] is not None else 0,
            int(row["Pair_Product_Size"]) if row["Pair_Product_Size"] is not None else 0,
            str(row["gene"]),

            str(row["tagged_FW"]),
            str(row["tagged_RV"]),
            str(row["tag"]),
            str(row["notes"]),
            str(row["passval"]),
            str(row["archive"]),
            None,  # p_mix
            None,  # p_dilution_date
            None,  # p_tray
            None,  # p_freezer
            None,  # p_grid_fw
            None,   # p_grid_rv
            None,  # p_manufacturer
            str(row["designer"])
        )

        cursor.execute(query, values)
        upi, primer_id = cursor.fetchone()
        inserted_ids.append([upi, primer_id])

    connection.commit()
    cursor.close()
    connection.close()
    return inserted_ids
