import argparse
import json
from psycopg2 import sql
from primer_design.helper_function import get_postgres_connection, liftover_37to38


def get_arguments() -> argparse.Namespace:
    """
    Uses argparse to define and handle command line input arguments
    and help menu
        Return argparse.Namespace (object): Contains the parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t",
        "--table",
        required=True,
        help="table name",
        )
    parser.add_argument(
        "-db",
        "--database",
        required=True,
        help="database name",
        )
    parser.add_argument(
        "-c1",
        "--query_col1",
        help="col1 to search",
        )
    parser.add_argument(
        "-c2",
        "--query_col2",
        help="col2 to search",
        )
    parser.add_argument(
        "-c3",
        "--query_col3",
        help="col3 to search",
        )
    parser.add_argument(
        "-v1",
        "--query_value1",
        help="value1 to search",
        )
    parser.add_argument(
        "-v2",
        "--query_value2",
        help="value2 to search",
        )
    parser.add_argument(
        "-v3",
        "--query_value3",
        help="value3 to search",
        )
    parser.add_argument(
        "--filter_multiple",
        action="store_true",
        default=False,
        help="add this args to query with multiple filters",
        )
    parser.add_argument(
        "-f",
        "--filters",
        type=argparse.FileType('r'),
        default="./multiple_filters.json",
        help="json file path",
        )
    parser.add_argument(
        "-u",
        "--user",
        help="username",
        )
    parser.add_argument(
        "-pw",
        "--password",
        help="username",
        )
    return parser.parse_args()


def search_postgres(dbname, user, password, host,
                    chr_val, gene_val, variant_val,
                    validation_val, grch_val,
                    pos_start=None, pos_end=None):

    conditions = []
    values = []

    base_query = sql.SQL("""
        SELECT
            p.*,
            b.*
        FROM primer_tool.primers p
        LEFT JOIN primer_tool.primer_batches b
        ON p.primer_id = b.primer_id
    """)

    if chr_val:
        conditions.append(sql.SQL("p.chr = %s"))
        values.append(chr_val)

    if gene_val:
        conditions.append(sql.SQL("p.gene = %s"))
        values.append(gene_val)

    if variant_val:
        conditions.append(sql.SQL("p.variant = %s"))
        values.append(variant_val)

    if validation_val:
        conditions.append(sql.SQL("b.passed_validation = %s"))
        values.append(validation_val)

    try:
        pos_start = int(pos_start) if pos_start else None
        pos_end = int(pos_end) if pos_end else None
    except:
        pos_start = None
        pos_end = None

    def safe_liftover(chrom, pos):
        try:
            lifted = liftover_37to38(chrom, pos)
            return int(lifted) if lifted else None
        except:
            return None
    if grch_val:
        grch_val = str(grch_val)

        if grch_val == "38":
            conditions.append(sql.SQL("p.grch = %s"))
            values.append(38)

            if pos_start is not None and pos_end is not None:
                conditions.append(sql.SQL("""
                    p.left_primer_start >= %s
                    AND p.right_primer_end <= %s
                """))
                values.extend([pos_start, pos_end])

        # liftover if grch37 and pos is given
        elif grch_val == "37":

            chrom = f"chr{chr_val}" if chr_val else None

            if pos_start is not None and pos_end is not None and chrom:

                lifted_start = safe_liftover(chrom, pos_start)
                lifted_end = safe_liftover(chrom, pos_end)

                # If liftover works → dual query
                if lifted_start and lifted_end:
                    conditions.append(sql.SQL("""
                        (
                            (p.grch = %s AND p.left_primer_start >= %s AND p.right_primer_end <= %s)
                            OR
                            (p.grch = %s AND p.left_primer_start >= %s AND p.right_primer_end <= %s)
                        )
                    """))

                    values.extend([
                        37, pos_start, pos_end,
                        38, lifted_start, lifted_end
                    ])

                else:
                    # fallback if liftover fails
                    conditions.append(sql.SQL("""
                        p.grch = %s
                        AND p.left_primer_start >= %s
                        AND p.right_primer_end <= %s
                    """))
                    values.extend([37, pos_start, pos_end])

            else:
                conditions.append(sql.SQL("p.grch IN (37)"))

    query = base_query

    if conditions:
        query += sql.SQL(" WHERE ") + sql.SQL(" AND ").join(conditions)
    query += sql.SQL(" ORDER BY insert_time ASC, batch_id ASC")
    query += sql.SQL(" LIMIT 100")
    conn = get_postgres_connection(dbname, user, password, host)
    cursor = conn.cursor()
    cursor.execute(query, values)

    rows = cursor.fetchall()
    cols = [d[0] for d in cursor.description]
    result = [dict(zip(cols, row)) for row in rows]
    if "passed_validation" in cols:
        cols.insert(5, cols.pop(cols.index("passed_validation")))
        result = [{col: row[col] for col in cols} for row in result]
    conn.close()
    return result


def filter_multiple_postgres(dbname, user, password, schema, table, filters):
    """
    search query for if filter values are given in json file
    """
    # connect to db
    conn = get_postgres_connection(dbname, user, password)
    cursor = conn.cursor()

    # build condition
    conditions = [sql.SQL("{} = %s").format(sql.Identifier(col)) for col in filters.keys()]
    where_clause = sql.SQL(" AND ").join(conditions)
    params = list(filters.values())

    # build query
    query = sql.SQL("SELECT * FROM {} WHERE {}").format(
        sql.Identifier(schema, table),
        where_clause
    )
    cursor.execute(query, params)
    rows = cursor.fetchall()
    # Get column names
    col_names = [desc[0] for desc in cursor.description]
    # Convert to list of dicts
    result = [dict(zip(col_names, row)) for row in rows]

    print(result)
    conn.close()
    return result


if __name__ == "__main__":
    parsed_args = get_arguments()
    if not parsed_args.filter_multiple:
        search_postgres(parsed_args.database, parsed_args.user, parsed_args.password,
                        parsed_args.table, parsed_args.query_col1,
                        parsed_args.query_col2, parsed_args.query_col3, parsed_args.query_value1,
                        parsed_args.query_value2, parsed_args.query_value3)
    else:
        filters = json.load(parsed_args.filters)
        filter_multiple_postgres(parsed_args.database, parsed_args.user, parsed_args.password,
                                 parsed_args.table, filters)
