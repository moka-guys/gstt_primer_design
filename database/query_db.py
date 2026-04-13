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
                    schema, table, col1, col2, col3, col4,
                    value1, value2, value3, value4, col5, value7,
                    pos_start=None, pos_end=None):
    """
    Search query with exact filters + optional overlapping range
    """
    # Pair columns and values
    col_value_pairs = [
        (str(c).strip(), v)
        for c, v in zip([col1, col2, col3, col4, col5],
                        [value1, value2, value3, value4, value7])
        if c not in (None, '') and v not in (None, '')
    ]

    conditions = []
    value_list = []

    for c, v in col_value_pairs:
        if c == col5:  # grch input
            if str(v) == "37":
                conditions.append(sql.SQL("{} IN (%s, %s)").format(sql.Identifier(c)))
                value_list.extend(["37", "38"])
            else:
                conditions.append(sql.SQL("{} = %s").format(sql.Identifier(c)))
                value_list.append(v)
        else:
            conditions.append(sql.SQL("{} = %s").format(sql.Identifier(c)))
            value_list.append(v)

    # Only add range condition if pos_start/pos_end are valid integers
    try:
        pos_start_int = int(pos_start) if pos_start not in (None, '') else None
    except ValueError:
        pos_start_int = None

    try:
        pos_end_int = int(pos_end) if pos_end not in (None, '') else None
    except ValueError:
        pos_end_int = None

    if pos_start_int is not None and pos_end_int is not None:
        if str(value7) == "38":
            conditions.append(
                sql.SQL("left_primer_start >= %s AND right_primer_end <= %s")
            )
            value_list.extend([pos_start_int, pos_end_int])
        elif str(value7) == "37":
            # 37 uses original, 38 uses liftover
            lifted_start = liftover_37to38(f"chr{value1}", pos_start_int)
            lifted_end = liftover_37to38(f"chr{value1}", pos_end_int)
            conditions.append(sql.SQL("""
            (
                (grch = %s AND left_primer_start >= %s AND right_primer_end <= %s)
                OR
                (grch = %s AND left_primer_start >= %s AND right_primer_end <= %s)
            )
            """))

            value_list.extend([
                "37", pos_start_int, pos_end_int,
                "38", lifted_start, lifted_end
            ])

    # query
    base_query = sql.SQL("SELECT * FROM {}").format(sql.Identifier(schema, table))

    if conditions:
        query = base_query + sql.SQL(" WHERE ") + sql.SQL(" AND ").join(conditions)
    else:
        query = base_query

    query += sql.SQL(" LIMIT 100")
    conn = get_postgres_connection(dbname, user, password, host)
    cursor = conn.cursor()
    print("QUERY:", query.as_string(conn))
    print("VALUES:", value_list)

    cursor.execute(query, value_list)
    rows = cursor.fetchall()
    col_names = [desc[0] for desc in cursor.description]
    result = [dict(zip(col_names, row)) for row in rows]
    # reorder to show passedvalidation in earlier col
    col_names.insert(5, col_names.pop(col_names.index("passedvalidation")))
    result = [
                {col: row[col] for col in col_names}
                for row in result
            ]
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
