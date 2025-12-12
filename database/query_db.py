import argparse
import json
from psycopg2 import sql
from primer_design.helper_function import get_postgres_connection


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
                    schema, table, col1, col2, col3,
                    value1, value2, value3):
    """
    search query with max three filter values
    """
    # Combine input values
    value = [value1, value2, value3]
    col = [col1, col2, col3]

    # Filter out empty or None
    value_list = [x for x in value if x not in (None, '')]
    col_list = [x for x in col if x not in (None, '')]

    # Build WHERE conditions
    combined = sql.SQL(" AND ").join(
        sql.SQL("{} = %s").format(sql.Identifier(c)) for c in col_list
    )

    query = sql.SQL("SELECT * FROM {} WHERE ").format(sql.Identifier(schema, table)) + combined

    # Connect to db
    conn = get_postgres_connection(dbname, user, password, host)
    cursor = conn.cursor()
    cursor.execute(query, value_list)

    # Fetch results
    rows = cursor.fetchall()
    col_names = [desc[0] for desc in cursor.description]
    result = [dict(zip(col_names, row)) for row in rows]

    conn.close()
    print(result)
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
