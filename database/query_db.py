from datetime import datetime, timedelta
from typing import Any, Optional
from psycopg2 import sql
from primer_design.helper_function import (get_postgres_connection, liftover,
                                           liftover_crossmap)


def safe_liftover(chrom, pos, build) -> tuple[Optional[int], str]:
    """
    Function to liftover between GRCh37 and GRCh38.
    Two different python packages are used to do liftover and the results are
    compared and only the consistent results are taken
    Args:
        chrom: chromosome number
        pos: genomic position
        build: GRCh

    Returns:
        Tuple of liftover value (Optional) and message string
    """
    try:
        print(chrom, pos, build)
        lifted = liftover(chrom, pos, build)
        print(lifted)
        crossmap = liftover_crossmap(chrom, pos, pos, build)
        print(crossmap)

        if lifted is None or crossmap is None:
            msg_to_return = "No corresponding lift over region is found"
            return None, msg_to_return

        if int(lifted) == int(crossmap):
            msg_to_return = ("Lift over is checked by two different "
                             "tools and results are consistent")

            return int(lifted), msg_to_return

        else:
            msg_to_return = ("Lift over is not consistent for this region;"
                             " therefore, region is not lifted over")
            return None, msg_to_return

    except Exception as exc:
        print("error:", exc)
        return None, str(exc)


def validate_position_inputs(variant_pos, pos_start, pos_end, grch_val) -> None:
    """
    Check if required field(s) are provided for query
    Args:
        variant_pos: genomic position of variant
        pos_start: genomic start position
        pos_end: genomic end position
        grch_val: GRCh
    """
    if variant_pos or pos_start or pos_end:
        if grch_val not in ("37", "38"):
            raise ValueError(
                "GRCh must be 37 or 38 when searching by genomic or variant position."
            )


def get_base_query(db_schema) -> sql.SQL:
    """
    Get base SQL query
    Args:
        db_schema: database schema

    Returns:
        SQL query
    """
    base_query = sql.SQL("""
        SELECT
            b.primer_id,
            p.chr,
            p.grch,
            p.primer_name,
            b.passed_validation,
            b.archive,
            b.notes,
            b.tag,
            p.left_primer_seq,
            p.right_primer_seq,
            p.left_primer_start,
            p.left_primer_end,
            p.right_primer_start,
            p.right_primer_end,
            p.product_size,
            p.gene,
            b.mix,
            b.arrival_date,
            b.tray,
            b.freezer,
            b.grid_fw,
            b.grid_rv,
            b.manufacturer,
            b.designer,
            b.insert_time,
            p.upi,

            -- LEFT highlight
            (
                %s IS NOT NULL AND %s IS NOT NULL AND
                (
                    (p.grch = %s AND p.left_primer_start >= %s AND p.left_primer_end <= %s)
                    OR
                    (p.grch = %s AND p.left_primer_start >= %s AND p.left_primer_end <= %s)
                )
            ) AS lt_in_range,

            -- RIGHT highlight
            (
                %s IS NOT NULL AND %s IS NOT NULL AND
                (
                    (p.grch = %s AND p.right_primer_start >= %s AND p.right_primer_end <= %s)
                    OR
                    (p.grch = %s AND p.right_primer_start >= %s AND p.right_primer_end <= %s)
                )
            ) AS rt_in_range

        FROM {}.primers p
        LEFT JOIN {}.primer_batches b
            ON p.upi = b.upi
    """).format(
        sql.Identifier(db_schema),
        sql.Identifier(db_schema)
    )
    return base_query


def get_condition_value(chr_val, gene_val, validation_val, archive_val,
                        notes, primer_name, primer_id, tray, start_date,
                        end_date) -> tuple[list[sql.SQL], list[Any]]:
    """
    Get the filter conditions and values for primer query
    Args:
        chr_val: chromosome
        gene_val: gene
        validation_val: validation status
        archive_val: archiving status
        notes: notes from user
        primer_name: name of primer
        primer_id: primer ID
        tray: primer tray
        start_date: start date to search
        end_date: end date to search

    Returns:
        Tuple of list of SQL and values if any
    """
    conditions = []
    values = []
    if chr_val:
        conditions.append(sql.SQL("p.chr ILIKE %s"))
        values.append(chr_val)

    if gene_val:
        conditions.append(sql.SQL("p.gene ILIKE %s"))
        values.append(gene_val)

    if validation_val:
        conditions.append(sql.SQL("b.passed_validation = %s"))
        values.append(validation_val)

    if archive_val:
        conditions.append(sql.SQL("b.archive = %s::yes_no"))
        values.append(archive_val)
    if notes:
        conditions.append(sql.SQL("b.notes ILIKE %s"))
        values.append(f"%{notes}%")

    if primer_name:
        conditions.append(sql.SQL("p.primer_name ILIKE %s"))
        values.append(f"%{primer_name}%")

    if primer_id:
        conditions.append(sql.SQL("b.primer_id = %s"))
        values.append(int(primer_id))

    if tray:
        conditions.append(sql.SQL("b.tray ILIKE %s"))
        values.append(tray)

    if start_date and end_date:
        end_date_plus_one = (datetime.strptime(end_date, "%Y-%m-%d") +
                             timedelta(days=1)).strftime("%Y-%m-%d")
        conditions.append(sql.SQL("b.insert_time >= %s AND b.insert_time < %s"))
        values.extend([start_date, end_date_plus_one])
    return conditions, values


def check_liftover_required(pos_start, pos_end, msg_to_return,
                            variant_pos) -> tuple[int | None, int | None, str,
                                                  int | None]:
    """
    Check if liftover is required
    Args:
        pos_start: genomic start position
        pos_end: genomic end position
        msg_to_return: message to put on web app
        variant_pos: genomic position of variant
    Returns:
        genomic positions (if any) and message
    """

    try:
        pos_start = int(pos_start) if pos_start is not None else None
        pos_end = int(pos_end) if pos_end is not None else None
        msg_to_return = "Range lift over done."

    except (ValueError, TypeError):
        msg_to_return = "Range lift over is not done. "
        pos_start = None
        pos_end = None

    try:
        variant_pos = int(variant_pos) if variant_pos is not None else None
        msg_to_return += "Variant lift over done"
    except (ValueError, TypeError):
        msg_to_return += "Variant lift over is not done"
        variant_pos = None
    return pos_start, pos_end, msg_to_return, variant_pos


def get_highlighted_value(pos_start, pos_end, grch_val,
                          liftover_grch, lifted_start, lifted_end) -> list[Any]:
    """
    Get the values to highlight in the output table
    Args:
        pos_start: original genomic start position
        pos_end: original genomic end position
        grch_val: original GRCh
        liftover_grch: liftover GRCh
        lifted_start: liftover genomic start position
        lifted_end: liftover genomic end position
    Returns:
        list of genomic positions and GRCh
    """
    select_values = [
        # LEFT highlight
        pos_start, pos_end,
        grch_val, pos_start, pos_end,
        liftover_grch, lifted_start, lifted_end,

        # RIGHT highlight
        pos_start, pos_end,
        grch_val, pos_start, pos_end,
        liftover_grch, lifted_start, lifted_end,
    ]
    return select_values


def append_pos_liftover(use_liftover, conditions, values, grch_val,
                        pos_start, pos_end, liftover_grch, lifted_start,
                        lifted_end) -> tuple[list[sql.SQL], list[Any]]:
    """
    Append the range liftover conditions and values for query
    Args:
        use_liftover: bool to use liftover or not
        conditions: sql conditions for query
        values: input values for query
        grch_val: original GRCh
        pos_start: original genomic start position
        pos_end: original genomic end position
        liftover_grch: liftover GRCh
        lifted_start: liftover genomic start position
        lifted_end: liftover genomic end position

    Returns:
        Tuple of list of SQL and values if any
    """
    if use_liftover:
        conditions.append(sql.SQL("""
            (
                (
                    p.grch = %s
                    AND (
                        (p.left_primer_start >= %s AND p.left_primer_end <= %s)
                        OR
                        (p.right_primer_start >= %s AND p.right_primer_end <= %s)
                    )
                )
                OR
                (
                    p.grch = %s
                    AND (
                        (p.left_primer_start >= %s AND p.left_primer_end <= %s)
                        OR
                        (p.right_primer_start >= %s AND p.right_primer_end <= %s)
                    )
                )
            )
        """))

        values.extend([
            grch_val,
            pos_start, pos_end,
            pos_start, pos_end,
            liftover_grch,
            lifted_start, lifted_end,
            lifted_start, lifted_end
        ])

    else:
        conditions.append(sql.SQL("""
            (
                p.grch = %s
                AND (
                    (p.left_primer_start >= %s AND p.left_primer_end <= %s)
                    OR
                    (p.right_primer_start >= %s AND p.right_primer_end <= %s)
                )
            )
        """))

        values.extend([
            grch_val,
            pos_start, pos_end,
            pos_start, pos_end
        ])
    return conditions, values


def append_variant_liftover(use_variant_liftover, conditions, values,
                            grch_val, variant_pos, liftover_variant_grch,
                            lifted_variant_pos) -> tuple[list[sql.SQL], list[Any]]:
    """
    Append the variant liftover conditions and values for query
    Args:
        use_variant_liftover: bool to use liftover or not
        conditions: sql conditions for query
        values: input values for query
        grch_val: original GRCh
        variant_pos: original variant position
        liftover_variant_grch: liftover GRCh for variant
        lifted_variant_pos: liftover genomic position for variant

    Returns:
        Tuple of list of SQL and values if any
    """
    if use_variant_liftover:
        conditions.append(sql.SQL("""
                                (
                                    (p.grch = %s AND p.left_primer_end <= %s - 25  AND p.right_primer_start >= %s + 25)

                                    OR

                                    (p.grch = %s AND p.left_primer_end <= %s - 25 AND p.right_primer_start >= %s + 25)

                                    OR

                                    (p.grch = %s AND p.left_primer_end <= %s - 25 AND p.right_primer_start = 1 AND p.left_primer_end <> 1)

                                    OR

                                    (p.grch = %s AND p.left_primer_end = 1 AND p.right_primer_start >= %s + 25 AND p.right_primer_start <> 1)

                                    OR

                                    (p.grch = %s AND p.left_primer_end <= %s - 25 AND p.right_primer_start = 1 AND p.left_primer_end <> 1)

                                    OR

                                    (p.grch = %s AND p.left_primer_end = 1 AND p.right_primer_start >= %s + 25 AND p.right_primer_start <> 1)
                                )
                                """))
        values.extend([
                        grch_val, variant_pos, variant_pos,
                        liftover_variant_grch, lifted_variant_pos, lifted_variant_pos,

                        grch_val, variant_pos,
                        grch_val, variant_pos,

                        liftover_variant_grch, lifted_variant_pos,
                        liftover_variant_grch, lifted_variant_pos,
                    ])
    else:
        conditions.append(sql.SQL("""
                                    (
                                        (
                                            p.grch = %s
                                            AND p.left_primer_end <= %s - 25
                                            AND p.right_primer_start >= %s + 25
                                        )
                                        OR
                                        (
                                            p.grch = %s
                                            AND p.left_primer_end <= %s - 25
                                            AND p.right_primer_start = 1
                                            AND p.left_primer_end <> 1
                                        )
                                        OR
                                        (
                                            p.grch = %s
                                            AND p.left_primer_end = 1
                                            AND p.right_primer_start >= %s + 25
                                            AND p.right_primer_start <> 1
                                        )
                                    )
                                    """))
        values.extend([
                    grch_val, variant_pos, variant_pos,
                    grch_val, variant_pos,
                    grch_val, variant_pos
                ])
    return conditions, values


def search_prada_primer(db_schema, dbname, user, password, host,
                        chr_val, gene_val, primer_name, primer_id,
                        validation_val, grch_val, tray, archive_val,
                        notes, variant_pos=None, pos_start=None,
                        pos_end=None, start_date=None,
                        end_date=None) -> tuple[list[dict[str, Any]],
                                                str | None]:
    """
    Query primers in PRADA database with given inputs
    Various combinations could be provided to do query
    If the variant position, start or end position are given,
    the liftover is done to retrieve data from both GRCh37 and GRCh38.
    If the start and end Pos are given to do query, the primers that fall
    within start and end position are highlighted with green
    Args:
        db_schema: database schema
        dbname: database name
        user: username for database 
        password: pw for database
        host: database host
        chr_val: chromosome
        gene_val: gene
        primer_name: primer name
        primer_id: primer ID
        validation_val: validation status of primer
        grch_val: GRCh
        tray: primer tray
        archive_val: archiving status of primer
        notes: notes from user about primer
        variant_pos: position of variant
    Returns:
        query output results and message
    """
    # check if inputs are valid
    validate_position_inputs(variant_pos, pos_start, pos_end, grch_val)
    # get base query
    base_query = get_base_query(db_schema)
    # get the query conditions provided by user
    conditions, values = get_condition_value(chr_val, gene_val, validation_val,
                                             archive_val, notes, primer_name,
                                             primer_id, tray, start_date,
                                             end_date)
    # initiate params
    msg_to_return = None
    if grch_val == "":
        grch_val = None

    chrom = f"chr{chr_val}" if chr_val else None
    (
        liftover_grch,
        lifted_start,
        lifted_end,
        liftover_variant_grch,
        lifted_variant_pos,
    ) = (None, None, None, None, None)
    if grch_val:
        grch_val = str(grch_val)
        liftover_variant_grch = liftover_grch = "38" if grch_val == "37" else "37"
    # check if liftover is required or not
    (pos_start, pos_end,
     msg_to_return, variant_pos) = check_liftover_required(pos_start, pos_end,
                                                           msg_to_return,
                                                           variant_pos)

    # range position lift over
    if pos_start is not None and pos_end is not None and chrom:
        lifted_start, msg_to_return = safe_liftover(chrom, pos_start, grch_val)
        lifted_end, msg_to_return = safe_liftover(chrom, pos_end, grch_val)

    use_liftover = lifted_start is not None and lifted_end is not None
    if not use_liftover:
        liftover_grch = grch_val
        lifted_start = pos_start
        lifted_end = pos_end
    # append filters if range liftover is successful
    if pos_start is not None and pos_end is not None:

        conditions, values = append_pos_liftover(use_liftover, conditions,
                                                 values, grch_val, pos_start,
                                                 pos_end, liftover_grch,
                                                 lifted_start, lifted_end)

    #  variant pos lift over
    if variant_pos is not None and chrom:
        lifted_variant_pos, msg_to_return = safe_liftover(chrom, variant_pos,
                                                          grch_val)

    use_variant_liftover = lifted_variant_pos is not None
    # append filters if variant liftover is successful
    if variant_pos is not None:
        conditions, values = append_variant_liftover(use_variant_liftover,
                                                     conditions, values,
                                                     grch_val, variant_pos,
                                                     liftover_variant_grch,
                                                     lifted_variant_pos)

    # Only apply grch filter when not using liftover
    if not use_liftover and not use_variant_liftover:
        if grch_val:
            conditions.append(sql.SQL("p.grch = %s"))
            values.append(grch_val)
        else:
            conditions.append(sql.SQL("p.grch IN (%s, %s)"))
            values.extend(["37", "38"])
    # build SQL query
    query = base_query

    if conditions:
        query += sql.SQL(" WHERE ") + sql.SQL(" AND ").join(conditions)
    query += sql.SQL(" ORDER BY insert_time ASC, primer_id ASC LIMIT 10000")
    # Highlight values
    select_values = get_highlighted_value(pos_start, pos_end, grch_val,
                                          liftover_grch, lifted_start,
                                          lifted_end)
    with get_postgres_connection(dbname, user, password, host) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query, select_values + values)
            rows = cursor.fetchall()
            cols = [d[0] for d in cursor.description]
            result = [dict(zip(cols, row)) for row in rows]

    return result, msg_to_return


def query_moka_by_id(dbname, user, password, host,
                     primer_name) -> list[tuple[Any, ...]]:
    """
    Query MOKA legacy data with primer ID.
    Args:
        dbname: database name
        user: username for database
        password: password for database
        host: host of database
        primer_name: primer name
    Returns:
        Query results if any
    """
    query = """
    SELECT DISTINCT
        pa."PrimerName",
        c."Chr",
        pa."Start19",
        pa."Stop19",
        pa."ForwardSeq",
        pa."ReverseSeq",
        ftag."Item" AS "ForwardTag",
        rtag."Item" AS "ReverseTag",
        pa."Mix",
        pa."Notes",
        pa."TestResultNotes",
        pa."FFreezer",
        pa."FTray",
        pa."FGrid",
        pa."RFreezer",
        pa."RTray",
        pa."RGrid",
        pa."DateOrdered",
        pa."Manufacturer"
    FROM "moka_legacy"."PrimerAmplicon" pa
    INNER JOIN "moka_legacy"."Chromosome" c
        ON pa."ChromosomeID" = c."ChrID"
    INNER JOIN "moka_legacy"."Item" ftag
        ON pa."FTagName" = ftag."ItemID"
    INNER JOIN "moka_legacy"."Item" rtag
        ON pa."RTagName" = rtag."ItemID"
    WHERE pa."PrimerName" = %s
      AND pa."Status" = 1202218832;
    """

    with get_postgres_connection(dbname, user, password, host) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query, (primer_name,))
            return cursor.fetchall()


def query_moka_by_position(dbname, user, password, host,
                           chromosome, position) -> list[tuple[Any, ...]]:
    """
    Query MOKA legacy data with variant position.
    Args:
        dbname: database name
        user: username for database
        password: password for database
        host: host of database
        chromosome: chromosome number
        position: genomic position
    Returns:
        Query results if any
    """

    query = """
    SELECT DISTINCT
        pa."PrimerName",
        c."Chr",
        pa."Start19",
        pa."Stop19",
        pa."ForwardSeq",
        pa."ReverseSeq",
        ftag."Item" AS "ForwardTag",
        rtag."Item" AS "ReverseTag",
        pa."Mix",
        pa."Notes",
        pa."TestResultNotes",
        pa."FFreezer",
        pa."FTray",
        pa."FGrid",
        pa."RFreezer",
        pa."RTray",
        pa."RGrid",
        pa."DateOrdered",
        pa."Manufacturer"
    FROM "moka_legacy"."PrimerAmplicon" pa
    INNER JOIN "moka_legacy"."Chromosome" c
        ON pa."ChromosomeID" = c."ChrID"
    INNER JOIN "moka_legacy"."Item" ftag
        ON pa."FTagName" = ftag."ItemID"
    INNER JOIN "moka_legacy"."Item" rtag
        ON pa."RTagName" = rtag."ItemID"
    WHERE
        c."Chr" = %s
        AND (pa."Start19" + LENGTH(pa."ForwardSeq") + 25) < %s
        AND (pa."Stop19" - LENGTH(pa."ReverseSeq") - 25) > %s
        AND pa."Status" = 1202218832;
    """

    with get_postgres_connection(dbname, user, password, host) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query, (chromosome, position, position))
            return cursor.fetchall()


def query_moka_approved(dbname, user, password,
                        host, filters=None) -> list[tuple[Any, ...]]:
    """
    Query all approved primers from MOKA legacy data
    Filters can be applied to select specific primers
    Args:
        dbname: database name
        user: username for database
        password: password for database
        host: host of database
        filters: filter applied in output table

    Returns:
        Query results if any
    """

    query = """
    SELECT
        pa."PrimerName",
        c."Chr",
        pa."Start19",
        pa."Stop19",
        pa."ForwardSeq",
        pa."ReverseSeq",
        rtag."Item" AS "ReverseTag",
        ftag."Item" AS "ForwardTag",
        pa."Mix",
        pa."Notes",
        pa."TestResultNotes",
        pa."FFreezer",
        pa."FTray",
        pa."FGrid",
        pa."RFreezer",
        pa."RTray",
        pa."RGrid",
        pa."DateOrdered",
        pa."Manufacturer"

    FROM "moka_legacy"."PrimerAmplicon" pa

    INNER JOIN "moka_legacy"."Chromosome" c
        ON pa."ChromosomeID" = c."ChrID"

    INNER JOIN "moka_legacy"."Item" rtag
        ON pa."RTagName" = rtag."ItemID"

    INNER JOIN "moka_legacy"."Item" ftag
        ON pa."FTagName" = ftag."ItemID"

    WHERE pa."Status" = 1202218832
    """

    params = []

    if filters is None:
        filters = {}

    # add filters
    if filters.get("amplicon_id"):
        query += ' AND pa."AmpliconID" = %s'
        params.append(filters["amplicon_id"])

    if filters.get("chromosome"):
        query += ' AND c."Chr" = %s'
        params.append(filters["chromosome"])

    if filters.get("start"):
        query += ' AND pa."Start19" = %s'
        params.append(filters["start"])

    if filters.get("stop"):
        query += ' AND pa."Stop19" = %s'
        params.append(filters["stop"])

    if filters.get("primer_name"):
        query += ' AND pa."PrimerName" = %s'
        params.append(filters["primer_name"])

    if filters.get("manufacturer"):
        query += ' AND pa."Manufacturer" = %s'
        params.append(filters["manufacturer"])

    # non-specific filters

    if filters.get("notes"):
        query += ' AND pa."Notes" ILIKE %s'
        params.append("%" + filters["notes"] + "%")

    if filters.get("forward_seq"):
        query += ' AND pa."ForwardSeq" ILIKE %s'
        params.append("%" + filters["forward_seq"] + "%")

    if filters.get("reverse_seq"):
        query += ' AND pa."ReverseSeq" ILIKE %s'
        params.append("%" + filters["reverse_seq"] + "%")

    if filters.get("mix"):
        query += ' AND pa."Mix" ILIKE %s'
        params.append("%" + filters["mix"] + "%")

    if filters.get("r_tray"):
        query += ' AND pa."RTray" ILIKE %s'
        params.append("%" + filters["r_tray"] + "%")

    if filters.get("r_freezer"):
        query += ' AND pa."RFreezer" ILIKE %s'
        params.append("%" + filters["r_freezer"] + "%")

    if filters.get("f_grid"):
        query += ' AND pa."FGrid" ILIKE %s'
        params.append("%" + filters["f_grid"] + "%")

    if filters.get("f_tray"):
        query += ' AND pa."FTray" ILIKE %s'
        params.append("%" + filters["f_tray"] + "%")

    if filters.get("f_freezer"):
        query += ' AND pa."FFreezer" ILIKE %s'
        params.append("%" + filters["f_freezer"] + "%")

    if filters.get("test_result_notes"):
        query += ' AND pa."TestResultNotes" ILIKE %s'
        params.append("%" + filters["test_result_notes"] + "%")

    if filters.get("r_grid"):
        query += ' AND pa."RGrid" ILIKE %s'
        params.append("%" + filters["r_grid"] + "%")

    if filters.get("reverse_tag"):
        query += ' AND rtag."Item" ILIKE %s'
        params.append("%" + filters["reverse_tag"] + "%")

    if filters.get("forward_tag"):
        query += ' AND ftag."Item" ILIKE %s'
        params.append("%" + filters["forward_tag"] + "%")

    # date
    if filters.get("date_ordered"):
        query += ' AND CAST(pa."DateOrdered" AS TEXT) ILIKE %s'
        params.append("%" + filters["date_ordered"] + "%")

    with get_postgres_connection(dbname, user, password, host) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query, params)
            return cursor.fetchall()
