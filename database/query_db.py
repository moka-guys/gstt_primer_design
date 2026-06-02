from psycopg2 import sql
from primer_design.helper_function import get_postgres_connection, liftover, liftover_bed


def search_postgres(dbname, user, password, host,
                    chr_val, gene_val,  primer_name,
                    validation_val, grch_val, archive_val,
                    notes, pos_start=None, pos_end=None):

    conditions = []
    values = []
    msg_to_return = None
    base_query = sql.SQL("""
        SELECT
            p.unique_primer_id,
            b.primer_id,
            p.chr,
            p.grch,
            p.primer_name,
            b.passed_validation,
            b.archive,
            p.left_primer_seq,
            p.right_primer_seq,
            p.left_primer_start,
            p.left_primer_end,
            p.right_primer_start,
            p.right_primer_end,
            p.product_size,
            p.gene,
            b.tag,
            b.notes,
            b.mix,
            b.dilution_date,
            b.tray,
            b.freezer,
            b.grid_fw,
            b.grid_rv,
            b.manufacturer,
            b.insert_time,

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

        FROM primer_tool.primers p
        LEFT JOIN primer_tool.primer_batches b
        ON p.unique_primer_id = b.unique_primer_id
    """)

    try:
        pos_start = int(pos_start) if pos_start is not None else None
        pos_end = int(pos_end) if pos_end is not None else None
    except:
        msg_to_return = "Liftover is not done as start and end positions are not provided"
        pos_start = None
        pos_end = None

    if chr_val:
        conditions.append(sql.SQL("p.chr = %s"))
        values.append(chr_val)

    if gene_val:
        conditions.append(sql.SQL("p.gene = %s"))
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

    liftover_grch = None
    lifted_start = None
    lifted_end = None

    def safe_liftover(chrom, pos, build):
        try:
            lifted = liftover(chrom, pos, build)
            crossmap = liftover_bed(chrom, pos, pos, build)

            if lifted is None or crossmap is None:
                msg_to_return = "No corresponding lift over region is found"
                return None, msg_to_return

            if int(lifted) == int(crossmap):
                msg_to_return = "Lift over is checked by two different tools and results are consistent"
                return int(lifted), msg_to_return

            else:
                msg_to_return = "Lift over is not consistent for this region; therefore, region is not lifted over"
                return None, msg_to_return

        except Exception as e:
            print("error:", e)
            return None, str(e)
            
    if grch_val:
        grch_val = str(grch_val)
        liftover_grch = "38" if grch_val == "37" else "37"

        chrom = f"chr{chr_val}" if chr_val else None
        if pos_start is not None and pos_end is not None and chrom:
            lifted_start, msg_to_return = safe_liftover(chrom, pos_start, grch_val)
            lifted_end, msg_to_return = safe_liftover(chrom, pos_end, grch_val)

    use_liftover = lifted_start is not None and lifted_end is not None
    if not use_liftover:
        liftover_grch = grch_val
        lifted_start = pos_start
        lifted_end = pos_end
    # Highlight values
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
    if pos_start is not None and pos_end is not None:

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

    # Only apply grch filter when NOT using liftover
    if grch_val and not use_liftover:
        conditions.append(sql.SQL("p.grch = %s"))
        values.append(grch_val)

    query = base_query

    if conditions:
        query += sql.SQL(" WHERE ") + sql.SQL(" AND ").join(conditions)
    query += sql.SQL(" ORDER BY insert_time ASC, primer_id ASC LIMIT 100")
    conn = get_postgres_connection(dbname, user, password, host)
    cursor = conn.cursor()
    cursor.execute(query, select_values + values)

    rows = cursor.fetchall()
    cols = [d[0] for d in cursor.description]
    result = [dict(zip(cols, row)) for row in rows]
    conn.close()
    return result, msg_to_return

