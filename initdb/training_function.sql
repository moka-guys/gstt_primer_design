CREATE OR REPLACE FUNCTION prada_training.insert_primer_with_batch(
    -- Primer fields
    p_chr VARCHAR,
    p_start_pos INTEGER,
    p_end_pos INTEGER,
    p_grch INTEGER,
    p_primer_name VARCHAR,
    p_left_primer_seq VARCHAR,
    p_right_primer_seq VARCHAR,
    p_left_primer_start INTEGER,
    p_left_primer_end INTEGER,
    p_right_primer_start INTEGER,
    p_right_primer_end INTEGER,
    p_product_size INTEGER,
    p_gene VARCHAR,

    -- Batch fields
    p_tagged_left VARCHAR,
    p_tagged_right VARCHAR,
    p_tag VARCHAR,
    p_notes VARCHAR,
    p_passed_validation TEXT,
    p_archive TEXT, 
    p_mix VARCHAR,
    p_arrival_date VARCHAR,
    p_tray VARCHAR,
    p_freezer VARCHAR,
    p_grid_fw VARCHAR,
    p_grid_rv VARCHAR,
    p_manufacturer VARCHAR,
    p_designer VARCHAR
)
RETURNS TABLE (
    upi INTEGER,
    primer_id INTEGER
)
LANGUAGE plpgsql
AS $$
DECLARE
    v_upi INTEGER;
    v_primer_id INTEGER;
    v_passed yes_no_status;
    v_archive yes_no;
BEGIN

    IF p_passed_validation NOT IN ('Pass', 'Fail', 'Not_Done') THEN
        RAISE EXCEPTION 'passed_validation must be Pass or Fail or Not_Done';
    END IF;

    v_passed := p_passed_validation::yes_no_status;
    v_archive := p_archive::yes_no;


    -- Primer info (no duplicate, no update allowed)
    INSERT INTO prada_training.primers (
        chr,
        start_pos,
        end_pos,
        grch,
        primer_name,
        left_primer_seq,
        right_primer_seq,
        left_primer_start,
        left_primer_end,
        right_primer_start,
        right_primer_end,
        product_size,
        gene
    )
    VALUES (
        p_chr,
        p_start_pos,
        p_end_pos,
        p_grch,
        p_primer_name,
        p_left_primer_seq,
        p_right_primer_seq,
        p_left_primer_start,
        p_left_primer_end,
        p_right_primer_start,
        p_right_primer_end,
        p_product_size,
        p_gene
    )
    ON CONFLICT (
        chr,
        primer_name,
        grch,
        left_primer_seq,
        right_primer_seq,
        left_primer_start,
        left_primer_end,
        right_primer_start,
        right_primer_end,
        product_size,
        gene
    )
    DO NOTHING
    RETURNING primers.upi INTO v_upi;

    -- If conflict happened, fetch existing primer_id
    IF v_upi IS NULL THEN
        SELECT p.upi INTO v_upi
        FROM prada_training.primers AS p
        WHERE p.chr = p_chr
          AND p.start_pos = p_start_pos
          AND p.end_pos = p_end_pos
          AND p.grch = p_grch
          AND p.primer_name = p_primer_name
          AND p.left_primer_seq = p_left_primer_seq
          AND p.right_primer_seq = p_right_primer_seq
          AND p.left_primer_start = p_left_primer_start
          AND p.left_primer_end = p_left_primer_end
          AND p.right_primer_start = p_right_primer_start
          AND p.right_primer_end = p_right_primer_end
          AND p.product_size = p_product_size
          AND p.gene = p_gene;
    END IF;

    -- =========================
    -- Batch info (can duplicate, update allowed for some)
    -- =========================
    INSERT INTO prada_training.primer_batches (
        upi,
        tagged_left,
        tagged_right,
        tag,
        notes,
        passed_validation,
        archive,
        mix,
        arrival_date,
        tray,
        freezer,
        grid_fw,
        grid_rv,
        manufacturer,
        designer
    )
    VALUES (
        v_upi,
        p_tagged_left,
        p_tagged_right,
        p_tag,
        p_notes,
        v_passed,
        v_archive,
        p_mix,
        p_arrival_date,
        p_tray,
        p_freezer,
        p_grid_fw,
        p_grid_rv,
        p_manufacturer,
        p_designer
    )

    RETURNING primer_batches.primer_id INTO v_primer_id;
    RETURN QUERY
    SELECT v_upi, v_primer_id;

END;
$$;