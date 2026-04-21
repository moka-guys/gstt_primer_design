CREATE OR REPLACE FUNCTION primer_tool.insert_primer_with_batch(
    -- Primer fields
    p_chr VARCHAR,
    p_start_pos INTEGER,
    p_end_pos INTEGER,
    p_grch INTEGER,
    p_variant VARCHAR,
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
    p_passedvalidation TEXT,
    p_mix VARCHAR,
    p_dilute_time VARCHAR,
    p_tray VARCHAR,
    p_freezer VARCHAR,
    p_grid_fw VARCHAR,
    p_grid_rv VARCHAR
)
RETURNS INTEGER
LANGUAGE plpgsql
AS $$
DECLARE
    v_primer_id INTEGER;
    v_passed yes_no;
BEGIN

    IF p_passedvalidation NOT IN ('Yes', 'No') THEN
        RAISE EXCEPTION 'passed_validation must be Yes or No';
    END IF;

    v_passed := p_passedvalidation::yes_no;


    -- Primer info (no duplicate, no update allowed)
    INSERT INTO primer_tool.primers (
        chr,
        start_pos,
        end_pos,
        grch,
        variant,
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
        p_variant,
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
        start_pos,
        end_pos,
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
    RETURNING primer_id INTO v_primer_id;

    -- If conflict happened, fetch existing primer_id
    IF v_primer_id IS NULL THEN
        SELECT primer_id INTO v_primer_id
        FROM primer_tool.primers
        WHERE chr = p_chr
          AND start_pos = p_start_pos
          AND end_pos = p_end_pos
          AND grch = p_grch
          AND left_primer_seq = p_left_primer_seq
          AND right_primer_seq = p_right_primer_seq
          AND left_primer_start = p_left_primer_start
          AND left_primer_end = p_left_primer_end
          AND right_primer_start = p_right_primer_start
          AND right_primer_end = p_right_primer_end
          AND product_size = p_product_size
          AND gene = p_gene;
    END IF;

    -- =========================
    -- Batch info (can duplicate, update allowed for some)
    -- =========================
    INSERT INTO primer_tool.primer_batches (
        primer_id,
        tagged_left,
        tagged_right,
        tag,
        notes,
        passed_validation,
        mix,
        dilute_time,
        tray,
        freezer,
        grid_fw,
        grid_rv
    )
    VALUES (
        v_primer_id,
        p_tagged_left,
        p_tagged_right,
        p_tag,
        p_notes,
        v_passed,
        p_mix,
        p_dilute_time,
        p_tray,
        p_freezer,
        p_grid_fw,
        p_grid_rv
    );

    RETURN v_primer_id;

END;
$$;