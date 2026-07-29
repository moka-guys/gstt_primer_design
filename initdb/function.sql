CREATE OR REPLACE FUNCTION primer_tool.insert_primer_with_batch(
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
RETURNS INTEGER
LANGUAGE plpgsql
AS $$
DECLARE
    v_unique_primer_id INTEGER;
    v_passed yes_no_status;
    v_archive yes_no;
BEGIN

    IF p_passed_validation NOT IN ('Pass', 'Fail', 'Not_Done') THEN
        RAISE EXCEPTION 'passed_validation must be Pass or Fail or Not_Done';
    END IF;

    v_passed := p_passed_validation::yes_no_status;
    v_archive := p_archive::yes_no;


    -- Primer info (no duplicate, no update allowed)
    INSERT INTO primer_tool.primers (
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
    RETURNING unique_primer_id INTO v_unique_primer_id;

    -- If conflict happened, fetch existing primer_id
    IF v_unique_primer_id IS NULL THEN
        SELECT unique_primer_id INTO v_unique_primer_id
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
        unique_primer_id,
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
        v_unique_primer_id,
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
    );

    RETURN v_unique_primer_id;

END;
$$;