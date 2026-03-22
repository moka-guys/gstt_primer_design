-- Create app_users table
CREATE TABLE IF NOT EXISTS app_users (
    id SERIAL PRIMARY KEY,
    username VARCHAR(50) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL
);

CREATE SCHEMA IF NOT EXISTS primer_tool;
CREATE TYPE yes_no AS ENUM ('Yes', 'No');
-- Create primer table within the schema if not exist
CREATE TABLE IF NOT EXISTS primer_tool.ordered_primers (
    primerid SERIAL PRIMARY KEY,
    chr VARCHAR(3) NOT NULL
    CHECK (
    chr IN ('X', 'Y')
    OR (chr ~ '^[0-9]+$' AND chr::int BETWEEN 1 AND 22)
    ),
    start_pos INTEGER CHECK (start_pos > 0 AND start_pos <= 300000000) NOT NULL,
    end_pos INTEGER CHECK (end_pos > 0 AND end_pos <= 300000000) NOT NULL,
    variant VARCHAR(15), 
    left_primer_seq VARCHAR(50),
    right_primer_seq VARCHAR(50),
    left_primer_start INTEGER CHECK (left_primer_start >= 0),
    left_primer_end INTEGER CHECK (left_primer_end >= 0),
    right_primer_start INTEGER CHECK (right_primer_start >= 0),
    right_primer_end INTEGER CHECK (right_primer_end >= 0),
    product_size INTEGER CHECK (product_size > 0 AND product_size < 10000),
    gene VARCHAR(20),    
    tagged_left VARCHAR(50),
    tagged_right VARCHAR(50),
    tag VARCHAR(10),
    GRCh INTEGER NOT NULL CHECK (GRCh IN (37, 38)),
    Notes VARCHAR(100),
    PassedValidation yes_no,
    insert_time TIMESTAMPTZ DEFAULT NOW(),
    CONSTRAINT ordered_primers_unique_all UNIQUE (
        chr,
        start_pos,
        end_pos,
        variant,
        left_primer_seq,
        right_primer_seq,
        left_primer_start,
        left_primer_end,
        right_primer_start,
        right_primer_end,
        product_size,
        gene,
        tagged_left,
        tagged_right,
        tag,
        GRCh
    )
);