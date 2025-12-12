-- Create app_users table
CREATE TABLE IF NOT EXISTS app_users (
    id SERIAL PRIMARY KEY,
    username VARCHAR(50) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL
);

CREATE SCHEMA IF NOT EXISTS primer_tool_test;

-- Create primer table within the schema if not exist
CREATE TABLE IF NOT EXISTS primer_tool_test.ordered_primers (
    primerid SERIAL PRIMARY KEY,
    chr TEXT,
    position TEXT,
    variant TEXT,
    left_primer_seq TEXT,
    right_primer_seq TEXT,
    left_primer_start INTEGER,
    left_primer_end INTEGER,
    right_primer_start INTEGER,
    right_primer_end INTEGER,
    product_size INTEGER,
    gene TEXT,
    tagged_left TEXT,
    tagged_right TEXT,
    tag TEXT,
    GRCh INTEGER,
    insert_time TIMESTAMPTZ DEFAULT NOW(),
    CONSTRAINT ordered_primers_unique_all UNIQUE (
        chr,
        position,
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