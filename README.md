# PRADA (PRimer and Amplicon Design Application)

This is a web application used to design primer(s) and insert into database

## How to set up locally

1. Clone github repo locally
2. Set .env file appropriately
3. Run `docker-compose up`. Docker image used should point to the latest tag (7030c31 as of Sept 2026)
4. Load the web app at `http://127.0.0.1:5000/`
5. To build DB locally, run sql files inside initdb inside postgres docker
6. Create username and password with 02_insert_default_users.sql
7. Log in to the web application and apply relevant features. 


# How primers are designed

Primer3 python package is used to design primers using given parameters. The designed primers were checked for specificity with Bowtie2 and the presence of SNP at designed primer regions are checked using the SNP vcf files downloaded from gnomad. 

# Features in PRADA
1. Design primers using primer3 automatically and insert designed primers into database
2. Insert manually designed primers into database
3. Query primers from PRADA
4. Query primers from MOKA legacy database

PRADA has two modes 1. Production and 2. Training. All features use the same source code but the database schema is different. Therefore, any primers designed/inserted during Training mode wil not affect on the Production database.
