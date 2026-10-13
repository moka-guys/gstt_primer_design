import os
import subprocess
import uuid
import pysam
import json
import primer3
import pandas as pd
import gffutils
from primer_design.helper_function import *
from dataclasses import dataclass


@dataclass
class InputParam:
    """Input parameters used for primer design."""
    chrom: str
    primer_name: str
    pos_start: int
    pos_end: int
    build: int
    tag: str
    opt_tm: float
    min_tm: float
    max_tm: float
    opt_gc: float
    min_gc: float
    max_gc: float
    min_p_size: int
    max_p_size: int
    opt_primer_size: int
    min_primer_size: int
    max_primer_size: int


class DesignPrimer:
    """Design primers for genomic regions."""
    def __init__(self, config_file, datetimestr,
                 input_param: InputParam) -> None:
        """Initialize primer design settings and reference resources.

        Args:
            config_file (str): Path to the configuration file.
            datetimestr (str): Timestamp used for logging.
            input_param (InputParam): Primer design parameters.
        """
        with open(config_file, "r") as file:
            self.config = json.load(file)
        self.input_param = input_param
        self.chr = self.input_param.chrom
        self.primer_name = self.input_param.primer_name
        self.pos_start = self.input_param.pos_start
        self.pos_end = self.input_param.pos_end
        self.build = self.input_param.build
        self.max_dist = self.config["design_param"]["max_dist"]
        self.tag = get_value(self.input_param.tag, "T1")
        self.opt_tm = get_value(self.input_param.opt_tm, self.config["design_param"]["primer_opt_tm"])
        self.min_tm = get_value(self.input_param.min_tm, self.config["design_param"]["primer_min_tm"])
        self.max_tm = get_value(self.input_param.max_tm, self.config["design_param"]["primer_max_tm"])
        self.opt_gc = get_value(self.input_param.opt_gc, self.config["design_param"]["primer_opt_gc"])
        self.min_gc = get_value(self.input_param.min_gc, self.config["design_param"]["primer_min_gc"])
        self.max_gc = get_value(self.input_param.max_gc, self.config["design_param"]["primer_max_gc"])
        self.min_product_size = get_value(self.input_param.min_p_size, self.config["design_param"]["min_product_size"])
        self.max_product_size = get_value(self.input_param.max_p_size, self.config["design_param"]["max_product_size"])
        self.min_primer_size = get_value(self.input_param.min_primer_size, self.config["design_param"]["primer_min_size"])
        self.max_primer_size = get_value(self.input_param.max_primer_size, self.config["design_param"]["primer_max_size"])
        self.opt_primer_size = get_value(self.input_param.opt_primer_size, self.config["design_param"]["primer_opt_size"])
        self.log_folder = self.config["directory"]["log_folder"]
        self.logger = get_log(self.log_folder, datetimestr)

        if self.build == 37:
            self.bowtie_ref = self.config["ref_b37"]["bowtie_ref"]
            self.ref_genome = self.config["ref_b37"]["genome_ref"]
            self.common_snp_exome = self.config["ref_b37"]["snp_ref_exome"]
            self.common_snp_genome = self.config["ref_b37"]["snp_ref_genome"]
            self.exon_db = self.config["ref_b37"]["exon_ref"]
            self.nc_pair = self.config["nc_data_b37"]

        elif self.build == 38:
            self.bowtie_ref = self.config["ref_b38"]["bowtie_ref"]
            self.ref_genome = self.config["ref_b38"]["genome_ref"]
            self.common_snp_exome = self.config["ref_b38"]["snp_ref_exome"]
            self.common_snp_genome = self.config["ref_b38"]["snp_ref_genome"]
            self.exon_db = self.config["ref_b38"]["exon_ref"]
            self.nc_pair = self.config["nc_data_b38"]

        else:
            self.logger.info("Invalid build is used in input file")
        self.logger.info(f"Primer design for {self.chr} {self.pos_start}-{self.pos_end}")
        self.logger.info(f"Design for Build {self.build} using {self.ref_genome}, "
                         f"{self.bowtie_ref}, and {self.exon_db}")

    def map_chr(self) -> str:
        """Map a chromosome name to its corresponding RefSeq NC accession.

        Returns:
            str: RefSeq NC accession number corresponding to the chromosome.
        """
        reverse_data = {v: k for k, v in self.nc_pair.items()}
        nc_number = reverse_data.get(self.chr)
        self.logger.info(f"nc_number for chr {self.chr} is {nc_number}")
        return nc_number

    def get_seq(self, nc_number, upstream_start, downstream_end) -> str:
        """Retrieve genomic sequence from the reference genome.

        Args:
            nc_number (str): RefSeq chromosome accession number.
            upstream_start (int): Start coordinate of the query region.
            downstream_end (int): End coordinate of the query region.

        Returns:
            str: Genomic sequence for the specified region.
        """
        fasta = pysam.FastaFile(self.ref_genome)
        # pysam is inclusive for given start
        # so no -1 here
        sequence = fasta.fetch(nc_number, upstream_start, downstream_end)
        return sequence

    def get_unique_exon(self, nc_number, start, end) -> tuple:
        """Retrieve a unique exon overlapping the specified region.

        Args:
            nc_number (str): RefSeq chromosome accession number.
            start (int): Start position of the query region.
            end (int): End position of the query region.

        Returns:
            tuple: Exon object and a boolean indicating whether
            exon was found.
        """
        db = gffutils.FeatureDB(self.exon_db, keep_order=True)
        exons = list(db.region(region=(nc_number, start, end),
                     completely_within=False, featuretype="exon"))
        unique_gene = []
        unique_exons = {}
        # get unique exons from found exons

        for exon in exons:
            gene_id = exon.attributes["gene_id"]
            key = (exon.start, exon.end, exon.strand, exon["gene_id"][0])
            if gene_id not in unique_gene:
                unique_gene.append(gene_id)
                unique_exons[key] = exon
        if len(unique_gene) == 1:
            exon = list(unique_exons.values())[0]
            exon_found = True
        else:
            exon = None
            exon_found = False
        return exon, exon_found

    def get_exon(self, nc_number) -> tuple:
        """Retrieve exon, gene, and transcript information for the target region.
        If the target position is not located within a unique exon, nearby
        positions (+/-15bp) are checked. If no unique exon is found, a dummy
        exon (i.e pos_start - 10, pos_end + 10) is returned.

        Args:
            nc_number (str): RefSeq chromosome accession number.

        Returns:
            tuple: Exon start position, exon end position, exon size, exon
            number, gene ID, and transcript ID.
        """
        exon_found = False
        offsets = [0, -15, +15]

        for offset in offsets:
            exon, exon_found = self.get_unique_exon(
                nc_number,
                self.pos_start + offset,
                self.pos_end + offset,
            )
            if exon_found:
                break
        else:
            # enter only if "for" loop does not break
            exon_start = self.pos_start - 10
            exon_end = self.pos_end + 10
            exon_size = exon_end - exon_start
            self.logger.info("0 or more than 1 exon found, so use dummy exon")
            return (exon_start, exon_end, exon_size, "not_defined",
                    "not_defined", "not_defined")
        if exon.end - exon.start > 450:
            exon_start = self.pos_start - 10
            exon_end = self.pos_end + 10
            exon_size = exon_end - exon_start
        else:
            exon_start = exon.start - 1  # to make 0 based
            exon_end = exon.end
            exon_size = exon_end - exon_start
        self.logger.info(f"Exon {exon['exon_number'][0]} on  "
                         f"{exon['gene_id'][0]}, {exon['transcript_id'][0]},  "
                         f"from {exon_start} to {exon_end}")
        return (exon_start, exon_end, exon_size, exon['exon_number'][0],
                exon['gene_id'][0], exon['transcript_id'][0])

    def get_snp(self, start, end) -> list:
        """Retrieve common SNPs from the ref SNP files.

        Args:
            start (int): Start position of the query region.
            end (int): End position of the query region.

        Returns:
            list: SNP records found within the specified region.
        """
        if self.build == 38:
            chrom = f"chr{self.chr}"
        else:
            chrom = self.chr
        snp_list = []
        for snp_file in [self.common_snp_genome, self.common_snp_exome]:
            with pysam.TabixFile(snp_file) as vcf:
                snp_list.extend(vcf.fetch(chrom, start, end))

        return snp_list

    def get_primer_value(self, primers, start, end) -> dict:
        """Extract primer-related key-value pairs from Primer3 output.

        Args:
            primers (dict): Dictionary of Primer3 results.
            start (str): Prefix used to select keys.
            end (str): Suffix used to select keys.

        Returns:
            dict: Filtered dictionary containing matching primer values.
        """

        key_value = {
            k: v for k, v in primers.items() if k.startswith(start) and k.endswith(end)
        }
        return key_value

    def run_primer3(self, sequence, gene_name, primer_space, exon_plus,
                    upstream_start, job_id) -> tuple:
        """Design PCR primers using Primer3.
        Primer candidates are generated and filtered to remove repetitive
        primer sequences before being written to a FASTA file.

        Args:
            sequence (str): Template sequence used for primer design.
            gene_name (str): Gene identifier
            primer_space (int): Distance allowed for primer placement.
            exon_plus (int): Target region size including flanking sequence.
            upstream_start (int): Genomic start coordinate of the sequence.
            job_id (str): Unique identifier for the design run.

        Returns:
            tuple: FASTA file handle, primer DataFrame, and a boolean
            indicating whether primers were successfully designed.
        """
        seq = {
            "SEQUENCE_TEMPLATE": sequence,
            "SEQUENCE_ID": gene_name,
        }
        self.logger.info(f"Exon used has length {exon_plus - 31} (including flanking regions)")
        self.logger.info(f"Seq used has length {len(sequence)}")
        self.logger.info(f"primer is designed for {seq}")

        param = {
            "PRIMER_MAX_SIZE": self.max_primer_size,
            "PRIMER_OPT_SIZE": self.opt_primer_size,
            "PRIMER_MIN_SIZE": self.min_primer_size,
            "PRIMER_TASK": "pick_pcr_primers",
            "PRIMER_EXPLAIN_FLAG": 1,
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PRODUCT_SIZE_RANGE": [[self.min_product_size, self.max_product_size]],
            "PRIMER_LOWERCASE_MASKING": 1,
            "PRIMER_MAX_NS_ACCEPTED": 0,
            "PRIMER_MAX_SELF_ANY": 8,
            "PRIMER_MAX_SELF_END": 3,
            "PRIMER_MIN_THREE_PRIME_DISTANCE": 3,
            "SEQUENCE_TARGET": [primer_space, exon_plus],
            "PRIMER_NUM_RETURN": self.config["design_param"]["primer_num_return"],
            "PRIMER_MIN_TM": self.min_tm,
            "PRIMER_OPT_TM": self.opt_tm,
            "PRIMER_MAX_TM": self.max_tm,
            "PRIMER_MIN_GC": self.min_gc,
            "PRIMER_OPT_GC": self.opt_gc,
            "PRIMER_MAX_GC": self.max_gc,
            "PRIMER_MAX_POLY_X": self.config["design_param"]["primer_max_poly_x"],
            "PRIMER_GC_CLAMP": self.config["design_param"]["primer_gc_clamp"],
         }
        self.logger.info(f"primer is designed with parameters {param}")
        primers = primer3.bindings.design_primers(seq, param)
        primer_key_value = {}
        # extract primer info from primer3 outputs
        for end in [
            "_SEQUENCE",
            "_PENALTY",
            "_TM",
            "_GC_PERCENT",
            "_ANY_TH",
            "_SELF_ANY_TH",
            "_HAIRPIN_TH",
            "_END_STABILITY",
            "_PRODUCT_SIZE",
        ]:
            primer_key_value.update(self.get_primer_value(primers, "PRIMER_", end))
        # extract primer positions
        primers_pos = {
            k: v for k, v in primers.items() if k.startswith("PRIMER_") and k[-1].isdigit()
        }
        primer_key_value.update(primers_pos)
        # check how many primers are designed
        primer_indices = set()
        for key in primer_key_value:
            if key.startswith("PRIMER_LEFT_") or key.startswith("PRIMER_RIGHT_"):
                idx = int(key.split("_")[2])
                primer_indices.add(idx)

        primer_indices = sorted(primer_indices)
        # put primer outputs into dict and then to df
        rows = []
        df = pd.DataFrame()
        for i in primer_indices:
            row = {
                "Primer_Pair": i,
                "Left_Sequence": primer_key_value.get(f"PRIMER_LEFT_{i}_SEQUENCE"),
                "Right_Sequence": primer_key_value.get(f"PRIMER_RIGHT_{i}_SEQUENCE"),
                "Left_Penalty": primer_key_value.get(f"PRIMER_LEFT_{i}_PENALTY"),
                "Right_Penalty": primer_key_value.get(f"PRIMER_RIGHT_{i}_PENALTY"),
                "Left_Tm": primer_key_value.get(f"PRIMER_LEFT_{i}_TM"),
                "Right_Tm": primer_key_value.get(f"PRIMER_RIGHT_{i}_TM"),
                "Left_GC": primer_key_value.get(f"PRIMER_LEFT_{i}_GC_PERCENT"),
                "Right_GC": primer_key_value.get(f"PRIMER_RIGHT_{i}_GC_PERCENT"),
                "Left_End_Stability": primer_key_value.get(
                    f"PRIMER_LEFT_{i}_END_STABILITY"
                ),
                "Right_End_Stability": primer_key_value.get(
                    f"PRIMER_RIGHT_{i}_END_STABILITY"
                ),
                "Left_ANY_TH": primer_key_value.get(f"PRIMER_LEFT_{i}_ANY_TH"),
                "Pair_Product_Size": primer_key_value.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
            }
            rows.append(row)
            df = pd.DataFrame(rows)

            df["Right_Start"] = ""
            df["Right_End"] = ""
            df["Left_Start"] = ""
            df["Left_End"] = ""

            for k, v in primers_pos.items():
                primer_num = int(k.split("_")[-1])
                if "RIGHT" in k:
                    absolute_end = upstream_start + v[0] + 1  # as in zippy
                    absolute_start = absolute_end - v[1]
                    df.loc[primer_num, "Right_Start"] = absolute_start
                    df.loc[primer_num, "Right_End"] = absolute_end
                else:
                    absolute_start = upstream_start + v[0]
                    absolute_end = absolute_start + v[1]
                    df.loc[primer_num, "Left_Start"] = absolute_start
                    df.loc[primer_num, "Left_End"] = absolute_end

        df_filtered = df[df.apply(self.keep_primer_pair, axis=1)].copy().reset_index(drop=True)
        df_filtered["Primer_Pair"] = range(len(df_filtered))

        if df_filtered.empty:
            primer_found = False
            self.logger.info("xxxxx 0 primer is designed by primer3 xxxxx")
            return None, None, primer_found

        with open(f"designed_primer_{job_id}.fa", "a") as f:
            for i in range(df_filtered.shape[0]):
                name = "PRIMER_Left_" + str(df_filtered["Primer_Pair"][i])
                f.write(">" + name + "|" + str(self.chr) + ":" + str(df_filtered["Left_Start"][i])
                        + "-" + str(df_filtered["Left_End"][i]) + "\n")
                f.write(df_filtered["Left_Sequence"][i] + "\n")
                name = "PRIMER_Right_" + str(df_filtered["Primer_Pair"][i])
                f.write(">" + name + "|" + str(self.chr) + ":" + str(df_filtered["Right_Start"][i])
                        + "-" + str(df_filtered["Right_End"][i]) + "\n")
                f.write(df_filtered["Right_Sequence"][i] + "\n")
        primer_found = True

        return f, df_filtered, primer_found

    def has_tandem_repeat(self, seq, unit_size, min_repeats) -> bool:
        """Determine whether a sequence contains tandem repeats.

        Args:
            seq (str): Primer sequence.
            unit_size (int): Size of the repeating motif.
            min_repeats (int): Minimum number of repeat units required.

        Returns:
            bool: True if a tandem repeat is detected, otherwise False.
        """
        seq = seq.upper()

        for i in range(len(seq) - unit_size * min_repeats + 1):
            motif = seq[i:i + unit_size]

            if motif * min_repeats in seq:
                return True

        return False

    def passes_repeat_filter(self, seq) -> bool:
        """Check whether a primer sequence passes repeat filtering.

        Args:
            seq (str): Primer sequence.

        Returns:
            bool: True if the sequence passes all repeat filters, otherwise
            False.
        """
        # Reject dinucleotide repeats (e.g. ACACACAC), min_repeat num inclusive
        if self.has_tandem_repeat(seq, unit_size=2, min_repeats=4):
            return False

        # Reject trinucleotide repeats (e.g. CAGCAGCAG)
        if self.has_tandem_repeat(seq, unit_size=3, min_repeats=3):
            return False

        # Reject tetra nucleotide repeats
        if self.has_tandem_repeat(seq, unit_size=4, min_repeats=3):
            return False

        return True

    def keep_primer_pair(self, row) -> bool:
        """Determine whether a primer pair passes repeat filtering.

        Args:
            row (pd.Series): Primer pair record.

        Returns:
            bool: True if both primers pass repeat filtering, otherwise
            False.
        """
        return (
            self.passes_repeat_filter(row["Left_Sequence"]) and
            self.passes_repeat_filter(row["Right_Sequence"])
        )

    def bowtie_mapping(self, primer_fa) -> pd.DataFrame:
        """Assess primer specificity using Bowtie2 alignments.

        Args:
            primer_fa (str): Path to the FASTA file containing designed
            primers.

        Returns:
            pd.DataFrame: DataFrame containing Bowtie2 alignment results.
        """
        bowtie_output = primer_fa + ".sam"

        if not os.path.exists(bowtie_output):
            with open(bowtie_output, "w") as outfile:
                proc = subprocess.check_call(
                    [
                        "bowtie2", "-f", "--end-to-end", "-p", "2",
                        "-k", str(50), "-L", "10", "-N", "1", "-D", "20",
                        "-R", "3", "-x", self.bowtie_ref, "-U", primer_fa,
                    ],
                    stdout=outfile,
                )
        rows = []
        with open(bowtie_output) as f:
            for line in f:
                if line.startswith("@"):
                    continue  # skip header lines
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue  # skip malformed lines
                row = {
                        "primer_pair": int(parts[0].split("|")[0].split("_")[-1]),
                        "name": parts[0].split("|")[0],
                        "strand": parts[1],
                        "chr": parts[2],
                        "pos_start": parts[3],
                        "cigar": parts[5],
                        "primer_seq": parts[9]
                    }
                rows.append(row)

            df = pd.DataFrame(rows)

        return df

    def update_df(self, primer_df, bowtie_df, nc_number) -> pd.DataFrame:
        """Update primer results with specificity information.

        Args:
            primer_df (pd.DataFrame): Designed primer information.
            bowtie_df (pd.DataFrame): Bowtie2 alignment results.
            nc_number (str): Expected chromosome accession number.

        Returns:
            pd.DataFrame: Updated primer DataFrame containing specificity
            classifications.
        """
        primer_df = primer_df[["Primer_Pair", "Left_Sequence",
                               "Right_Sequence", "Pair_Product_Size",
                               "Left_Start", "Left_End",
                               "Right_Start", "Right_End"]]

        primer_pair = list(primer_df["Primer_Pair"].unique())
        primer_df['Specificity'] = None
        # loop for each primer pair
        for i in range(len(primer_pair)):
            idx = primer_df.index[primer_df["Primer_Pair"] == i][0]
            temp_df = bowtie_df[bowtie_df["primer_pair"] == i]
            temp_df = temp_df.reset_index(drop=True)
            amplicon = []
            temp_df["pos_start"] = pd.to_numeric(temp_df["pos_start"], errors="coerce")

            result_list = []
            left = temp_df[temp_df["name"].str.contains("Left", case=False)]
            right = temp_df[temp_df["name"].str.contains("Right", case=False)]
            # check how many primer pair can generate amplicon
            # if left and right primers are not on the same chr
            # they can generate amplicon
            # if left and right primers are on the same chr, but
            # distance is greater than 2kb, cannot generate the amplicon
            # checking how many pair (between lt and rt) can generate amplicon
            for _, l_row in left.iterrows():
                for _, r_row in right.iterrows():
                    if l_row["chr"] == r_row["chr"]:  # only same chromosome
                        size = abs(r_row["pos_start"] + len(r_row["primer_seq"]) - l_row["pos_start"])
                        if size < self.max_dist:
                            result_list.append({
                                "primer_pair": i,
                                "chr": l_row["chr"],
                                "left_pos": l_row["pos_start"],
                                "right_pos": r_row["pos_start"],
                                "size": size
                            })

            amplicon = pd.DataFrame(result_list)
            row = primer_df.loc[primer_df["Primer_Pair"] == i].iloc[0]
            if len(amplicon) == 1:
                # check if the unique primer matches with req parameter
                amp = amplicon.iloc[0]

                if (
                    amp["chr"] == nc_number
                    and amp["size"] == row["Pair_Product_Size"]
                    and amp["left_pos"] - 1 == row["Left_Start"]
                    and amp["right_pos"] - 1 == row["Right_Start"]
                ):
                    primer_df['Specificity'][idx] = "valid"
                else:
                    primer_df['Specificity'][idx] = "invalid"

            else:
                primer_df['Specificity'][idx] = "invalid"

        return primer_df

    def classify_variant(self, ref, alt) -> list:
        """Classify variants as SNPs, INDELs, or complex variants.

        Args:
            ref (str): Reference allele.
            alt (str): Alternate alleles.

        Returns:
            list: Classification for each alternate allele.
        """
        alts = alt.split(',')
        classifications = []

        for alt_allele in alts:
            if len(ref) == 1 and len(alt_allele) == 1:
                cls = "SNP"
            elif len(ref) != len(alt_allele):
                cls = "INDEL"
            else:
                cls = "complex"
            classifications.append(cls)

        return classifications

    def has_snp(self, new_col, updated_df, p_start, p_end) -> pd.DataFrame:
        """Identify common variants within primer binding regions.
        Variants with allele frequency greater than or equal to 1% are
        reported for each primer.

        Args:
            new_col (str): Output column name.
            updated_df (pd.DataFrame): Primer DataFrame.
            p_start (str): Start position column name.
            p_end (str): End position column name.

        Returns:
            pd.DataFrame: Updated DataFrame containing SNP information.
        """
        updated_df[new_col] = None
        for i in range(updated_df.shape[0]):
            if updated_df["Specificity"][i] == "valid":  # check for primer with bowtie specificity passed
                snp_list = self.get_snp(updated_df[p_start][i], updated_df[p_end][i])
                snps = []
                if len(snp_list) > 0:
                    for s in snp_list:
                        cols = s.strip().split('\t')
                        info_field = cols[7]
                        info_dict = dict(item.split("=", 1) for item in info_field.split(";") if "=" in item)
                        af = info_dict.get("AF")
                        if af is not None:
                            try:
                                af_value = float(af)
                                if af_value >= 0.01:
                                    ref = cols[3]
                                    alt = cols[4]
                                    # check type of variant SNP or INDEL or else
                                    variant_class = self.classify_variant(ref, alt)
                                    snps.append(f"{cols[2]}-{variant_class} with AF:{af_value}")
                            except ValueError:
                                continue  # skip malformed AF
                    if snps:
                        updated_df[new_col][i] = snps
                    else:
                        updated_df[new_col][i] = "no common snp"
                else:
                    updated_df[new_col][i] = "no snp"

            else:
                updated_df[new_col][i] = "not checked"

        return updated_df

    def classify_snp(self, df) -> pd.DataFrame:
        """Determine primer validity based on detected variants.
        Primer pairs containing common variants are marked as invalid.

        Args:
            df (pd.DataFrame): Primer DataFrame containing SNP annotations.

        Returns:
            pd.DataFrame: Updated DataFrame containing SNP validity status.
        """
        df["snp_validity"] = None
        for i in range(df.shape[0]):
            if df["Specificity"][i] == "valid":
                fw_snp = make_list(df["FW_primer_snp"][i])
                rv_snp = make_list(df["RV_primer_snp"][i])
                total_snp = fw_snp + rv_snp
                if any("INDEL" in s for s in total_snp):
                    df["snp_validity"][i] = "invalid"

                elif len(total_snp) >= 1:
                    df["snp_validity"][i] = "invalid"

                else:
                    df["snp_validity"][i] = "valid"

        return df

    def design_primer(self) -> pd.DataFrame:
        """Design primers for the specified genomic region.
        Primer design is performed iteratively with increasing padding until
        a valid primer pair is identified or the maximum padding limit is
        reached.

        Returns:
            pd.DataFrame: DataFrame containing designed primer information.
        """
        # get nc number mapped with chr
        nc_number = self.map_chr()
        # get exon, gene and transcript info for given POS
        (exon_start, exon_end, exon_size,
         exon_num, gene, transcript) = self.get_exon(nc_number)
        # add flanking region
        exon_size = exon_size + (2*self.config["design_param"]["flank"])
        exon_plus = exon_size + 31  # +1 as last number is not inclusive in 0 based
        padding = self.config["design_param"]["padding"]
        primer_found = False
        valid_primer = False
        #order_df = pd.DataFrame()
        designed_primer = pd.DataFrame()

        while (padding <= self.config["design_param"]["max_padding"] and
                (not primer_found or not valid_primer)):
            self.logger.info(f"******Running primer3 with padding {padding}*****")
            primer_space = padding - 30
            upstream_start = exon_start - padding
            downstream_end = exon_end + padding
            if self.build == 38:
                sequence = self.get_seq(nc_number, upstream_start, downstream_end)
            else:
                sequence = self.get_seq(str(self.chr), upstream_start, downstream_end)
            job_id = uuid.uuid4().hex
            primer_fa, primer_df, primer_found = self.run_primer3(
                                                    sequence, gene,
                                                    primer_space, exon_plus,
                                                    upstream_start, job_id)
            padding = padding + 30
            # if primer3 generates any primer
            if primer_found:
                bowtie_df = self.bowtie_mapping(primer_fa.name)
                if self.build == 38:
                    updated_df = self.update_df(primer_df, bowtie_df, nc_number)
                else:
                    updated_df = self.update_df(primer_df, bowtie_df, str(self.chr))
                # check snp on FW and RV primers
                updated_df = self.has_snp("FW_primer_snp", updated_df, "Left_Start", "Left_End")
                updated_df = self.has_snp("RV_primer_snp", updated_df, "Right_Start", "Right_End")
                updated_df = self.classify_snp(updated_df)
                valid_df = updated_df[(updated_df["Specificity"] == "valid") & 
                                      (updated_df["snp_validity"] == "valid")]
                # if any designed primer is valid
                if valid_df.shape[0] > 0:
                    valid_primer = True
                    designed_primer = updated_df.copy()
                    designed_primer["gene"] = gene
                    designed_primer["exon_num"] = exon_num
                    designed_primer["transcript"] = transcript
                    designed_primer["order_tag"] = self.tag
                    designed_primer["chr"] = self.chr
                    designed_primer["GRCh"] = self.build
                    designed_primer["start_POS"] = int(self.pos_start)
                    designed_primer["end_POS"] = int(self.pos_end)
                    designed_primer["primer_name"] = self.primer_name
                    self.logger.info(f"****primer found with padding {padding - 30}****")
                    del_file([f"designed_primer_{job_id}.fa", f"designed_primer_{job_id}.fa.sam"])

                # if none of designed primer is valid and max padding not reach yet
                elif valid_df.shape[0] == 0 and padding <= self.config["design_param"]["max_padding"]:
                    valid_primer = False
                    del_file([f"designed_primer_{job_id}.fa", f"designed_primer_{job_id}.fa.sam"])
                    continue
                else:
                    valid_primer = False
                    self.logger.info("Valid primer not found till max padding. "
                                     "Try with different config for primer design")
                    del_file([f"designed_primer_{job_id}.fa", f"designed_primer_{job_id}.fa.sam"])

        return designed_primer


class GeneratePrimer:
    """Generate primer designs from input files."""
    def __init__(self, app_datetimestr, config_path=None) -> None:
        """Initialize primer generation settings.

        Args:
            app_datetimestr (str): Timestamp used for logging.
            config_path (str | None): Path to the configuration file.
        """
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))

        if config_path is None:
            config_path = os.path.join(BASE_DIR, "config.json")

        self.config_path = config_path
        self.datetimestr = app_datetimestr

    def parse_input(self, input_file) -> tuple[pd.DataFrame, str | None]:
        """Generate primer designs from an input CSV file.
        Each row of the input file is processed independently and all
        successful primer designs are combined into a single output
        DataFrame.

        Args:
            input_file (str): Path to the input CSV file.

        Returns:
            tuple[pd.DataFrame, str | None]: Generated primer DataFrame and
            an error message if processing fails, otherwise None.
        """
        input_param = parse_csv(input_file)
        if "error" in input_param:
            return pd.DataFrame(), input_param["error"]
        dfs = []
        for i in range(len(input_param["chrom"])):
            param_i = {k: v[i] for k, v in input_param.items() if k in InputParam.__annotations__}
            param_dict = InputParam(**param_i)
            start_design = DesignPrimer(self.config_path, self.datetimestr, param_dict)
            designed_primer = start_design.design_primer()
            dfs.append(designed_primer)

        dfs = [df for df in dfs if not df.empty]

        if dfs:
            df_all = pd.concat(dfs, ignore_index=True)
        else:
            df_all = pd.DataFrame()

        return df_all, None
