import pandas as pd
import pytest
from primer_design.primer3 import GeneratePrimer
import glob
import os


@pytest.fixture
def generator():
    return GeneratePrimer(
        app_datetimestr="20260806120000",
        config_path="test/data/test_config.json"
    )



def test_run_multiple_row_build(generator):
    input_file = "test/data/multiple_row_build.csv"

    result, error = generator.parse_input(input_file)

    assert error is None
    assert isinstance(result, pd.DataFrame)
    assert not result.empty

    expected_columns = {
        "Primer_Pair",
        "Left_Sequence",
        "Right_Sequence",
        "Pair_Product_Size",
        "Left_Start",
        "Left_End",
        "Right_Start",
        "Right_End",
        "Specificity",
        "FW_primer_snp",
        "RV_primer_snp",
        "snp_validity",
        "gene",
        "exon_num",
        "transcript",
        "order_tag",
        "chr",
        "GRCh",
        "start_POS",
        "end_POS",
        "primer_name",
    }

    assert expected_columns.issubset(result.columns)
    assert (result["Specificity"] == "valid").any()
    assert (result["snp_validity"] == "valid").any()
    assert result["primer_name"].nunique() == 2
    assert set(result["primer_name"]) == {"test1", "test2"}




def test_no_primer_output(generator):
    input_file = "test/data/no_primer_output.csv"

    result, error = generator.parse_input(input_file)
    assert error is None     
    assert result.empty


def test_build37(generator):
    input_file = "test/data/build37.csv"

    result, error = generator.parse_input(input_file)
    assert error is None     
    assert isinstance(result, pd.DataFrame)
    assert not result.empty

    expected_columns = {
                "Primer_Pair",
                "Left_Sequence",
                "Right_Sequence",
                "Pair_Product_Size",
                "Left_Start",
                "Left_End",
                "Right_Start",
                "Right_End",
                "Specificity",
                "FW_primer_snp",
                "RV_primer_snp",
                "snp_validity",
                "gene",
                "exon_num",
                "transcript",
                "order_tag",
                "chr",
                "GRCh",
                "start_POS",
                "end_POS",
                "primer_name",
            }

    assert expected_columns.issubset(result.columns)
    assert result["primer_name"].nunique() == 1
    snp_valid = result[result["snp_validity"] == "valid"]
    assert len(snp_valid) == 2
    snp_invalid = result[result["snp_validity"] == "invalid"]
    assert len(snp_invalid) == 1
    spe = result[result["Specificity"] == "valid"]
    assert len(spe) == 3
    assert result["gene"].iloc[0] == "POLG"
    assert int(result["exon_num"].iloc[0]) == 3
    assert result["transcript"].iloc[0] == "NM_001126131.2"


def test_build38(generator):
    input_file = "test/data/build38.csv"

    result, error = generator.parse_input(input_file)
    assert error is None     
    assert isinstance(result, pd.DataFrame)
    assert not result.empty

    expected_columns = {
                "Primer_Pair",
                "Left_Sequence",
                "Right_Sequence",
                "Pair_Product_Size",
                "Left_Start",
                "Left_End",
                "Right_Start",
                "Right_End",
                "Specificity",
                "FW_primer_snp",
                "RV_primer_snp",
                "snp_validity",
                "gene",
                "exon_num",
                "transcript",
                "order_tag",
                "chr",
                "GRCh",
                "start_POS",
                "end_POS",
                "primer_name",
            }

    assert expected_columns.issubset(result.columns)
    assert result["primer_name"].nunique() == 1
    snp_valid = result[result["snp_validity"] == "valid"]
    assert len(snp_valid) == 1
    snp_invalid = result[result["snp_validity"] == "invalid"]
    assert len(snp_invalid) == 1
    spe = result[result["Specificity"] == "valid"]
    assert len(spe) == 2
    assert result["gene"].iloc[0] == "ZFHX4"
    assert int(result["exon_num"].iloc[0]) == 4
    assert result["transcript"].iloc[0] == "NM_024721.5"


def test_gene_not_defined(generator):
    input_file = "test/data/gene_not_defined.csv"

    result, error = generator.parse_input(input_file)
    assert error is None     
    assert isinstance(result, pd.DataFrame)
    assert not result.empty

    expected_columns = {
                "Primer_Pair",
                "Left_Sequence",
                "Right_Sequence",
                "Pair_Product_Size",
                "Left_Start",
                "Left_End",
                "Right_Start",
                "Right_End",
                "Specificity",
                "FW_primer_snp",
                "RV_primer_snp",
                "snp_validity",
                "gene",
                "exon_num",
                "transcript",
                "order_tag",
                "chr",
                "GRCh",
                "start_POS",
                "end_POS",
                "primer_name",
            }

    assert expected_columns.issubset(result.columns)
    assert result["primer_name"].nunique() == 1
    assert result["gene"].iloc[0] == "not_defined"
    assert result["exon_num"].iloc[0] == "not_defined"
    assert result["transcript"].iloc[0] == "not_defined"

@pytest.fixture(autouse=True)
def cleanup_primer_files():
    yield

    print("Cleaning:", os.getcwd())

    for pattern in [
        "designed_primer_*.fa",
        "designed_primer_*.fa.sam",
    ]:
        for file in glob.glob(pattern):
            print("Deleting:", file)
            os.remove(file)
