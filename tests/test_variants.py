from pathlib import Path

import pytest

import pandas as pd

from protein_vis.variants import (
    VariantParseError,
    VariantValidationError,
    load_variant_table,
    parse_variant,
    pivot_long_to_wide,
    validate_against_sequence,
)

FIXTURES = Path(__file__).parent / "fixtures"
REFERENCE_SEQ = "EEEEMKTAYIAKQREEEE"  # positions 5-14 = MKTAYIAKQR


def test_parse_variant_ok():
    v = parse_variant("L2484P")
    assert v.wt == "L"
    assert v.pos == 2484
    assert v.mut == "P"
    assert v.raw == "L2484P"


def test_parse_variant_lowercase_normalized():
    v = parse_variant("l10a")
    assert v.wt == "L"
    assert v.mut == "A"


@pytest.mark.parametrize("bad", ["fs2484", "L2484*", "L2484fs", "123", "LP"])
def test_parse_variant_rejects_malformed(bad):
    with pytest.raises(VariantParseError):
        parse_variant(bad)


def test_load_variant_table_ragged_csv():
    df = load_variant_table(FIXTURES / "tiny_variants.csv")
    assert set(df["class_name"]) == {"classA", "classB"}
    assert len(df) == 5  # classA has 3, classB has 2 (one blank cell dropped)
    assert df.loc[df["raw"] == "M5V", "pos"].iloc[0] == 5


def test_load_variant_table_rejects_bad_variant(tmp_path):
    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("classA\nL10fs\n")
    with pytest.raises(VariantValidationError):
        load_variant_table(bad_csv)


def test_validate_against_sequence_strict_drops_mismatch():
    df = load_variant_table(FIXTURES / "tiny_variants.csv")
    # Corrupt one row's wt to force a mismatch.
    df.loc[df["raw"] == "M5V", "wt"] = "X"
    validated, messages = validate_against_sequence(df, REFERENCE_SEQ, strict=True)
    assert "M5V" not in validated["raw"].values
    assert any("M5V" in m for m in messages)


def test_validate_against_sequence_nonstrict_keeps_mismatch():
    df = load_variant_table(FIXTURES / "tiny_variants.csv")
    df.loc[df["raw"] == "M5V", "wt"] = "X"
    validated, _ = validate_against_sequence(df, REFERENCE_SEQ, strict=False)
    row = validated[validated["raw"] == "M5V"].iloc[0]
    assert bool(row["wt_mismatch"]) is True


def test_validate_against_sequence_out_of_range_dropped():
    df = load_variant_table(FIXTURES / "tiny_variants.csv")
    df.loc[df["raw"] == "M5V", "pos"] = 999
    validated, messages = validate_against_sequence(df, REFERENCE_SEQ, strict=True)
    assert 999 not in validated["pos"].values
    assert any("out of range" in m for m in messages)


def test_pivot_long_to_wide_pads_ragged_and_dedupes():
    long_df = pd.DataFrame(
        {
            "raw": ["L2484P", "S2516A", "S2516A", "K2551I"],
            "class_name": ["benign", "benign", "benign", "pathogenic"],
        }
    )
    wide = pivot_long_to_wide(long_df, variant_col="raw", class_col="class_name")
    assert set(wide.columns) == {"benign", "pathogenic"}
    assert wide["benign"].dropna().tolist() == ["L2484P", "S2516A"]  # deduped + sorted
    assert wide["pathogenic"].dropna().tolist() == ["K2551I"]
    assert pd.isna(wide["pathogenic"].iloc[1])  # ragged padding
