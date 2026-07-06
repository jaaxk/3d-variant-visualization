#!/usr/bin/env python3
"""Reformat pkd1_dataset.xlsx into the wide-format CSV protein_vis expects,
using the 'Variant' and 'Pathogenicity mechanism' columns.

Non-missense entries (in-frame deletions like 'L4132del') and rows with no
pathogenicity-mechanism annotation (including the supplemental-table
footnote row that pandas reads as a data row) are excluded.

Usage:
    python scripts/build_pkd1_variants.py \\
        --input /scratch/jv2807/pkd1/data/pkd1_dataset.xlsx \\
        --output /scratch/jv2807/pkd1/data/pkd1_variants_labeled.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from protein_vis.variants import VariantParseError, parse_variant, pivot_long_to_wide  # noqa: E402


def is_missense(raw) -> bool:
    if not isinstance(raw, str):
        return False
    try:
        v = parse_variant(raw)
    except VariantParseError:
        return False
    return v.wt != v.mut


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--sheet", default="Sheet1")
    ap.add_argument("--header-row", type=int, default=3, help="0-indexed row containing column headers")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    df = pd.read_excel(args.input, sheet_name=args.sheet, header=args.header_row)
    print(f"read {len(df)} row(s) from {args.input}")

    keep = df["Variant"].apply(is_missense) & df["Pathogenicity mechanism"].notna()
    dropped = df[~keep]
    if len(dropped):
        print(f"dropping {len(dropped)} row(s) (non-missense variant or no classification):")
        print(dropped[["Variant", "Pathogenicity mechanism"]].to_string(index=False))

    long_df = df[keep].rename(
        columns={"Variant": "raw", "Pathogenicity mechanism": "class_name"}
    )[["raw", "class_name"]]
    print(f"kept {len(long_df)} classified missense variant(s):")
    print(long_df["class_name"].value_counts().to_string())

    wide_df = pivot_long_to_wide(long_df, variant_col="raw", class_col="class_name")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    wide_df.to_csv(args.output, index=False)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
