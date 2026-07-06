#!/usr/bin/env python3
"""Turn a single-column variant CSV (no real class name) into a properly
headered wide-format CSV consumable by protein_vis.variants.load_variant_table.

Generic helper -- works for any single-column input, not BRCA2-specific.

Usage:
    python scripts/prepare_variant_csv.py \\
        --input /scratch/jv2807/sounak_brca2/dataset/brca2_panelC_missense_variants.csv \\
        --class-name hypomorphic \\
        --output /scratch/jv2807/sounak_brca2/dataset/brca2_panelC_hypomorphic.csv
"""

from __future__ import annotations

import argparse

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to single-column variant CSV.")
    parser.add_argument("--class-name", required=True, help="Class label for this variant set.")
    parser.add_argument("--output", required=True, help="Path to write the wide-format CSV.")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if df.shape[1] != 1:
        raise SystemExit(
            f"expected a single-column CSV, got columns: {list(df.columns)}"
        )

    out = pd.DataFrame({args.class_name: df.iloc[:, 0]})
    out.to_csv(args.output, index=False)
    print(f"wrote {len(out)} variant(s) under column '{args.class_name}' -> {args.output}")


if __name__ == "__main__":
    main()
