#!/usr/bin/env python3
"""Merge BRCA2 panel-C hypomorphic variants with 6k_dms_dataset.xlsx's
benign/pathogenic classifications into a single 3-class wide-format CSV.

Precedence: hypomorphic wins. If a variant is labeled hypomorphic (from the
panel-C source) it is excluded from the benign/pathogenic classes even if
6k_dms_dataset.xlsx also classifies it.

6k_dms_dataset.xlsx handling (see README for full rationale):
- 'AA.change' uses a 'p.' prefixed single-letter notation (e.g. 'p.L2480S').
  Non-missense entries (Intronic, synonymous wt==mut) are excluded.
- 'Classification' uses graded ACMG terms (Benign Very Strong/Strong/
  Moderate/Supporting, same for Pathogenic, plus Uncertain). Benign/
  Pathogenic variants are collapsed to their base label; Uncertain rows are
  dropped (doesn't fit a 3-class benign/pathogenic/hypomorphic scheme).
- Duplicate genomic-level entries occasionally give the SAME protein change
  conflicting benign/pathogenic calls; these variants are dropped entirely
  and reported, rather than guessing which call is right.

Usage:
    python scripts/build_brca2_3class_variants.py \\
        --hypomorphic-csv /scratch/jv2807/sounak_brca2/dataset/brca2_panelC_missense_variants.csv \\
        --dms-xlsx /scratch/jv2807/sounak_brca2/dataset/6k_dms_dataset.xlsx \\
        --output /scratch/jv2807/sounak_brca2/dataset/brca2_panelC_and_6k_3class.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from protein_vis.variants import VariantParseError, parse_variant, pivot_long_to_wide  # noqa: E402


def parse_aa_change(raw) -> object | None:
    """Parse 6k_dms_dataset's 'p.LxxxY' notation, or None if not a simple
    missense change (Intronic, synonymous, indel, etc.)."""
    if not isinstance(raw, str):
        return None
    s = raw.strip()
    if not s.startswith("p."):
        return None
    try:
        v = parse_variant(s[2:])
    except VariantParseError:
        return None
    if v.wt == v.mut:
        return None
    return v


def collapse_classification(classification) -> str | None:
    if not isinstance(classification, str):
        return None
    if classification.startswith("Benign"):
        return "benign"
    if classification.startswith("Pathogenic"):
        return "pathogenic"
    return None  # Uncertain or anything else -- dropped


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hypomorphic-csv", required=True)
    ap.add_argument("--dms-xlsx", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    hypo_df = pd.read_csv(args.hypomorphic_csv)
    hypo_variants = set(hypo_df.iloc[:, 0].dropna().astype(str).str.strip().str.upper())
    print(f"hypomorphic source: {len(hypo_variants)} unique variant(s)")

    dms = pd.read_excel(args.dms_xlsx, sheet_name="Sheet 1", header=1)
    dms["parsed"] = dms["AA.change"].apply(parse_aa_change)
    missense = dms[dms["parsed"].notna()].copy()
    missense["variant_str"] = missense["parsed"].apply(lambda v: f"{v.wt}{v.pos}{v.mut}")
    print(f"6k_dms_dataset.xlsx: {len(dms)} total row(s), {len(missense)} missense row(s)")

    missense["collapsed"] = missense["Classification"].apply(collapse_classification)
    classified = missense[missense["collapsed"].notna()]
    print(f"dropped {len(missense) - len(classified)} Uncertain/unclassified missense row(s)")

    conflict_counts = classified.groupby("variant_str")["collapsed"].nunique()
    conflicting = set(conflict_counts[conflict_counts > 1].index)
    print(f"dropping {len(conflicting)} variant(s) with conflicting benign/pathogenic calls "
          f"across duplicate entries")

    resolved = (
        classified[~classified["variant_str"].isin(conflicting)]
        .drop_duplicates("variant_str")
    )

    overridden = resolved[resolved["variant_str"].isin(hypo_variants)]
    if len(overridden):
        print(f"{len(overridden)} variant(s) present in both hypomorphic and 6k_dms_dataset "
              f"sources -- keeping hypomorphic label only")
    resolved = resolved[~resolved["variant_str"].isin(hypo_variants)]

    long_rows = [{"class_name": "hypomorphic", "raw": v} for v in sorted(hypo_variants)] + [
        {"class_name": row["collapsed"], "raw": row["variant_str"]}
        for _, row in resolved.iterrows()
    ]
    long_df = pd.DataFrame(long_rows)
    wide_df = pivot_long_to_wide(long_df, variant_col="raw", class_col="class_name")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    wide_df.to_csv(args.output, index=False)
    print(f"wrote {args.output}")
    print(wide_df.count().to_dict())


if __name__ == "__main__":
    main()
