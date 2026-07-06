"""Parsing and validation of protein variant lists.

Input CSVs are "wide" format: one column per class, column headers are class
names, and cell values are missense variants in compact [WT][pos][MUT]
notation (e.g. "L2484P"). Columns may be ragged (different lengths) --
pandas NaN-pads short columns automatically on read_csv.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

VARIANT_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")

# Standard single-letter amino acid alphabet.
AA1 = set("ACDEFGHIKLMNPQRSTVWY")


class VariantParseError(ValueError):
    pass


class VariantValidationError(ValueError):
    pass


@dataclass(frozen=True)
class Variant:
    raw: str
    wt: str
    pos: int
    mut: str


def parse_variant(s: str) -> Variant:
    """Parse a single missense variant string like 'L2484P'.

    Raises VariantParseError with the offending string on malformed input
    (e.g. frameshift/deletion notations like 'fs' or stop-codon '*', which
    some real variant-calling exports contain but this pipeline doesn't
    support) so failures are actionable rather than a downstream KeyError.
    """
    raw = s.strip()
    match = VARIANT_RE.match(raw)
    if not match:
        raise VariantParseError(
            f"could not parse variant {raw!r} as [WT][position][MUT] "
            f"(e.g. 'L2484P') -- frameshift/indel/stop-codon notations are "
            f"not supported"
        )
    wt, pos_str, mut = match.groups()
    wt, mut = wt.upper(), mut.upper()
    if wt not in AA1 or mut not in AA1:
        raise VariantParseError(
            f"variant {raw!r} has a non-standard amino acid letter "
            f"(wt={wt!r}, mut={mut!r})"
        )
    return Variant(raw=raw, wt=wt, pos=int(pos_str), mut=mut)


def load_variant_table(csv_path: str | Path) -> pd.DataFrame:
    """Load a wide-format variant CSV into a long DataFrame.

    Returns columns: class_name, raw, wt, pos, mut.
    """
    df = pd.read_csv(csv_path)
    if df.shape[1] == 0:
        raise VariantValidationError(f"{csv_path}: no columns found")

    long_df = df.melt(var_name="class_name", value_name="raw").dropna(subset=["raw"])
    long_df["raw"] = long_df["raw"].astype(str).str.strip()
    long_df = long_df[long_df["raw"] != ""]

    errors: list[str] = []
    records = []
    for class_name, raw in zip(long_df["class_name"], long_df["raw"]):
        try:
            v = parse_variant(raw)
        except VariantParseError as exc:
            errors.append(f"[{class_name}] {exc}")
            continue
        records.append(
            {"class_name": class_name, "raw": v.raw, "wt": v.wt, "pos": v.pos, "mut": v.mut}
        )

    if errors:
        raise VariantValidationError(
            f"{csv_path}: {len(errors)} unparseable variant(s):\n" + "\n".join(errors)
        )

    result = pd.DataFrame.from_records(
        records, columns=["class_name", "raw", "wt", "pos", "mut"]
    )
    result = result.drop_duplicates(subset=["class_name", "raw"]).sort_values("pos")
    return result.reset_index(drop=True)


def pivot_long_to_wide(
    long_df: pd.DataFrame, *, variant_col: str, class_col: str
) -> pd.DataFrame:
    """Inverse of load_variant_table's melt step.

    Turns (class, variant) pairs into the one-column-per-class wide CSV
    protein_vis expects, deduping within each class and padding ragged
    columns with blank cells. Used by data-prep scripts that build a
    variants CSV from external sources (e.g. merging multiple annotation
    files) rather than starting from one already in the wide format.
    """
    columns = {
        class_name: pd.Series(sorted(sub[variant_col].unique()))
        for class_name, sub in long_df.groupby(class_col)
    }
    return pd.DataFrame(columns)


def validate_against_sequence(
    df: pd.DataFrame, reference_seq: str, *, strict: bool = True
) -> tuple[pd.DataFrame, list[str]]:
    """Check each variant's WT residue against the canonical reference sequence.

    This checks the CSV against the protein's canonical (e.g. UniProt)
    sequence -- a data-entry sanity check -- and is distinct from checking a
    variant against a *structure's* sequence, which may legitimately differ
    (engineered constructs, orthologs, unresolved regions); see
    structure.align_to_reference for that.

    Positions out of range are always an error (row dropped). WT mismatches
    are an error (row dropped) if strict, else a warning (row kept, flagged
    via a new 'wt_mismatch' column).
    """
    messages: list[str] = []
    seq_len = len(reference_seq)
    keep_mask = []
    mismatch_flags = []

    for _, row in df.iterrows():
        pos, wt, raw, class_name = row["pos"], row["wt"], row["raw"], row["class_name"]
        if pos < 1 or pos > seq_len:
            messages.append(
                f"[{class_name}] {raw}: position {pos} out of range for "
                f"reference sequence of length {seq_len} -- dropped"
            )
            keep_mask.append(False)
            mismatch_flags.append(False)
            continue

        expected = reference_seq[pos - 1]
        if expected != wt:
            msg = (
                f"[{class_name}] {raw}: reference sequence has {expected!r} "
                f"at position {pos}, variant states wt={wt!r}"
            )
            if strict:
                messages.append(msg + " -- dropped (strict mode)")
                keep_mask.append(False)
                mismatch_flags.append(True)
                continue
            else:
                messages.append(msg + " -- kept (non-strict mode)")
                keep_mask.append(True)
                mismatch_flags.append(True)
                continue

        keep_mask.append(True)
        mismatch_flags.append(False)

    result = df.copy()
    result["wt_mismatch"] = mismatch_flags
    result = result[keep_mask].reset_index(drop=True)
    return result, messages
