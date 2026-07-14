#!/usr/bin/env python3
"""Assemble the multi-chain FASTA AlphaFold2-multimer expects for the
PKD1/PKD2 complex: one PKD1 record + N identical PKD2 records (the known
1:3 PKD1:PKD2 stoichiometry), all in a single FASTA file.

--pkd1-start lets the same script produce either the full-length PKD1
attempt or the fallback (trimmed-with-overlap) input if the full 7,207-
residue job doesn't finish in a reasonable time -- no separate script
needed, just a different flag value.

Usage:
    python build_complex_fasta.py --output ../fasta/pkd1_pkd2_complex_full.fasta
    python build_complex_fasta.py --pkd1-start 1900 \\
        --output ../fasta/pkd1_pkd2_complex_partial.fasta
"""

from __future__ import annotations

import argparse
from pathlib import Path

DEFAULT_PKD1_FASTA = "/scratch/jv2807/pkd1/structure_cache/uniprot/P98161.fasta"
DEFAULT_PKD2_FASTA = "/scratch/jv2807/pkd1/structure_cache/uniprot/Q13563.fasta"


def read_fasta_sequence(path: str | Path) -> str:
    lines = Path(path).read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pkd1-fasta", default=DEFAULT_PKD1_FASTA)
    ap.add_argument("--pkd2-fasta", default=DEFAULT_PKD2_FASTA)
    ap.add_argument(
        "--pkd1-start", type=int, default=1,
        help="1-based UniProt position to start PKD1's sequence from (trim "
             "everything before it). Default 1 = full-length, untrimmed. "
             "Use e.g. 1900 for the fallback run, giving generous overlap "
             "into 6A70's ~2221-4303 resolved region.",
    )
    ap.add_argument("--pkd2-copies", type=int, default=3, help="PKD1:PKD2 stoichiometry is 1:3.")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    pkd1_seq = read_fasta_sequence(args.pkd1_fasta)
    pkd2_seq = read_fasta_sequence(args.pkd2_fasta)

    if not (1 <= args.pkd1_start <= len(pkd1_seq)):
        raise ValueError(f"--pkd1-start {args.pkd1_start} out of range for a {len(pkd1_seq)}aa sequence")
    pkd1_trimmed = pkd1_seq[args.pkd1_start - 1:]

    records = [(f"PKD1_P98161_from{args.pkd1_start}", pkd1_trimmed)]
    records += [(f"PKD2_Q13563_copy{i + 1}", pkd2_seq) for i in range(args.pkd2_copies)]

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")

    total = sum(len(seq) for _, seq in records)
    print(f"wrote {out_path}")
    print(f"PKD1: residues {args.pkd1_start}-{len(pkd1_seq)} ({len(pkd1_trimmed)} aa)")
    print(f"PKD2: {args.pkd2_copies} copies x {len(pkd2_seq)} aa")
    print(f"total residues across all chains: {total}")


if __name__ == "__main__":
    main()
