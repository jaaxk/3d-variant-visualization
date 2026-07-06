#!/bin/bash
# PKD1 (Polycystin-1, UniProt P98161) variants classified by pathogenicity
# mechanism (Benign / Nontrafficking / Function) -- see
# scripts/build_pkd1_variants.py for how the input CSV is built.
#
# Structure: PDB 6A70 chain B (human PKD1-PKD2 complex cryo-EM structure,
# Su et al. 2018 Science), covering only residues 3049-4169 of PKD1's 4303
# aa -- no AlphaFold model exists for full-length PKD1 either (excluded,
# likely due to size, same situation as BRCA2). ~42% of variants fall
# outside this range and will show as unmapped in run_report.json,
# disproportionately affecting the Nontrafficking class. Unlike BRCA2's rat
# ortholog, 6A70's PKD1 chain uses native human UniProt numbering directly.
set -euo pipefail
cd /home/jv2807/protein_vis

sbatch slurm/run_visualize.sh \
    /scratch/jv2807/pkd1/data/pkd1_variants_labeled.csv \
    pdb:6A70:B \
    P98161 \
    configs/domains/P98161.yaml \
    /scratch/jv2807/pkd1/structure_cache \
    /scratch/jv2807/pkd1/results/protein_vis_pkd1 \
    pkd1
