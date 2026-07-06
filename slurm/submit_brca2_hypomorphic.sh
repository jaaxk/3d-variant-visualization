#!/bin/bash
# Reproduces the original BRCA2 panel-C hypomorphic-only run.
# Structure: PDB 1IYJ chain B (rat BRCA2 DBD-DSS1-ssDNA complex, Yang et al.
# 2002 Science) -- alphafold:P51587 does NOT work, AlphaFold DB has no model
# for BRCA2 (excluded from the bulk release, likely due to size). Variant
# positions are mapped onto 1IYJ's numbering via a real sequence alignment
# (~78% identity to human), not raw residue-number equality -- see
# configs/domains/P51587.yaml and src/protein_vis/structure.py for detail.
set -euo pipefail
cd /home/jv2807/protein_vis

sbatch slurm/run_visualize.sh \
    /scratch/jv2807/sounak_brca2/dataset/brca2_panelC_hypomorphic.csv \
    pdb:1IYJ:B \
    P51587 \
    configs/domains/P51587.yaml \
    /scratch/jv2807/sounak_brca2/structure_cache \
    /scratch/jv2807/sounak_brca2/results/protein_vis_panelC \
    brca2_panelC_hypomorphic
