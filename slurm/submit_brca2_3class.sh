#!/bin/bash
# BRCA2 panel-C hypomorphic variants merged with 6k_dms_dataset.xlsx's
# benign/pathogenic classifications (hypomorphic takes precedence on
# overlap) -- see scripts/build_brca2_3class_variants.py for how the input
# CSV is built. Same structure/domain config as the hypomorphic-only run
# (PDB 1IYJ chain B, configs/domains/P51587.yaml).
set -euo pipefail
cd /home/jv2807/protein_vis

sbatch slurm/run_visualize.sh \
    /scratch/jv2807/sounak_brca2/dataset/brca2_panelC_and_6k_3class.csv \
    pdb:1IYJ:B \
    P51587 \
    configs/domains/P51587.yaml \
    /scratch/jv2807/sounak_brca2/structure_cache \
    /scratch/jv2807/sounak_brca2/results/protein_vis_panelC_3class \
    brca2_panelC_3class
