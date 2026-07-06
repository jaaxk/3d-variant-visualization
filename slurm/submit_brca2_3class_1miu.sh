#!/bin/bash
# Same BRCA2 panel-C hypomorphic + 6k_dms_dataset.xlsx 3-class variants as
# submit_brca2_3class.sh, but against a different human BRCA2 DBD structure:
# PDB 1MIU chain A (Yang et al. 2002 Science) instead of 1IYJ chain B.
#
# 1MIU is the BRCA2(2378-3115)-DSS1 complex expressed from MOUSE BRCA2
# (Mus musculus, chain A), whereas 1IYJ is the RAT ortholog (Rattus
# norvegicus) -- both are the same crystallographic study, different
# expression constructs/species, covering an overlapping but not identical
# span of the DNA-binding domain (1MIU: 2378-3115 vs. 1IYJ's construct).
# protein_vis maps human variant positions onto either via a real sequence
# alignment (never raw residue-number equality), and reports the resulting
# identity/coverage in run_report.json -- compare that against the 1IYJ run
# (78.1% identity / 17.2% coverage) to see which structure gives better
# coverage for this variant set.
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis

sbatch slurm/run_visualize.sh \
    /scratch/jv2807/sounak_brca2/dataset/brca2_panelC_and_6k_3class.csv \
    pdb:1MIU:A \
    P51587 \
    configs/domains/P51587.yaml \
    /scratch/jv2807/sounak_brca2/structure_cache \
    /scratch/jv2807/sounak_brca2/results/protein_vis_panelC_3class_1MIU \
    brca2_panelC_3class_1MIU
