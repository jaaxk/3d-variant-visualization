#!/bin/bash
# PKD1 variant visualization against the full ~7,207-residue PKD1/PKD2
# complex: the already-grafted PKD1-PKD2 6A70 structure
# (pkd1_pkd2_grafted_model0.pdb) with a standalone full-length PKD1 monomer
# AlphaFold Server prediction (fold_pkd1_human) grafted on top, filling in
# PKD1's N-terminal ectodomain (~residues 1-2220) that neither 6A70 nor the
# original complex prediction covered (see
# af2_modeling/scripts/graft_6a70_onto_prediction.py's stage-2 usage in its
# docstring, and pkd1_full_pkd2_grafted.pdb/.provenance.json).
#
# Chain layout of the new structure (see stage-2 graft's console output):
# A = PKD1 (full-length, 4303 aa -- kept the monomer prediction's own id),
# B/C/D = the 3 PKD2 chains carried through unchanged from the complex
# (originally A/B/C there, renamed on id collision with the incoming PKD1
# chain -- see --include-unmatched-template-chains).
#
# --provenance .../pkd1_full_pkd2_grafted.provenance.json now drives the
# 3-way "EM/AF" mode: 6A70 / AlphaFold2: Complex / AlphaFold2: PKD1 monomer
# (see src/protein_vis/colors.py PROVENANCE_COLORS and pipeline.py's
# provenance block, both label-count-agnostic).
#
# --interface-json configs/domains/pkd1_pkd2_interface.json: reused as-is --
# a static, structure-independent artifact (see
# scripts/compute_pkd1_pkd2_interface.py), no recomputation needed.
#
# chain_labels / class_colors / variant_class_overrides: identical scheme to
# submit_pkd1_pkd2_grafted.sh, for visual consistency across runs.
#
# confidence_enabled ("1"): adds a "Confidence" mode -- per-residue pLDDT
# from each CA atom's B-factor (both grafted-in AlphaFold Server predictions
# write pLDDT there), bucketed into AlphaFold's own confidence bands.
# Residues the provenance JSON labels "6A70" show as "Experimentally
# resolved" instead of a confidence band (see colors.py CONFIDENCE_COLORS,
# pipeline.py's confidence block).
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis

sbatch slurm/run_visualize.sh \
    /scratch/jv2807/pkd1/data/pkd1_variants_labeled.csv \
    file:/scratch/jv2807/pkd1/af2_predictions/pkd1_full_pkd2_grafted/pkd1_full_pkd2_grafted.pdb \
    P98161 \
    auto \
    /scratch/jv2807/pkd1/structure_cache \
    /scratch/jv2807/pkd1/results/protein_vis_pkd1_full_pkd2_grafted \
    pkd1_full_pkd2_grafted \
    "A=PKD1,B=PKD2" \
    "Function=#E53935,Nontrafficking=#FB8C00,Benign=#43A047,Temperature_recovered=#1E88E5" \
    "R2215W=Temperature_recovered,R2220W=Temperature_recovered" \
    "B=Q13563" \
    "Q13563=auto" \
    /scratch/jv2807/pkd1/af2_predictions/pkd1_full_pkd2_grafted/pkd1_full_pkd2_grafted.provenance.json \
    configs/domains/pkd1_pkd2_interface.json \
    1
