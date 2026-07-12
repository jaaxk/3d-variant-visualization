#!/bin/bash
# Same PKD1 labeled-variant dataset as submit_pkd1.sh / submit_pkd1_alphafold.sh,
# but against a newer AlphaFold Server (AF3) prediction of the PKD1-PKD2
# complex: the last ~2100 residues of PKD1 (trimmed) plus 3 full-length
# PKD2 chains (the known 1:3 PKD1:PKD2 stoichiometry).
#
# Structure: file:.../fold_pkd1_2_complex_human_trimmed_model_0.cif (top-
# ranked of 5 models -- this time a clear winner, ranking_score=0.71 vs.
# 0.68/0.69/0.63/0.62 for the others, unlike the tied full-length PKD1-only
# run). The file has 4 chains: A/B/C are the 3 PKD2 copies (968 residues
# each, native numbering, starting at Met1), D is the trimmed PKD1 fragment
# (2083 residues, renumbered 1-2083 internally -- NOT native UniProt
# numbering). No explicit chain is passed in the structure spec (protein_vis's
# file: source kind doesn't support a :chain suffix the way pdb: does) --
# structure.load_structure's auto-select picks the chain with the most
# standard residues, which is chain D here (2083 > 968), so this resolves
# correctly without one. Verified via srun before this run (see
# run_report.json's "structure_chain" field, added specifically so this is
# auditable for multi-chain files like this one).
#
# Because PKD1 is trimmed, some variants (positions before the fragment's
# start) will be unmapped -- protein_vis's real BLOSUM62 sequence alignment
# (structure.align_to_reference) finds where this fragment's sequence
# matches within the full P98161 reference regardless of the structure's
# own internal renumbering, so the exact trim boundary never needs to be
# known ahead of time; unmapped variants are reported in run_report.json,
# never silently dropped.
#
# Domains: --domains auto, consistent with submit_pkd1_alphafold.sh.
#
# chain_labels: "D=PKD1,A=PKD2" -- labels the chain_overview legend with the
# real protein names. This information isn't in the structure file itself
# (an AlphaFold Server mmCIF's _entity.pdbx_description is empty for every
# entity), so it's supplied here from what we know built this complex.
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis

sbatch slurm/run_visualize.sh \
    /scratch/jv2807/pkd1/data/pkd1_variants_labeled.csv \
    file:/scratch/jv2807/pkd1/data/fold_pkd1_2_complex_human_trimmed/fold_pkd1_2_complex_human_trimmed_model_0.cif \
    P98161 \
    auto \
    /scratch/jv2807/pkd1/structure_cache \
    /scratch/jv2807/pkd1/results/protein_vis_pkd1_pkd2_complex \
    pkd1_pkd2_complex \
    "D=PKD1,A=PKD2"
