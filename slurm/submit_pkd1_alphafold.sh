#!/bin/bash
# Same PKD1 labeled-variant dataset as submit_pkd1.sh, but against a
# full-length AlphaFold Server (AF3) prediction the user generated directly
# -- not the bulk AlphaFold DB, which has no PKD1 model (excluded, likely
# due to size) -- instead of PDB 6A70's 3049-4169 fragment (~42% of variants
# unmapped there, concentrated in the Nontrafficking class).
#
# Structure: file:.../fold_pkd1_human_model_0.cif, the top-ranked of 5
# output models (all 5 tie at ranking_score=0.47/ptm=0.41, so "top-ranked by
# convention" is as good as any). It's an mmCIF file covering all 4303
# residues in native 1-based UniProt numbering (verified via
# structure.align_to_reference, not assumed) -- protein_vis already handles
# .cif via the same `file:<path>` source kind as any other structure
# (structure.load_structure picks MMCIFParser vs. PDBParser by extension),
# no format-specific code was needed.
#
# Domains: --domains auto (31 UniProt Domain/Region/Repeat features across
# the full protein) rather than the curated P98161.yaml, which only covers
# 6A70's fragment -- this also exercises the domain_overview visualization
# (whole structure colored by domain, all domains at once, with a legend)
# meaningfully for the first time.
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis

sbatch slurm/run_visualize.sh \
    /scratch/jv2807/pkd1/data/pkd1_variants_labeled.csv \
    file:/scratch/jv2807/pkd1/data/fold_pkd1_human/fold_pkd1_human_model_0.cif \
    P98161 \
    auto \
    /scratch/jv2807/pkd1/structure_cache \
    /scratch/jv2807/pkd1/results/protein_vis_pkd1_alphafold \
    pkd1_alphafold
