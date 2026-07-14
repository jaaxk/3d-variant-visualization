#!/bin/bash
# PKD1 variant visualization against the newest structure: the ~5k-residue
# AlphaFold Server PKD1-PKD2 complex prediction with PDB 6A70's real
# cryo-EM coordinates grafted onto the overlapping region (see
# af2_modeling/scripts/graft_6a70_onto_prediction.py). Same 4-chain layout
# as submit_pkd1_pkd2_complex.sh's structure (A/B/C = 3x PKD2, D = PKD1
# fragment, internally renumbered -- real UniProt positions recovered via
# protein_vis's own sequence alignment, not assumed from the file).
#
# This is also the first run to exercise overview.html's new "Color by:"
# toggle (Chain / EM-AF / Topology / Domain -- see
# render.render_multi_mode_overview_html), replacing the old separate
# domain_overview.html/chain_overview.html:
#
# --domains auto / --domains-for "Q13563=auto": Domain mode shows real
#   UniProt domain names (PKD1's LRRNT/LRR repeats/PKD 1-17/REJ/GAIN-B/
#   PLAT/GPS, PKD2's EF-hand, etc.), not just the coarse 3-domain curated
#   P98161.yaml.
# --chain-uniprot "A=Q13563": lets Topology/Domain modes align PKD2's own
#   chain-group to its own UniProt reference (Q13563) instead of only
#   ever coloring the primary PKD1 chain.
# --provenance .../pkd1_pkd2_grafted_model0.provenance.json: adds the
#   EM/AF mode (which residues are 6A70's real deposited coordinates vs.
#   AlphaFold-predicted) -- computed by graft_6a70_onto_prediction.py's
#   --provenance-only mode directly against the already-grafted structure
#   (no re-graft needed, see that script's docstring).
# --interface-json configs/domains/pkd1_pkd2_interface.json: merges the
#   real, contact-computed PKD1-PKD2 "Interface" domain (see
#   scripts/compute_pkd1_pkd2_interface.py) into the Domain mode for both
#   accessions.
#
# chain_labels / class_colors / variant_class_overrides: identical scheme to
# submit_pkd1_pkd2_complex.sh, for visual consistency across the two runs.
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis

sbatch slurm/run_visualize.sh \
    /scratch/jv2807/pkd1/data/pkd1_variants_labeled.csv \
    file:/scratch/jv2807/pkd1/af2_predictions/6a70_af2_server_5k/pkd1_pkd2_grafted_model0.pdb \
    P98161 \
    auto \
    /scratch/jv2807/pkd1/structure_cache \
    /scratch/jv2807/pkd1/results/protein_vis_pkd1_pkd2_grafted \
    pkd1_pkd2_grafted \
    "D=PKD1,A=PKD2" \
    "Function=#E53935,Nontrafficking=#FB8C00,Benign=#43A047,Temperature_recovered=#1E88E5" \
    "R2215W=Temperature_recovered,R2220W=Temperature_recovered" \
    "A=Q13563" \
    "Q13563=auto" \
    /scratch/jv2807/pkd1/af2_predictions/6a70_af2_server_5k/pkd1_pkd2_grafted_model0.provenance.json \
    configs/domains/pkd1_pkd2_interface.json
