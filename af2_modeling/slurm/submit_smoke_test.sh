#!/bin/bash
# Smoke test before committing a long GPU allocation to the full 7,207-residue
# attempt: run the exact same invocation chain (model_preset=multimer, same
# flags, same databases) against just PKD2 alone (Q13563, 968 aa -- already
# fetched, no new FASTA needed) to confirm the flags/paths/singularity chain
# actually work end-to-end on this cluster, and to observe real AF2 log
# output so watch_af2_progress.py's stage-detection markers can be
# calibrated/extended if needed.
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis/af2_modeling

sbatch slurm/run_af2_predict.sh \
    /scratch/jv2807/pkd1/structure_cache/uniprot/Q13563.fasta \
    /scratch/jv2807/pkd1/af2_predictions/smoke_test_pkd2 \
    smoke_test_pkd2 \
    1 \
    best
