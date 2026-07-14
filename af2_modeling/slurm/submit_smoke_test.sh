#!/bin/bash
# Smoke test before committing a long GPU allocation to the full 7,207-residue
# attempt: run the exact same invocation chain (model_preset=multimer, same
# flags, same databases) against just PKD2 alone (Q13563, 968 aa -- already
# fetched, no new FASTA needed) to confirm the flags/paths/singularity chain
# actually work end-to-end on this cluster, and to observe real AF2 log
# output so watch_af2_progress.py's stage-detection markers can be
# calibrated/extended if needed.
#
# First attempt (job 13627928, on H200) got through the full ~57min MSA
# search (jackhmmer/hhblits) fine, then crashed the instant GPU model
# inference started: "CUDA_ERROR_NOT_FOUND: named symbol not found" looking
# up a PTX kernel. This AF2 install's container is CUDA 11.4
# (cuda11.4.2-cudnn8.2.4-devel-ubuntu20.04.3.sif per run-alphafold-all.bash)
# -- H200 is Hopper (compute capability 9.0), which predates/exceeds what
# CUDA 11.4 can JIT-compile for. Retrying on A100 (Ampere, sm_80 -- squarely
# within CUDA 11.4's support) instead, and reusing the MSAs/features already
# cached to output_dir from that first run (use_precomputed_msas=true) so
# this retry only re-does the fast GPU step, not the slow CPU MSA search.
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis/af2_modeling

sbatch --partition=a100_cds --gres=gpu:a100:1 slurm/run_af2_predict.sh \
    /scratch/jv2807/pkd1/structure_cache/uniprot/Q13563.fasta \
    /scratch/jv2807/pkd1/af2_predictions/smoke_test_pkd2 \
    smoke_test_pkd2_a100 \
    1 \
    best \
    "$(date +%F)" \
    true
