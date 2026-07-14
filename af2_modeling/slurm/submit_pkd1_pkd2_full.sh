#!/bin/bash
# Feasibility experiment: fold the FULL, untrimmed PKD1/PKD2 complex --
# 1x PKD1 (P98161, 4303 aa) + 3x PKD2 (Q13563, 968 aa each) = 7,207 total
# residues -- using the cluster's official AlphaFold2 2.3.2 install.
#
# GPU choice: A100, not H200. The smoke test on H200 (job 13627928) crashed
# immediately at GPU model inference ("CUDA_ERROR_NOT_FOUND: named symbol
# not found" resolving a PTX kernel) -- this AF2 install's container is
# CUDA 11.4 (cuda11.4.2-cudnn8.2.4-devel-ubuntu20.04.3.sif), which predates
# Hopper (H200, compute capability 9.0) support. A100 is Ampere (sm_80),
# squarely within CUDA 11.4's support. Tradeoff: A100 has less VRAM than
# H200 (check the actual per-GPU memory on a100_cds before assuming it's
# enough) -- a real feasibility risk for 7,207 residues on top of the
# already-open question of whether this scale works at all; there's no
# newer-CUDA AF2 install on this cluster to fall back to instead.
#
# This is well beyond any published AlphaFold2/AlphaFold-Multimer benchmark
# (typically ~2-3k residues), so whether it finishes in a reasonable time
# (or fits in memory at all) is a genuine open question -- that's the point
# of this run. Watch the live progress via the job's .log file (tee'd to
# stdout by run_af2_predict.sh, wrapped in watch_af2_progress.py's
# elapsed-time + stage-change output). If it's clearly not converging (or
# OOMs), fall back to a partial-PKD1 run instead (rebuild the FASTA via
# build_complex_fasta.py --pkd1-start <N> with generous overlap into 6A70's
# ~2221-4303 resolved region, then resubmit through this same generic
# run_af2_predict.sh with the new FASTA).
#
# num_predictions_per_model=1 (not the default 5) to get a faster read on
# feasibility -- 5 models x 1 prediction each = 5 total predictions, not 25.
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis/af2_modeling

sbatch --partition=a100_cds --gres=gpu:a100:1 slurm/run_af2_predict.sh \
    fasta/pkd1_pkd2_complex_full.fasta \
    /scratch/jv2807/pkd1/af2_predictions/pkd1_pkd2_complex_full \
    pkd1_pkd2_complex_full \
    1 \
    best
