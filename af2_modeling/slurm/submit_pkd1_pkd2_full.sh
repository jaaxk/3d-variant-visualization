#!/bin/bash
# Feasibility experiment: fold the FULL, untrimmed PKD1/PKD2 complex --
# 1x PKD1 (P98161, 4303 aa) + 3x PKD2 (Q13563, 968 aa each) = 7,207 total
# residues -- using the cluster's official AlphaFold2 2.3.2 install on a
# single H200 GPU (141GB, 2TB node RAM).
#
# This is well beyond any published AlphaFold2/AlphaFold-Multimer benchmark
# (typically ~2-3k residues), so whether it finishes in a reasonable time
# is a genuine open question -- that's the point of this run. Watch the
# live progress via the job's .log file (tee'd to stdout by
# run_af2_predict.sh, wrapped in watch_af2_progress.py's elapsed-time +
# stage-change output). If it's clearly not converging in a reasonable
# time, scancel it and fall back to a partial-PKD1 run instead (rebuild the
# FASTA via build_complex_fasta.py --pkd1-start <N> with generous overlap
# into 6A70's ~2221-4303 resolved region, then resubmit through this same
# generic run_af2_predict.sh with the new FASTA).
#
# num_predictions_per_model=1 (not the default 5) to get a faster read on
# feasibility -- 5 models x 1 prediction each = 5 total predictions, not 25.
set -euo pipefail
cd /home/jv2807/dms_side_projects/protein_vis/af2_modeling

sbatch slurm/run_af2_predict.sh \
    fasta/pkd1_pkd2_complex_full.fasta \
    /scratch/jv2807/pkd1/af2_predictions/pkd1_pkd2_complex_full \
    pkd1_pkd2_complex_full \
    1 \
    best
