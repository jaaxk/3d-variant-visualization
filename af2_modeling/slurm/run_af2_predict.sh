#!/bin/bash
#SBATCH --account=torch_pr_800_cds
#SBATCH --partition=h200_cds
#SBATCH --gres=gpu:h200:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --chdir=/home/jv2807/dms_side_projects/protein_vis/af2_modeling
#SBATCH --output=/home/jv2807/dms_side_projects/protein_vis/af2_modeling/slurm/logs/%j.out
#SBATCH --job-name=af2_predict
#SBATCH --mail-user=jv2807@nyu.edu
#SBATCH --mail-type=ALL

# Generic AlphaFold2-multimer job -- calls the cluster's official,
# pre-provisioned AlphaFold2 2.3.2 install directly
# (/share/apps/af2/2.3.2-20231225/run-alphafold.py). That script itself
# handles the Singularity/database-mounting chain internally (see
# af2_modeling/README.md for how that chain works) -- we never touch
# Singularity ourselves, just call it like a normal CLI from within a
# GPU-allocated job. TF_FORCE_UNIFIED_MEMORY/XLA_PYTHON_CLIENT_MEM_FRACTION
# (letting JAX address beyond the GPU's own VRAM) are already set inside
# that script -- nothing to add here.
#
# Wrapped in watch_af2_progress.py so a live elapsed-time + stage-change
# log streams to stdout instead of a silent multi-hour/day black box --
# this is genuinely untested territory (7,207 residues is well beyond any
# published AF2-multimer benchmark), so visibility into whether it's
# converging matters as much as the result itself.
#
# Usage:
#   sbatch slurm/run_af2_predict.sh <fasta_path> <output_dir> [job_label] \
#       [num_predictions_per_model] [models_to_relax] [max_template_date]
#
# Don't call sbatch directly with raw paths -- use a thin submit_*.sh
# wrapper (see submit_pkd1_pkd2_full.sh) that documents exactly how each
# run was invoked, for reproducibility.

set -euo pipefail

FASTA_PATH="$1"
OUTPUT_DIR="$2"
JOB_LABEL="${3:-af2_predict}"
NUM_PREDICTIONS_PER_MODEL="${4:-1}"
MODELS_TO_RELAX="${5:-best}"
MAX_TEMPLATE_DATE="${6:-$(date +%F)}"

exec > >(tee "/home/jv2807/dms_side_projects/protein_vis/af2_modeling/slurm/logs/${JOB_LABEL}.log") 2>&1

AF2_PYTHON="/share/apps/af2/pyenv/bin/python3"
AF2_SCRIPT="/share/apps/af2/2.3.2-20231225/run-alphafold.py"
WATCH_SCRIPT="/home/jv2807/dms_side_projects/protein_vis/af2_modeling/scripts/watch_af2_progress.py"

mkdir -p "${OUTPUT_DIR}"

echo "[${JOB_LABEL}] fasta=${FASTA_PATH}"
echo "[${JOB_LABEL}] output_dir=${OUTPUT_DIR}"
echo "[${JOB_LABEL}] num_predictions_per_model=${NUM_PREDICTIONS_PER_MODEL} models_to_relax=${MODELS_TO_RELAX} max_template_date=${MAX_TEMPLATE_DATE}"
echo "[${JOB_LABEL}] node=$(hostname) gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>&1 || echo unknown)"

"${AF2_PYTHON}" "${WATCH_SCRIPT}" \
    --log-file "/home/jv2807/dms_side_projects/protein_vis/af2_modeling/slurm/logs/${JOB_LABEL}_af2_raw.log" \
    -- \
    "${AF2_PYTHON}" "${AF2_SCRIPT}" \
    --fasta_paths="${FASTA_PATH}" \
    --output_dir="${OUTPUT_DIR}" \
    --max_template_date="${MAX_TEMPLATE_DATE}" \
    --db_preset=full_dbs \
    --model_preset=multimer \
    --num_multimer_predictions_per_model="${NUM_PREDICTIONS_PER_MODEL}" \
    --models_to_relax="${MODELS_TO_RELAX}" \
    --use_gpu=true \
    --enable_gpu_relax=true \
    || { echo "[${JOB_LABEL}] AlphaFold2 prediction FAILED"; exit 1; }

echo ""
echo "[${JOB_LABEL}] Done. Outputs in ${OUTPUT_DIR}"
