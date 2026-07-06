#!/bin/bash
#SBATCH --account=torch_pr_800_cds
#SBATCH --time=0:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --chdir=/home/jv2807/dms_side_projects/protein_vis
#SBATCH --output=/home/jv2807/dms_side_projects/protein_vis/slurm/logs/%j.out
#SBATCH --job-name=protein_vis_render
#SBATCH --mail-user=jv2807@nyu.edu
#SBATCH --mail-type=ALL

# Generic render job -- takes all configuration as positional args so the
# same script works for any protein/dataset. Renders variants against a
# structure that was already fetched (`protein-vis fetch`, login node) --
# NO network calls happen here.
#
# Usage:
#   sbatch slurm/run_visualize.sh <variants_csv> <structure_spec> \
#       <uniprot_accession> <domains_config> <cache_dir> <output_dir> [job_label]
#
# Don't call sbatch directly with raw paths from the shell -- use one of the
# thin per-run wrapper scripts (submit_*.sh) instead, which document exactly
# how each historical run was invoked for reproducibility.

set -euo pipefail

VARIANTS_CSV="$1"
STRUCTURE_SPEC="$2"
UNIPROT_ACCESSION="$3"
DOMAINS_CONFIG="$4"
CACHE_DIR="$5"
OUTPUT_DIR="$6"
JOB_LABEL="${7:-protein_vis_render}"

exec > >(tee "/home/jv2807/dms_side_projects/protein_vis/slurm/logs/${JOB_LABEL}.log") 2>&1

SINGULARITY_OVERLAY="/scratch/jv2807/protein_vis_singularity/protein_vis.ext3"
SINGULARITY_IMAGE="/share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif"
SIF_CPU="singularity exec --overlay ${SINGULARITY_OVERLAY}:ro ${SINGULARITY_IMAGE} /bin/bash -c"

mkdir -p "${OUTPUT_DIR}"

echo "[${JOB_LABEL}] variants=${VARIANTS_CSV}"
echo "[${JOB_LABEL}] structure=${STRUCTURE_SPEC} uniprot=${UNIPROT_ACCESSION}"
echo "[${JOB_LABEL}] domains=${DOMAINS_CONFIG}"
echo "[${JOB_LABEL}] output_dir=${OUTPUT_DIR}"

$SIF_CPU "source /ext3/env.sh && cd /home/jv2807/dms_side_projects/protein_vis && python -m protein_vis.cli render \
    --variants ${VARIANTS_CSV} \
    --structure ${STRUCTURE_SPEC} \
    --uniprot ${UNIPROT_ACCESSION} \
    --cache-dir ${CACHE_DIR} \
    --domains ${DOMAINS_CONFIG} \
    --output-dir ${OUTPUT_DIR}" \
    || { echo "[${JOB_LABEL}] Render FAILED"; exit 1; }

echo ""
echo "[${JOB_LABEL}] Done. Outputs (overview + per-domain .html/.png, run_report.json) in ${OUTPUT_DIR}"
