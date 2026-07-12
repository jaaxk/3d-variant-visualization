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
#       <uniprot_accession> <domains_config> <cache_dir> <output_dir> \
#       [job_label] [chain_labels] [class_colors] [variant_class_overrides]
#
# chain_labels (optional) -- comma-separated chain_id=Name pairs (e.g.
# "D=PKD1,A=PKD2"), only used to label the chain_overview legend for
# multi-chain structures with real protein names instead of raw chain
# letters. Naming any one chain in a group of identical-sequence chains
# labels the whole group.
#
# class_colors (optional) -- comma-separated ClassName=#hex pairs (e.g.
# "Function=#E53935,Benign=#43A047") to override this run's variant-class
# color scheme. Run-specific only -- colors.py's defaults are untouched.
#
# variant_class_overrides (optional) -- comma-separated VariantName=ClassName
# pairs (e.g. "R2215W=Temperature_recovered,R2220W=Temperature_recovered") to
# reassign specific variants to a (possibly synthetic) class before
# rendering, so that class -- and its color, via class_colors above --
# takes precedence over whatever class the variant CSV originally assigned.
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
CHAIN_LABELS="${8:-}"
CLASS_COLORS="${9:-}"
VARIANT_CLASS_OVERRIDES="${10:-}"

exec > >(tee "/home/jv2807/dms_side_projects/protein_vis/slurm/logs/${JOB_LABEL}.log") 2>&1

SINGULARITY_OVERLAY="/scratch/jv2807/protein_vis_singularity/protein_vis.ext3"
SINGULARITY_IMAGE="/share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif"
SIF_CPU="singularity exec --overlay ${SINGULARITY_OVERLAY}:ro ${SINGULARITY_IMAGE} /bin/bash -c"

mkdir -p "${OUTPUT_DIR}"

echo "[${JOB_LABEL}] variants=${VARIANTS_CSV}"
echo "[${JOB_LABEL}] structure=${STRUCTURE_SPEC} uniprot=${UNIPROT_ACCESSION}"
echo "[${JOB_LABEL}] domains=${DOMAINS_CONFIG}"
echo "[${JOB_LABEL}] output_dir=${OUTPUT_DIR}"

CHAIN_LABEL_FLAGS=""
if [ -n "${CHAIN_LABELS}" ]; then
    echo "[${JOB_LABEL}] chain_labels=${CHAIN_LABELS}"
    IFS=',' read -ra PAIRS <<< "${CHAIN_LABELS}"
    for pair in "${PAIRS[@]}"; do
        CHAIN_LABEL_FLAGS="${CHAIN_LABEL_FLAGS} --chain-label ${pair}"
    done
fi

CLASS_COLOR_FLAGS=""
if [ -n "${CLASS_COLORS}" ]; then
    echo "[${JOB_LABEL}] class_colors=${CLASS_COLORS}"
    IFS=',' read -ra PAIRS <<< "${CLASS_COLORS}"
    for pair in "${PAIRS[@]}"; do
        CLASS_COLOR_FLAGS="${CLASS_COLOR_FLAGS} --class-color ${pair}"
    done
fi

VARIANT_CLASS_FLAGS=""
if [ -n "${VARIANT_CLASS_OVERRIDES}" ]; then
    echo "[${JOB_LABEL}] variant_class_overrides=${VARIANT_CLASS_OVERRIDES}"
    IFS=',' read -ra PAIRS <<< "${VARIANT_CLASS_OVERRIDES}"
    for pair in "${PAIRS[@]}"; do
        VARIANT_CLASS_FLAGS="${VARIANT_CLASS_FLAGS} --variant-class ${pair}"
    done
fi

$SIF_CPU "source /ext3/env.sh && cd /home/jv2807/dms_side_projects/protein_vis && python -m protein_vis.cli render \
    --variants ${VARIANTS_CSV} \
    --structure ${STRUCTURE_SPEC} \
    --uniprot ${UNIPROT_ACCESSION} \
    --cache-dir ${CACHE_DIR} \
    --domains ${DOMAINS_CONFIG} \
    --output-dir ${OUTPUT_DIR}${CHAIN_LABEL_FLAGS}${CLASS_COLOR_FLAGS}${VARIANT_CLASS_FLAGS}" \
    || { echo "[${JOB_LABEL}] Render FAILED"; exit 1; }

echo ""
echo "[${JOB_LABEL}] Done. Outputs (overview + per-domain .html/.png, run_report.json) in ${OUTPUT_DIR}"
