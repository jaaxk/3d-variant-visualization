#!/bin/bash
#SBATCH --account=torch_pr_800_cds
#SBATCH --time=0:20:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --chdir=/home/jv2807/protein_vis
#SBATCH --output=/home/jv2807/protein_vis/slurm/logs/%j.out
#SBATCH --job-name=protein_vis_brca2_panelC
#SBATCH --mail-user=jv2807@nyu.edu
#SBATCH --mail-type=ALL

exec > >(tee /home/jv2807/protein_vis/slurm/logs/protein_vis_brca2_panelC.log) 2>&1

# Render-only step -- NO network calls happen here. `protein-vis fetch` MUST
# have already been run interactively on the LOGIN NODE before this is
# submitted (see README "BRCA2 test run" section), populating
# STRUCTURE_CACHE with the structure file, UniProt JSON/FASTA, and
# 3Dmol.min.js.
#
# Structure source: PDB 1IYJ chain B (rat BRCA2 DBD-DSS1-ssDNA complex,
# Yang et al. 2002 Science) -- alphafold:P51587 does NOT work, AlphaFold DB
# has no model for BRCA2 (verified: the /api/prediction/P51587 endpoint
# returns no entries, likely excluded from the bulk release due to size).
# Variant positions are mapped onto 1IYJ's numbering via a real sequence
# alignment (~79% identity to human), not raw residue-number equality --
# see configs/domains/P51587.yaml and src/protein_vis/structure.py for detail.

VARIANTS_CSV="/scratch/jv2807/sounak_brca2/dataset/brca2_panelC_hypomorphic.csv"
STRUCTURE_CACHE="/scratch/jv2807/sounak_brca2/structure_cache"
DOMAINS_CONFIG="/home/jv2807/protein_vis/configs/domains/P51587.yaml"
OUTPUT_DIR="/scratch/jv2807/sounak_brca2/results/protein_vis_panelC"

SINGULARITY_OVERLAY="/scratch/jv2807/protein_vis_singularity/protein_vis.ext3"
SINGULARITY_IMAGE="/share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif"
SIF_CPU="singularity exec --overlay ${SINGULARITY_OVERLAY}:ro ${SINGULARITY_IMAGE} /bin/bash -c"

mkdir -p "${OUTPUT_DIR}"

echo "Rendering BRCA2 panelC hypomorphic variants -> ${OUTPUT_DIR}"
$SIF_CPU "source /ext3/env.sh && cd /home/jv2807/protein_vis && python -m protein_vis.cli render \
    --variants ${VARIANTS_CSV} \
    --structure pdb:1IYJ:B \
    --uniprot P51587 \
    --cache-dir ${STRUCTURE_CACHE} \
    --domains ${DOMAINS_CONFIG} \
    --output-dir ${OUTPUT_DIR}" \
    || { echo "Render FAILED"; exit 1; }

echo ""
echo "Done. Outputs (overview + per-domain .html/.png, run_report.json) in ${OUTPUT_DIR}"
