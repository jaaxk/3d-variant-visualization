# protein_vis

Interactive + static 3D visualization of protein variants, colored by class
and split per structural domain, built for reproducible use on NYU Torch HPC.

Given (1) a CSV of variants grouped into classes and (2) a 3D protein
structure, this produces:
- One **interactive, self-contained HTML** file per structural domain that
  has variants, plus one whole-structure overview -- rotatable/zoomable in
  any browser, no server or live Python kernel required.
- A matching lightweight static PNG per visualization (quick previews,
  README thumbnails).
- A `run_report.json` recording variant counts, alignment quality, and any
  variants that couldn't be mapped onto the structure (never silently dropped).

## Why the `fetch` / `render` split

Torch compute nodes may lack outbound internet access, so the CLI is split:

- `protein-vis fetch` -- **network**. Downloads and caches a structure file,
  a reference sequence, and (optionally) 3Dmol.js. Run this on the **login node**.
- `protein-vis render` -- **offline**. Reads only from the local cache
  populated by `fetch`. This is what runs inside the SLURM job. If something
  wasn't fetched yet, `render` fails with an explicit message telling you
  what to fetch, rather than trying (and failing) to reach the network itself.

## Input format

**Variants CSV** -- "wide" format: one column per class, values are missense
variants in `[WT][pos][MUT]` notation (e.g. `L2484P`), position numbered
against a canonical reference sequence (e.g. UniProt). Columns may be ragged
(different lengths) -- short columns are just padded with blank cells.

```csv
pathogenic,benign,hypomorphic
L2484P,S2500A,K2551I
D2566V,,N2553D
```

If you only have a single unlabeled column (e.g. a plain variant list), use
`scripts/prepare_variant_csv.py` to attach a class name:

```bash
python scripts/prepare_variant_csv.py \
  --input variants.csv --class-name hypomorphic --output variants_labeled.csv
```

**Structure** -- specified as one of:
- `alphafold:<UniProt accession>` -- fetched from AlphaFold DB. Numbering
  matches the reference sequence exactly, so this is the preferred default
  *when a model exists*. Not every protein has one (notably: BRCA2 does not
  -- see below).
- `pdb:<PDB ID>` or `pdb:<PDB ID>:<CHAIN>` -- fetched from the RCSB PDB.
- `file:<local path>` -- your own PDB/mmCIF file (e.g. a ColabFold prediction).

A structure's own residue numbering is **never** assumed to match the
reference sequence's numbering. `protein-vis` always performs a real pairwise
sequence alignment (BLOSUM62, Biopython) to build that correspondence, and
reports identity/coverage in the output -- this matters even for "normal"
cases (crystal structures routinely have gaps/offsets), not just unusual ones.

**Domains** -- a YAML config of named residue ranges (see
`configs/domains/P51587.yaml` for BRCA2), or `--domains auto` to derive a
coarser split from UniProt's own Domain/Region/Repeat feature annotations
when no curated config exists for your protein. A variant is rendered in
every domain visualization whose range contains it (domains may overlap by
design, e.g. a sub-domain nested inside a larger one).

## Usage

```bash
# One-time, per machine: create + provision the Singularity overlay (see
# "Container setup" below, or the `environments` skill for the general
# pattern). Fetch is network-dependent -- do it on the login node, or via a
# quick srun allocation:

protein-vis fetch \
  --structure pdb:1IYJ --uniprot P51587 \
  --cache-dir /scratch/$USER/my_project/structure_cache --bootstrap-js

# Rendering (offline, no network) always goes through SLURM. There's one
# generic job script -- slurm/run_visualize.sh -- that takes all
# configuration as positional args:
#
#   sbatch slurm/run_visualize.sh <variants_csv> <structure_spec> \
#       <uniprot_accession> <domains_config> <cache_dir> <output_dir> [job_label]
#
# Don't call it with raw paths from the shell though -- write a thin
# `slurm/submit_<run_name>.sh` wrapper (see the existing submit_*.sh files)
# that documents exactly how that run was invoked, for reproducibility.
```

## Example runs

Four runs exist so far, each with its own `submit_*.sh` wrapper:

**BRCA2, panel-C hypomorphic only** (`slurm/submit_brca2_hypomorphic.sh`).
BRCA2 (UniProt P51587, 3418 aa) has no AlphaFold DB model (excluded from the
bulk release, likely due to size) and no human PDB structure of its
DNA-binding domain. The only structural coverage of that region is **PDB
1IYJ** (Yang et al. 2002, *Science*) -- the **rat** BRCA2 DBD-DSS1-ssDNA
complex, ~78% identity to human. `protein-vis` maps human variant positions
onto it via the real sequence alignment described above (not raw
residue-number equality), and reports the resulting identity/coverage in
`run_report.json` and in each HTML's header so nobody mistakes it for a
human structure. The 135 panel-C variants (positions 2484-3180) span five
curated sub-domains of the DBD (`Helical`, `OB1`, `Tower`, `OB2`, `OB3` --
see `configs/domains/P51587.yaml`), so this run produces one HTML/PNG pair
per sub-domain plus a whole-structure overview.

**BRCA2, 3-class (benign / pathogenic / hypomorphic)**
(`slurm/submit_brca2_3class.sh`). Same structure/domains as above, but the
input CSV merges the panel-C hypomorphic list with `6k_dms_dataset.xlsx`'s
graded ACMG classifications (`scripts/build_brca2_3class_variants.py`):
`Benign */Pathogenic *` collapse to `benign`/`pathogenic`, `Uncertain` rows
are dropped, protein-level variants with conflicting benign/pathogenic
calls across duplicate genomic-level entries are dropped and logged, and
**hypomorphic takes precedence** over benign/pathogenic on overlap.

**BRCA2, 3-class, PDB 1MIU** (`slurm/submit_brca2_3class_1miu.sh`). Same
3980-variant dataset and domain config as the 1IYJ 3-class run above, but
against **PDB 1MIU chain A** (Yang et al. 2002, *Science*) instead -- the
same crystallographic study's **mouse** BRCA2(2378-3115)-DSS1 complex,
rather than 1IYJ's **rat** construct. Alignment identity/coverage is
similar (~78.2% identity, 19.6% coverage vs. 1IYJ's 78.1%/17.2%), but 1MIU's
longer resolved span yields more mapped variants (3598 vs. 3141 of 3980).
Results are written to a separate directory
(`/scratch/jv2807/sounak_brca2/results/protein_vis_panelC_3class_1MIU`) so
they can be compared side-by-side with the 1IYJ run rather than overwriting
it.

**PKD1 (Polycystin-1)** (`slurm/submit_pkd1.sh`). PKD1 (UniProt P98161, 4303
aa) also has no AlphaFold model. Structure is **PDB 6A70** (Su et al. 2018,
*Science*) -- the human PKD1-PKD2 complex cryo-EM structure -- but its PKD1
chain only covers residues 3049-4169 (the transmembrane/pore-forming
region), so ~42% of variants (mostly in the large N-terminal extracellular
domain) show as unmapped; this disproportionately affects the
`Nontrafficking` class. Unlike 1IYJ, 6A70's PKD1 chain uses native human
UniProt numbering directly. Classes come straight from the `Pathogenicity
mechanism` column (`Benign` / `Nontrafficking` / `Function`) via
`scripts/build_pkd1_variants.py`. Domains are `PLAT_domain`, `Five_TM_domain`,
`VGIC_pore_module` (see `configs/domains/P98161.yaml`).

## Container setup (one-time, manual)

See the `environments` skill for the general Torch HPC pattern (why
`apptainer overlay create` and not the gzipped `/share/apps/overlay-fs-ext3/`
templates, the Miniforge URL, the NFS stale-lock gotcha, etc). Concretely,
for this project:

```bash
mkdir -p /scratch/$USER/protein_vis_singularity
apptainer overlay create --size 5000 --sparse --create-dir ext3 \
  /scratch/$USER/protein_vis_singularity/protein_vis.ext3

singularity exec --overlay protein_vis.ext3:rw \
  /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
  /bin/bash -c '
    export PATH=/ext3/miniforge3/bin:$PATH
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /ext3/miniforge.sh
    bash /ext3/miniforge.sh -b -p /ext3/miniforge3
    source /ext3/miniforge3/etc/profile.d/conda.sh
    conda create -y -p /ext3/miniforge3/envs/pv python=3.11
    conda activate pv
    pip install -e /home/$USER/protein_vis
    cat > /ext3/env.sh << "EOF"
#!/bin/bash
source /ext3/miniforge3/etc/profile.d/conda.sh
conda activate pv
export PYTHONNOUSERSITE=1
EOF
  '
```

No GPU is required for this workload (pandas/biopython/py3Dmol/matplotlib
are all CPU-only, lightweight); `slurm/run_visualize.sh` accordingly
requests no `--gres`/`--nv`. All of the above -- including package installs
and dataset inspection, not just rendering -- should go through `srun`, not
run directly on the login node (only `protein-vis fetch`'s actual network
calls are a login-node exception, being genuinely negligible).

## Repository layout

```
src/protein_vis/
  cli.py         click CLI: `fetch` (network) / `render` (offline)
  variants.py    variant CSV parsing + validation against a reference sequence
  structure.py   structure fetch/cache/load + sequence alignment
  domains.py     domain config loading + variant-to-domain assignment
  colors.py      class -> color mapping
  render.py      interactive HTML (py3Dmol, offline-JS-inlined) + static PNG
  pipeline.py    orchestration
configs/domains/  per-protein curated domain boundary YAMLs
scripts/          data-prep helpers -- prepare_variant_csv.py (single-column
                  -> labeled CSV), build_brca2_3class_variants.py, and
                  build_pkd1_variants.py (bespoke per-dataset reformatters;
                  see "Example runs" above)
slurm/            run_visualize.sh (generic, parameterized) + one
                  submit_<run_name>.sh wrapper per historical run
tests/            pytest, fully offline (no network/cluster required)
```

## Adding a new protein

1. `protein-vis fetch --structure <spec> --uniprot <accession> --cache-dir ...`
2. (Optional) curate `configs/domains/<accession>.yaml` if you want a
   finer-grained domain split than `--domains auto` (UniProt's own
   Domain/Region/Repeat features) provides.
3. Format your variant CSV (wide, one column per class) and run `render`.

No other code changes needed -- this is the intended extension point.

## Testing

```bash
pip install -e ".[dev]"
pytest tests/
```

All tests run fully offline against small synthetic fixtures in
`tests/fixtures/` -- no network or Singularity container required, though
running inside the container (with its installed dependencies) works too.
