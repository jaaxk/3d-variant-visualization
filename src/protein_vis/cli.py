"""CLI entrypoint: `fetch` (network, login node) and `render` (offline, compute node).

This split exists because Torch HPC compute nodes may lack outbound internet
access: `fetch` downloads and caches everything `render` needs, and `render`
never makes a network call (enforced by structure.load_structure /
load_uniprot_sequence raising a clear error if something wasn't fetched).
"""

from __future__ import annotations

from pathlib import Path

import click

from . import pipeline
from . import structure as structure_mod


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--structure", "structures", multiple=True,
    help="Repeatable. e.g. 'alphafold:P51587', 'pdb:1IYJ', 'pdb:1IYJ:B'.",
)
@click.option(
    "--uniprot", "accessions", multiple=True,
    help="Repeatable UniProt accession(s) to fetch a reference FASTA/JSON for.",
)
@click.option("--cache-dir", required=True, type=click.Path(), help="Local cache directory.")
@click.option("--bootstrap-js", is_flag=True, help="Also cache 3Dmol.min.js for offline HTML.")
@click.option("--force", is_flag=True, help="Re-fetch even if already cached.")
def fetch(structures, accessions, cache_dir, bootstrap_js, force):
    """Run on the LOGIN NODE. The only subcommand that touches the network."""
    cache_dir = Path(cache_dir)
    for spec in structures:
        path = structure_mod.fetch_structure(spec, cache_dir, force=force)
        click.echo(f"fetched {spec} -> {path}")
    for accession in accessions:
        path = structure_mod.fetch_uniprot_reference(accession, cache_dir)
        click.echo(f"fetched UniProt {accession} -> {path}")
    if bootstrap_js:
        path = structure_mod.fetch_3dmol_js(cache_dir)
        click.echo(f"fetched 3Dmol.js -> {path}")


@cli.command()
@click.option("--variants", "variants_csv", required=True, type=click.Path(exists=True))
@click.option("--structure", "structure_spec", required=True,
              help="Must match a spec already fetched into the cache manifest.")
@click.option("--uniprot", "uniprot_accession", required=True,
              help="UniProt accession the variant CSV positions are numbered against.")
@click.option("--cache-dir", required=True, type=click.Path(exists=True))
@click.option("--domains", "domains_arg", required=True,
              help="Path to a domain YAML config, or the literal 'auto'.")
@click.option("--output-dir", required=True, type=click.Path())
@click.option("--no-strict-wt", is_flag=True,
              help="Warn (instead of drop) on WT-residue mismatches vs. the reference sequence.")
@click.option("--min-identity", default=0.4, show_default=True,
              help="Minimum alignment identity required between structure and reference.")
@click.option(
    "--chain-label", "chain_label_pairs", multiple=True,
    help="Repeatable 'chain_id=Name' (e.g. 'D=PKD1'), only used to label the "
         "chain_overview legend for multi-chain structures -- structure files "
         "carry no protein names themselves, so this must come from the caller. "
         "Naming any one chain in a group of identical-sequence chains labels "
         "the whole group.",
)
@click.option(
    "--class-color", "class_color_pairs", multiple=True,
    help="Repeatable 'ClassName=#hex' (e.g. 'Benign=#43A047') to override this "
         "run's variant-class color scheme. Applies only to this invocation --"
         "colors.py's defaults are untouched, so other proteins/runs are unaffected.",
)
@click.option(
    "--variant-class", "variant_class_pairs", multiple=True,
    help="Repeatable 'VariantName=ClassName' (e.g. 'R2215W=Temperature_recovered') "
         "to reassign specific variants to a (possibly synthetic) class before "
         "rendering -- this class then takes precedence over whatever class the "
         "variant CSV originally assigned it, since the reassignment happens "
         "before any coloring/grouping logic runs. Combine with --class-color to "
         "give the synthetic class its own color.",
)
def render(variants_csv, structure_spec, uniprot_accession, cache_dir, domains_arg,
           output_dir, no_strict_wt, min_identity, chain_label_pairs, class_color_pairs,
           variant_class_pairs):
    """Run inside the SLURM job / compute node. Never calls the network."""
    chain_labels = dict(pair.split("=", 1) for pair in chain_label_pairs)
    class_color_overrides = dict(pair.split("=", 1) for pair in class_color_pairs)
    variant_class_overrides = dict(pair.split("=", 1) for pair in variant_class_pairs)
    pipeline.run_render(
        variants_csv=variants_csv,
        structure_spec=structure_spec,
        uniprot_accession=uniprot_accession,
        cache_dir=cache_dir,
        domains_arg=domains_arg,
        output_dir=output_dir,
        strict_wt=not no_strict_wt,
        min_identity=min_identity,
        chain_labels=chain_labels,
        class_color_overrides=class_color_overrides,
        variant_class_overrides=variant_class_overrides,
    )


if __name__ == "__main__":
    cli()
