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
@click.option(
    "--chain-uniprot", "chain_uniprot_pairs", multiple=True,
    help="Repeatable 'chain_id=ACCESSION' (e.g. 'A=Q13563'), only needed for "
         "multi-protein complexes -- names which UniProt accession each "
         "OTHER chain-group (besides the primary --uniprot/--structure chain) "
         "belongs to, so the overview's Domain/Topology color modes can be "
         "computed for it too via its own alignment. The primary chain "
         "already gets this from --uniprot and never needs an entry here.",
)
@click.option(
    "--domains-for", "domains_for_pairs", multiple=True,
    help="Repeatable 'ACCESSION=path_or_auto', paired with --chain-uniprot -- "
         "which domain config to use for a secondary accession's Domain mode "
         "(defaults to 'auto' if the accession has no entry here).",
)
@click.option(
    "--provenance", "provenance_path", type=click.Path(exists=True), default=None,
    help="Path to an EM/AF provenance JSON (chain_id -> {resnum_str: label}, see "
         "af2_modeling/scripts/graft_6a70_onto_prediction.py's --provenance-output) "
         "-- when given, the overview gets an extra 'EM/AF' color mode with one "
         "region per distinct label actually present (e.g. '6A70', 'AlphaFold2: "
         "Complex', 'AlphaFold2: PKD1 monomer' for a chained graft). Omit for "
         "structures that weren't produced by grafting a real structure onto a "
         "prediction.",
)
@click.option(
    "--interface-json", "interface_json", type=click.Path(exists=True), default=None,
    help="Path to a {ACCESSION: [uniprot_pos, ...]} JSON (see "
         "scripts/compute_pkd1_pkd2_interface.py) -- when given, an extra "
         "domain (named by --interface-domain-name) is merged into every "
         "accession's Domain mode wherever that accession has an entry.",
)
@click.option(
    "--interface-domain-name", default="Interface", show_default=True,
    help="Display name for the domain added by --interface-json.",
)
@click.option(
    "--confidence", "confidence_enabled", is_flag=True,
    help="Add a 'Confidence' color mode, bucketing each residue's CA B-factor "
         "into AlphaFold's own pLDDT bands (Very high/Confident/Low/Very low). "
         "Only meaningful for AlphaFold-derived structures -- a real "
         "crystallographic/EM structure's B-factor is NOT pLDDT, so don't pass "
         "this for e.g. a plain PDB-only run. If --provenance is also given, "
         "any residue it labels '6A70' is shown as 'Experimentally resolved' "
         "instead of a confidence band (it has no pLDDT to report).",
)
def render(variants_csv, structure_spec, uniprot_accession, cache_dir, domains_arg,
           output_dir, no_strict_wt, min_identity, chain_label_pairs, class_color_pairs,
           variant_class_pairs, chain_uniprot_pairs, domains_for_pairs, provenance_path,
           interface_json, interface_domain_name, confidence_enabled):
    """Run inside the SLURM job / compute node. Never calls the network."""
    chain_labels = dict(pair.split("=", 1) for pair in chain_label_pairs)
    class_color_overrides = dict(pair.split("=", 1) for pair in class_color_pairs)
    variant_class_overrides = dict(pair.split("=", 1) for pair in variant_class_pairs)
    chain_uniprot = dict(pair.split("=", 1) for pair in chain_uniprot_pairs)
    domains_for = dict(pair.split("=", 1) for pair in domains_for_pairs)
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
        chain_uniprot=chain_uniprot,
        domains_for=domains_for,
        provenance_path=provenance_path,
        interface_json=interface_json,
        interface_domain_name=interface_domain_name,
        confidence_enabled=confidence_enabled,
    )


if __name__ == "__main__":
    cli()
