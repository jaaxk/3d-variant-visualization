"""End-to-end offline rendering pipeline.

Reads a variant CSV, a structure previously cached by `protein-vis fetch`,
and a domain config; validates and aligns everything; and renders one
interactive HTML + static PNG pair per domain that has >=1 assigned variant,
plus two always-on whole-structure renders: a plain overview (variants
colored by class) and a domain overview (backbone colored by domain, all
domains at once, variants still colored by class) -- see
render.render_domain_overview_html/_png. If the structure file has more
than one distinct chain sequence (e.g. a multimeric complex), a third chain
overview (backbone colored by chain, identical-sequence chains grouped
under one color) is also rendered -- see
render.render_chain_overview_html/_png. Never touches the network -- that
boundary is enforced by structure.load_structure / load_uniprot_sequence
raising if the required fetch hasn't happened yet.
"""

from __future__ import annotations

import json
from pathlib import Path

from . import domains as domains_mod
from . import render as render_mod
from . import structure as structure_mod
from . import variants as variants_mod
from .colors import ColorMap, generate_categorical_palette
from .domains import Domain


def run_render(
    *,
    variants_csv: str | Path,
    structure_spec: str,
    uniprot_accession: str,
    cache_dir: str | Path,
    domains_arg: str,
    output_dir: str | Path,
    strict_wt: bool = True,
    min_identity: float = 0.4,
    chain_labels: dict[str, str] | None = None,
) -> Path:
    cache_dir = Path(cache_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/6] loading variants from {variants_csv}")
    long_df = variants_mod.load_variant_table(variants_csv)
    print(f"      {len(long_df)} variant(s) across classes: "
          f"{sorted(long_df['class_name'].unique())}")

    print(f"[2/6] loading UniProt reference sequence for {uniprot_accession}")
    reference_seq = structure_mod.load_uniprot_sequence(cache_dir, uniprot_accession)
    print(f"      reference length: {len(reference_seq)} aa")

    print(f"[3/6] validating variants against reference sequence (strict={strict_wt})")
    long_df, messages = variants_mod.validate_against_sequence(
        long_df, reference_seq, strict=strict_wt
    )
    for msg in messages:
        print(f"      {msg}")
    print(f"      {len(long_df)} variant(s) remain after validation")

    print(f"[4/6] loading structure {structure_spec!r} and aligning to reference")
    struct = structure_mod.load_structure(structure_spec, cache_dir)
    alignment = structure_mod.align_to_reference(struct, reference_seq, min_identity=min_identity)
    print(f"      structure chain {struct.chain_id!r}, {len(struct.resnums)} residues; "
          f"alignment identity={alignment.identity:.1%} coverage={alignment.coverage:.1%}")

    print(f"[5/6] loading domain config ({domains_arg})")
    if domains_arg == "auto":
        uniprot_json = cache_dir / "uniprot" / f"{uniprot_accession}.json"
        domains_list = domains_mod.auto_domains_from_uniprot(uniprot_json)
    else:
        domains_list = domains_mod.load_domain_config(domains_arg)
    print(f"      {len(domains_list)} domain(s): {[d.name for d in domains_list]}")

    groups = domains_mod.group_variants_by_domain(long_df, domains_list)
    print(f"      variants land in {len(groups)} domain(s) with >=1 hit: {sorted(groups)}")

    domains_by_name = {d.name: d for d in domains_list}

    print(f"[6/6] rendering {2 + len(groups)} visualization(s) to {output_dir}")
    colors = ColorMap()
    domain_colors = ColorMap(fallback_cycle=generate_categorical_palette(len(domains_list)))
    for domain in domains_list:
        domain_colors.get(domain.name)  # pre-seed so the legend follows domain_list's N->C order
    report: dict = {
        "structure_spec": structure_spec,
        "structure_chain": struct.chain_id,
        "uniprot_accession": uniprot_accession,
        "alignment_identity": alignment.identity,
        "alignment_coverage": alignment.coverage,
        "total_variants": len(long_df),
        "class_counts": long_df["class_name"].value_counts().to_dict(),
        "domains_rendered": {},
    }

    def _render_one(name: str, sub_df, title: str, domain: Domain | None = None) -> None:
        highlight_resnums = domains_mod.resnums_for_domain(domain, alignment) if domain else None
        html_path = output_dir / f"{name}.html"
        png_path = output_dir / f"{name}.png"
        render_mod.render_interactive_html(
            struct, sub_df, alignment, colors, html_path, title=title, cache_dir=cache_dir,
            highlight_resnums=highlight_resnums,
        )
        render_mod.render_static_png(
            struct, sub_df, alignment, colors, png_path, title=title,
            highlight_resnums=highlight_resnums,
        )
        mapped, n_unmapped = render_mod._variant_positions_with_coords(sub_df, struct, alignment)
        report["domains_rendered"][name] = {
            "n_variants": len(sub_df),
            "n_mapped": len(mapped),
            "n_unmapped": n_unmapped,
            "n_structure_residues_highlighted": len(highlight_resnums) if highlight_resnums else None,
        }
        print(f"      wrote {html_path.name} / {png_path.name} "
              f"({len(mapped)} mapped, {n_unmapped} unmapped)")

    _render_one("overview", long_df, f"{uniprot_accession} — whole structure overview")

    domain_overview_html = output_dir / "domain_overview.html"
    domain_overview_png = output_dir / "domain_overview.png"
    render_mod.render_domain_overview_html(
        struct, long_df, alignment, colors, domain_colors, domains_list, domain_overview_html,
        title=f"{uniprot_accession} — domain architecture", cache_dir=cache_dir,
    )
    render_mod.render_domain_overview_png(
        struct, long_df, alignment, colors, domain_colors, domains_list, domain_overview_png,
        title=f"{uniprot_accession} — domain architecture",
    )
    mapped, n_unmapped = render_mod._variant_positions_with_coords(long_df, struct, alignment)
    report["domain_overview"] = {
        "n_variants": len(long_df),
        "n_mapped": len(mapped),
        "n_unmapped": n_unmapped,
        "n_domains_colored": len(domains_list),
    }
    print(f"      wrote {domain_overview_html.name} / {domain_overview_png.name} "
          f"({len(mapped)} mapped, {n_unmapped} unmapped, {len(domains_list)} domain(s) colored)")

    for name, sub_df in groups.items():
        _render_one(name, sub_df, f"{uniprot_accession} — domain: {name}", domain=domains_by_name[name])

    chain_groups = render_mod._group_chains_by_sequence(struct, chain_labels)
    if len(chain_groups) > 1:
        chain_colors = ColorMap()
        chain_overview_html = output_dir / "chain_overview.html"
        chain_overview_png = output_dir / "chain_overview.png"
        render_mod.render_chain_overview_html(
            struct, long_df, alignment, colors, chain_colors, chain_overview_html,
            title=f"{uniprot_accession} — chains", cache_dir=cache_dir, chain_labels=chain_labels,
        )
        render_mod.render_chain_overview_png(
            struct, long_df, alignment, colors, chain_colors, chain_overview_png,
            title=f"{uniprot_accession} — chains", chain_labels=chain_labels,
        )
        mapped, n_unmapped = render_mod._variant_positions_with_coords(long_df, struct, alignment)
        report["chain_overview"] = {
            "n_variants": len(long_df),
            "n_mapped": len(mapped),
            "n_unmapped": n_unmapped,
            "chains_colored": list(chain_groups),
        }
        print(f"      wrote {chain_overview_html.name} / {chain_overview_png.name} "
              f"({len(mapped)} mapped, {n_unmapped} unmapped, "
              f"{len(chain_groups)} chain group(s) colored)")

    report_path = output_dir / "run_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(f"done. wrote {report_path}")
    return output_dir
