"""End-to-end offline rendering pipeline.

Reads a variant CSV, a structure previously cached by `protein-vis fetch`,
and a domain config; validates and aligns everything; and renders one
interactive HTML + static PNG pair per domain that has >=1 assigned variant,
plus an always-on whole-structure "overview.html"/"overview_labeled.html" --
a single interactive page with a "Color by:" dropdown that switches the
backbone coloring between several precomputed schemes (Chain / Topology /
Domain, plus EM/AF when a --provenance file is given) without reloading --
see render.render_multi_mode_overview_html / render.ModeScheme. This
replaces the older separate domain_overview.html/chain_overview.html (those
two are now just entries in the same dropdown). "overview.png" stays a
plain static snapshot (variants only, no backbone coloring toggle -- a PNG
can't be interactive). Every interactive HTML is rendered twice -- once
plain and once as a "<name>_labeled.html" twin with each variant's name
floated next to its point. Never touches the network -- that boundary is
enforced by structure.load_structure / load_uniprot_sequence raising if the
required fetch hasn't happened yet.

class_color_overrides / variant_class_overrides let a caller retint a
specific run's variant classes (e.g. "Function=#E53935") and/or reassign
specific variants to a synthetic class that takes visual precedence (e.g.
{"R2215W": "Temperature_recovered"}) -- run-specific, passed at invocation
time (CLI --class-color / --variant-class), not hardcoded in colors.py, so
the tool stays reusable across proteins/runs with their own default colors.

chain_uniprot / domains_for generalize the Domain/Topology coloring modes to
multi-protein complexes: by default domains/topology are only known for the
single primary chain (`uniprot_accession`, aligned via `alignment` below).
chain_uniprot ("chain_id=ACCESSION") names which UniProt accession every
other chain-group's representative chain belongs to; domains_for
("ACCESSION=path_or_auto") says how to load *that* accession's domains
(defaults to "auto" if not given). Both are optional/additive -- a
single-chain run (BRCA2, PKD1-only) passes neither and behaves exactly as
before.
"""

from __future__ import annotations

import json
from pathlib import Path

from . import domains as domains_mod
from . import render as render_mod
from . import structure as structure_mod
from . import variants as variants_mod
from .colors import ColorMap, PROVENANCE_COLORS, TOPOLOGY_COLORS, generate_categorical_palette
from .domains import Domain
from .structure import StructureData


def _load_domains(accession: str, domains_arg: str, cache_dir: Path) -> list[Domain]:
    if domains_arg == "auto":
        uniprot_json = cache_dir / "uniprot" / f"{accession}.json"
        return domains_mod.auto_domains_from_uniprot(uniprot_json)
    return domains_mod.load_domain_config(domains_arg)


def _load_topology(accession: str, cache_dir: Path) -> list[Domain]:
    return domains_mod.auto_topology_from_uniprot(cache_dir / "uniprot" / f"{accession}.json")


def _chain_group_for(chain_groups: dict[str, list[str]], chain_id: str) -> list[str]:
    """Every chain id sharing chain_id's identical-sequence group -- a
    domain/topology mapping computed from one representative chain applies
    identically to every chain in its group (same sequence => same internal
    numbering, e.g. all 3 PKD2 copies)."""
    for chain_ids in chain_groups.values():
        if chain_id in chain_ids:
            return chain_ids
    return [chain_id]


def _sub_structure(struct: StructureData, chain_id: str) -> StructureData:
    """A throwaway StructureData standing in for one secondary chain, just
    enough for structure.align_to_reference (which only reads .sequence and
    .resnums) -- lets every chain-group reuse the exact same alignment
    machinery as the primary chain instead of a separate code path."""
    ca_coords = struct.all_chain_ca_coords[chain_id]
    return StructureData(
        chain_id=chain_id,
        resnums=sorted(ca_coords),
        sequence=struct.all_chain_sequences[chain_id],
        ca_coords=ca_coords,
        raw_text="",
        fmt=struct.fmt,
    )


def _regions_and_legend(
    domains_list: list[Domain],
    alignment: structure_mod.AlignmentResult,
    resolved_resnums: set[int],
    paint_chain_ids: list[str],
    color_map: ColorMap,
) -> tuple[list[tuple[str, list[int], str]], list[tuple[str, str]]]:
    """(regions, legend_items) for one chain-group's domain/topology
    coloring -- regions repeated across every chain id in `paint_chain_ids`
    (the chain-group, see _chain_group_for) since they all share the same
    resnum mapping."""
    regions: list[tuple[str, list[int], str]] = []
    legend: list[tuple[str, str]] = []
    seen_names: set[str] = set()
    for domain, resnums in render_mod._domain_resnums(domains_list, alignment, resolved_resnums):
        color = color_map.get(domain.name)
        for chain_id in paint_chain_ids:
            regions.append((chain_id, resnums, color))
        if domain.name not in seen_names:
            seen_names.add(domain.name)
            legend.append((domain.name, color))
    return regions, legend


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
    class_color_overrides: dict[str, str] | None = None,
    variant_class_overrides: dict[str, str] | None = None,
    chain_uniprot: dict[str, str] | None = None,
    domains_for: dict[str, str] | None = None,
    provenance_path: str | Path | None = None,
    interface_json: str | Path | None = None,
    interface_domain_name: str = "Interface",
) -> Path:
    cache_dir = Path(cache_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/6] loading variants from {variants_csv}")
    long_df = variants_mod.load_variant_table(variants_csv)
    if variant_class_overrides:
        for raw, new_class in variant_class_overrides.items():
            n = (long_df["raw"] == raw).sum()
            long_df.loc[long_df["raw"] == raw, "class_name"] = new_class
            print(f"      variant override: {raw} -> class {new_class!r} ({n} row(s))")
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
    domains_list = _load_domains(uniprot_accession, domains_arg, cache_dir)
    print(f"      {len(domains_list)} domain(s): {[d.name for d in domains_list]}")

    interface_by_accession: dict[str, list[int]] = {}
    if interface_json:
        interface_by_accession = json.loads(Path(interface_json).read_text())
        positions = interface_by_accession.get(uniprot_accession)
        if positions:
            domains_list.append(
                Domain(name=interface_domain_name, positions=frozenset(positions),
                       source=f"computed contact interface ({interface_json})")
            )
            print(f"      + {interface_domain_name} ({len(positions)} residue(s) from {interface_json})")

    groups = domains_mod.group_variants_by_domain(long_df, domains_list)
    print(f"      variants land in {len(groups)} domain(s) with >=1 hit: {sorted(groups)}")

    domains_by_name = {d.name: d for d in domains_list}

    print(f"[6/6] rendering {2 + len(groups)} visualization(s) to {output_dir}")
    colors = ColorMap(overrides=class_color_overrides)
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
        labeled_html_path = output_dir / f"{name}_labeled.html"
        png_path = output_dir / f"{name}.png"
        render_mod.render_interactive_html(
            struct, sub_df, alignment, colors, html_path, title=title, cache_dir=cache_dir,
            highlight_resnums=highlight_resnums,
        )
        render_mod.render_interactive_html(
            struct, sub_df, alignment, colors, labeled_html_path, title=title, cache_dir=cache_dir,
            highlight_resnums=highlight_resnums, show_variant_labels=True,
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
        print(f"      wrote {html_path.name} / {labeled_html_path.name} / {png_path.name} "
              f"({len(mapped)} mapped, {n_unmapped} unmapped)")

    # --- Chain-group bookkeeping shared by the Chain/Domain/Topology modes ---
    chain_groups = render_mod._group_chains_by_sequence(struct, chain_labels)
    accession_by_chain: dict[str, str] = {struct.chain_id: uniprot_accession, **(chain_uniprot or {})}
    domains_for = domains_for or {}

    # --- Chain mode ---
    chain_colors = ColorMap()
    chain_regions: list[tuple[str, list[int], str]] = []
    chain_legend: list[tuple[str, str]] = []
    for label, chain_ids in chain_groups.items():
        color = chain_colors.get(label)
        chain_legend.append((label, color))
        for chain_id in chain_ids:
            chain_regions.append((chain_id, sorted(struct.all_chain_ca_coords[chain_id]), color))

    # --- Domain + Topology modes (per chain-group, each aligned to its own accession) ---
    # First pass: resolve each chain-group's own alignment/domain-list/topology-list
    # so domain_colors' fallback palette can be sized to the TRUE total category
    # count across every accession (not just the primary one) before any color gets
    # assigned -- otherwise a secondary accession's domains could exhaust an
    # under-sized cycle and fall back to a generic gray.
    per_chain_domain_data: dict[str, tuple] = {}  # chain_id -> (alignment, resolved_resnums, chain_domains, topo_domains)
    for chain_id, accession in accession_by_chain.items():
        if accession == uniprot_accession and chain_id == struct.chain_id:
            chain_alignment, resolved_resnums, chain_domains = alignment, set(struct.ca_coords), domains_list
        else:
            sub_struct = _sub_structure(struct, chain_id)
            ref_seq = structure_mod.load_uniprot_sequence(cache_dir, accession)
            chain_alignment = structure_mod.align_to_reference(sub_struct, ref_seq, min_identity=min_identity)
            resolved_resnums = set(sub_struct.ca_coords)
            chain_domains = _load_domains(accession, domains_for.get(accession, "auto"), cache_dir)
            positions = interface_by_accession.get(accession)
            if positions:
                chain_domains.append(
                    Domain(name=interface_domain_name, positions=frozenset(positions),
                           source=f"computed contact interface ({interface_json})")
                )
        topo_domains = _load_topology(accession, cache_dir)
        per_chain_domain_data[chain_id] = (chain_alignment, resolved_resnums, chain_domains, topo_domains)

    total_domain_count = sum(len(v[2]) for v in per_chain_domain_data.values())
    domain_colors = ColorMap(fallback_cycle=generate_categorical_palette(total_domain_count))
    topology_colors = ColorMap(overrides=TOPOLOGY_COLORS)

    domain_regions: list[tuple[str, list[int], str]] = []
    domain_legend: list[tuple[str, str]] = []
    domain_seen: set[str] = set()
    topology_regions: list[tuple[str, list[int], str]] = []
    topology_legend: list[tuple[str, str]] = []
    topology_seen: set[str] = set()

    for chain_id in accession_by_chain:
        paint_chain_ids = _chain_group_for(chain_groups, chain_id)
        chain_alignment, resolved_resnums, chain_domains, topo_domains = per_chain_domain_data[chain_id]

        regions, legend = _regions_and_legend(
            chain_domains, chain_alignment, resolved_resnums, paint_chain_ids, domain_colors
        )
        domain_regions.extend(regions)
        for name, color in legend:
            if name not in domain_seen:
                domain_seen.add(name)
                domain_legend.append((name, color))

        regions, legend = _regions_and_legend(
            topo_domains, chain_alignment, resolved_resnums, paint_chain_ids, topology_colors
        )
        topology_regions.extend(regions)
        for name, color in legend:
            if name not in topology_seen:
                topology_seen.add(name)
                topology_legend.append((name, color))

    # --- EM/AF provenance mode (optional -- only when a grafted structure's
    # provenance sidecar was passed; see af2_modeling/scripts/graft_6a70_onto_prediction.py) ---
    provenance_regions: list[tuple[str, list[int], str]] = []
    provenance_legend: list[tuple[str, str]] = []
    if provenance_path:
        provenance_colors = ColorMap(overrides=PROVENANCE_COLORS)
        em_color = provenance_colors.get("6A70")
        af_color = provenance_colors.get("AlphaFold2")
        provenance_legend = provenance_colors.legend_items()
        provenance = json.loads(Path(provenance_path).read_text())
        for chain_id, ca_coords in struct.all_chain_ca_coords.items():
            em_resnums = sorted(set(provenance.get(chain_id, [])) & set(ca_coords))
            af_resnums = sorted(set(ca_coords) - set(em_resnums))
            if em_resnums:
                provenance_regions.append((chain_id, em_resnums, em_color))
            if af_resnums:
                provenance_regions.append((chain_id, af_resnums, af_color))

    modes = {
        "Chain": render_mod.ModeScheme(label="Chain", regions=chain_regions, legend_items=chain_legend),
        "Domain": render_mod.ModeScheme(label="Domain", regions=domain_regions, legend_items=domain_legend),
        "Topology": render_mod.ModeScheme(
            label="Topology", regions=topology_regions, legend_items=topology_legend
        ),
    }
    if provenance_path:
        modes["EM/AF"] = render_mod.ModeScheme(
            label="EM/AF", regions=provenance_regions, legend_items=provenance_legend
        )

    overview_html = output_dir / "overview.html"
    overview_labeled_html = output_dir / "overview_labeled.html"
    overview_png = output_dir / "overview.png"
    render_mod.render_multi_mode_overview_html(
        struct, long_df, alignment, colors, modes, overview_html,
        title=f"{uniprot_accession} — overview", cache_dir=cache_dir,
    )
    render_mod.render_multi_mode_overview_html(
        struct, long_df, alignment, colors, modes, overview_labeled_html,
        title=f"{uniprot_accession} — overview", cache_dir=cache_dir, show_variant_labels=True,
    )
    render_mod.render_static_png(
        struct, long_df, alignment, colors, overview_png,
        title=f"{uniprot_accession} — whole structure overview",
    )
    mapped, n_unmapped = render_mod._variant_positions_with_coords(long_df, struct, alignment)
    report["overview"] = {
        "n_variants": len(long_df),
        "n_mapped": len(mapped),
        "n_unmapped": n_unmapped,
        "modes": list(modes),
    }
    print(f"      wrote {overview_html.name} / {overview_labeled_html.name} / {overview_png.name} "
          f"({len(mapped)} mapped, {n_unmapped} unmapped, modes: {list(modes)})")

    for name, sub_df in groups.items():
        _render_one(name, sub_df, f"{uniprot_accession} — domain: {name}", domain=domains_by_name[name])

    report_path = output_dir / "run_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(f"done. wrote {report_path}")
    return output_dir
