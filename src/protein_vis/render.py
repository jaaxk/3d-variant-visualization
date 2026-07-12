"""Rendering: interactive self-contained HTML (py3Dmol) + static PNG (matplotlib).

py3Dmol's `view.write_html()` does NOT emit a plain `<script src="...">` CDN
tag -- verified by reading its source (py3Dmol 2.5.5): it loads 3Dmol.js via
a JS `loadScriptAsync(uri)` call, gated behind a `$3Dmolpromise` variable
that's only assigned `if(typeof $3Dmolpromise === 'undefined')`. So to make
the output fully self-contained/offline:
  1. Construct the view with `js=""` so no CDN URL string appears anywhere
     in the output at all.
  2. Prepend a `<script>` block that (a) inlines the real 3Dmol.min.js
     content and (b) pre-defines `$3Dmolpromise = Promise.resolve()` --
     since script tags execute in document order, py3Dmol's own generated
     `<script>` block then sees `$3Dmolpromise` already defined and skips
     the (now-empty, never-called) loadScriptAsync entirely.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import py3Dmol  # noqa: E402

from .colors import ColorMap
from .domains import Domain, resnums_for_domain
from .structure import AlignmentResult, StructureData


class RenderError(RuntimeError):
    pass


def _build_legend_html(items: list[tuple[str, str]], heading: str | None = None) -> str:
    rows = "\n".join(
        f'<li><span style="display:inline-block;width:12px;height:12px;'
        f'background:{color};margin-right:6px;border-radius:50%;"></span>{name}</li>'
        for name, color in items
    )
    heading_html = (
        f'<h4 style="margin:8px 8px 0;font-family:sans-serif;">{heading}</h4>' if heading else ""
    )
    return (
        f"{heading_html}"
        f'<ul style="list-style:none;padding:0;margin:8px;font-family:sans-serif;">{rows}</ul>'
    )


def _variant_positions_with_coords(
    variants_df: pd.DataFrame, struct: StructureData, alignment: AlignmentResult
) -> tuple[list[dict], int]:
    """Map each variant row to a 3D coordinate via the alignment's position map.

    Returns (mapped_rows, n_unmapped). Unmapped variants (outside the aligned
    region, or in a gap with no resolved CA atom) are counted, never silently
    dropped without a trace.
    """
    mapped = []
    n_unmapped = 0
    for _, row in variants_df.iterrows():
        resnum = alignment.pos_to_resnum.get(int(row["pos"]))
        coord = struct.ca_coords.get(resnum) if resnum is not None else None
        if coord is None:
            n_unmapped += 1
            continue
        mapped.append({"class_name": row["class_name"], "raw": row["raw"], "coord": coord})
    return mapped, n_unmapped


def render_static_png(
    struct: StructureData,
    variants_df: pd.DataFrame,
    alignment: AlignmentResult,
    colors: ColorMap,
    out_path: str | Path,
    *,
    title: str,
    highlight_resnums: set[int] | None = None,
    zoom_padding: float = 6.0,
) -> Path:
    """Render a static 3D backbone + variant scatter.

    highlight_resnums -- structure resnums belonging to the domain being
    visualized (if any). When given, that backbone segment is drawn in a
    distinct color on top of the (dimmed) full backbone, and the camera is
    cropped to its bounding box -- without this, every domain's PNG shows
    the same whole-structure view and is indistinguishable from the others
    except for which dots happen to be colored.
    """
    mapped, n_unmapped = _variant_positions_with_coords(variants_df, struct, alignment)

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    sorted_resnums = sorted(struct.resnums)
    base_color = "#DDE3E8" if highlight_resnums else "#B0BEC5"
    base_alpha = 0.5 if highlight_resnums else 0.6

    backbone = np.array([struct.ca_coords[r] for r in sorted_resnums if r in struct.ca_coords])
    if len(backbone) > 0:
        ax.plot(*backbone.T, "-", lw=0.6, color=base_color, alpha=base_alpha)

    highlight_coords = None
    if highlight_resnums:
        highlight_coords = np.array(
            [struct.ca_coords[r] for r in sorted_resnums if r in highlight_resnums and r in struct.ca_coords]
        )
        if len(highlight_coords) > 0:
            ax.plot(*highlight_coords.T, "-", lw=2.4, color="#4A7FBF", alpha=0.9)

    by_class: dict[str, list[np.ndarray]] = {}
    for item in mapped:
        by_class.setdefault(item["class_name"], []).append(item["coord"])

    for class_name, coords in by_class.items():
        arr = np.array(coords)
        ax.scatter(*arr.T, color=colors.get(class_name), s=40, edgecolors="white", label=class_name)

    if highlight_coords is not None and len(highlight_coords) > 0:
        mins = highlight_coords.min(axis=0) - zoom_padding
        maxs = highlight_coords.max(axis=0) + zoom_padding
        ax.set_xlim(mins[0], maxs[0])
        ax.set_ylim(mins[1], maxs[1])
        ax.set_zlim(mins[2], maxs[2])

    footnote = f"{len(mapped)} variant(s) plotted"
    if n_unmapped:
        footnote += f", {n_unmapped} unmapped (outside aligned/resolved structure region)"
    ax.set_title(f"{title}\n{footnote}", fontsize=10)
    if by_class:
        ax.legend(loc="upper left", fontsize=8)
    ax.set_axis_off()

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_interactive_html(
    struct: StructureData,
    variants_df: pd.DataFrame,
    alignment: AlignmentResult,
    colors: ColorMap,
    out_path: str | Path,
    *,
    title: str,
    cache_dir: str | Path,
    highlight_resnums: set[int] | None = None,
) -> Path:
    """Render an interactive, self-contained HTML viewer.

    highlight_resnums -- structure resnums belonging to the domain being
    visualized (if any). When given, that residue range is colored
    distinctly and the initial camera zooms to it instead of the whole
    structure -- without this, every domain's HTML shows the same
    whole-structure view and is indistinguishable from the others except
    for which spheres happen to be colored.
    """
    js_path = Path(cache_dir) / "js" / "3Dmol.min.js"
    if not js_path.exists():
        raise RenderError(
            f"no cached 3Dmol.min.js at {js_path} -- run "
            f"`protein-vis fetch --bootstrap-js` on the login node first"
        )
    js_text = js_path.read_text()

    mapped, n_unmapped = _variant_positions_with_coords(variants_df, struct, alignment)

    view = py3Dmol.view(width=900, height=650, js="")
    view.addModel(struct.raw_text, struct.fmt)
    view.setStyle({}, {"cartoon": {"color": "lightgray"}})
    if highlight_resnums:
        highlight_sel = {"chain": struct.chain_id, "resi": sorted(highlight_resnums)}
        view.setStyle(highlight_sel, {"cartoon": {"color": "steelblue"}})
    for item in mapped:
        x, y, z = (float(c) for c in item["coord"])
        view.addSphere(
            {
                "center": {"x": x, "y": y, "z": z},
                "radius": 1.2,
                "color": colors.get(item["class_name"]),
            }
        )
    if highlight_resnums:
        view.zoomTo(highlight_sel)
    else:
        view.zoomTo()
    viewer_html = view.write_html()

    legend_html = _build_legend_html(colors.legend_items())
    footnote = f"{len(mapped)} variant(s) shown"
    if n_unmapped:
        footnote += f", {n_unmapped} unmapped (outside aligned/resolved structure region)"

    full_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
</head>
<body style="margin:0;padding:0;font-family:sans-serif;">
<h2 style="margin:8px;">{title}</h2>
<p style="margin:8px;color:#555;font-size:13px;">
  Structure: {struct.chain_id} | alignment identity {alignment.identity:.1%},
  coverage {alignment.coverage:.1%} | {footnote}
</p>
{legend_html}
<script>{js_text}</script>
<script>var $3Dmolpromise = Promise.resolve();</script>
{viewer_html}
</body>
</html>"""

    if "cdn.jsdelivr" in full_html or "3dmol.org" in full_html.lower():
        raise RenderError(
            "generated HTML unexpectedly references an external CDN -- "
            "self-contained/offline guarantee violated"
        )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(full_html)
    return out_path


def _domain_resnums(
    domains: list[Domain], alignment: AlignmentResult, struct: StructureData
) -> list[tuple[Domain, list[int]]]:
    """(domain, sorted structure resnums) for every domain with >=1 resolved
    residue, in the given domain order -- shared by the PNG/HTML domain
    overview renderers so a domain's resnum set is computed once."""
    out = []
    for domain in domains:
        resnums = sorted(resnums_for_domain(domain, alignment) & set(struct.ca_coords))
        if resnums:
            out.append((domain, resnums))
    return out


def render_domain_overview_png(
    struct: StructureData,
    variants_df: pd.DataFrame,
    alignment: AlignmentResult,
    class_colors: ColorMap,
    domain_colors: ColorMap,
    domains: list[Domain],
    out_path: str | Path,
    *,
    title: str,
) -> Path:
    """Whole-structure static render with the backbone colored by domain --
    every domain in `domains` drawn simultaneously in its own color (later
    domains win on any overlapping residues), instead of the single
    highlighted/dimmed domain used by render_static_png's per-domain mode.
    Variants are still overlaid and colored by class exactly as every other
    render in this pipeline; only the backbone coloring differs. No
    zoom/crop -- this is meant to show the whole domain architecture at once.
    """
    mapped, n_unmapped = _variant_positions_with_coords(variants_df, struct, alignment)
    domain_segments = _domain_resnums(domains, alignment, struct)

    fig = plt.figure(figsize=(9, 7.5))
    ax = fig.add_subplot(111, projection="3d")

    backbone = np.array([struct.ca_coords[r] for r in sorted(struct.resnums) if r in struct.ca_coords])
    if len(backbone) > 0:
        ax.plot(*backbone.T, "-", lw=0.6, color="#DDE3E8", alpha=0.5)

    for domain, resnums in domain_segments:
        coords = np.array([struct.ca_coords[r] for r in resnums])
        ax.plot(*coords.T, "-", lw=2.0, color=domain_colors.get(domain.name), alpha=0.9)

    by_class: dict[str, list[np.ndarray]] = {}
    for item in mapped:
        by_class.setdefault(item["class_name"], []).append(item["coord"])
    for class_name, coords in by_class.items():
        arr = np.array(coords)
        ax.scatter(
            *arr.T, color=class_colors.get(class_name), s=40, edgecolors="white", label=class_name
        )

    footnote = f"{len(mapped)} variant(s) plotted"
    if n_unmapped:
        footnote += f", {n_unmapped} unmapped (outside aligned/resolved structure region)"
    ax.set_title(f"{title}\n{footnote}", fontsize=10)

    class_legend = None
    if by_class:
        class_legend = ax.legend(
            loc="upper left", fontsize=8, title="Variant class", title_fontsize=8
        )
    if domain_segments:
        domain_handles = [
            Line2D([0], [0], color=domain_colors.get(d.name), lw=3) for d, _ in domain_segments
        ]
        domain_labels = [d.name for d, _ in domain_segments]
        ax.legend(
            domain_handles,
            domain_labels,
            loc="upper right",
            fontsize=6,
            title="Domain",
            title_fontsize=7,
            ncol=2 if len(domain_handles) > 12 else 1,
        )
        if class_legend is not None:
            ax.add_artist(class_legend)
    ax.set_axis_off()

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_domain_overview_html(
    struct: StructureData,
    variants_df: pd.DataFrame,
    alignment: AlignmentResult,
    class_colors: ColorMap,
    domain_colors: ColorMap,
    domains: list[Domain],
    out_path: str | Path,
    *,
    title: str,
    cache_dir: str | Path,
) -> Path:
    """Whole-structure interactive render with the backbone colored by
    domain -- every domain in `domains` colored simultaneously (later
    domains win on any overlapping residues), instead of a single
    highlighted/zoomed domain. Variants are still colored by class exactly
    as every other render. No zoom/crop -- shows the whole structure."""
    js_path = Path(cache_dir) / "js" / "3Dmol.min.js"
    if not js_path.exists():
        raise RenderError(
            f"no cached 3Dmol.min.js at {js_path} -- run "
            f"`protein-vis fetch --bootstrap-js` on the login node first"
        )
    js_text = js_path.read_text()

    mapped, n_unmapped = _variant_positions_with_coords(variants_df, struct, alignment)
    domain_segments = _domain_resnums(domains, alignment, struct)

    view = py3Dmol.view(width=900, height=650, js="")
    view.addModel(struct.raw_text, struct.fmt)
    view.setStyle({}, {"cartoon": {"color": "lightgray"}})

    domain_legend_items: list[tuple[str, str]] = []
    for domain, resnums in domain_segments:
        color = domain_colors.get(domain.name)
        view.setStyle({"chain": struct.chain_id, "resi": resnums}, {"cartoon": {"color": color}})
        domain_legend_items.append((domain.name, color))

    for item in mapped:
        x, y, z = (float(c) for c in item["coord"])
        view.addSphere(
            {
                "center": {"x": x, "y": y, "z": z},
                "radius": 1.2,
                "color": class_colors.get(item["class_name"]),
            }
        )
    view.zoomTo()
    viewer_html = view.write_html()

    domain_legend_html = _build_legend_html(domain_legend_items, heading="Domains")
    class_legend_html = _build_legend_html(class_colors.legend_items(), heading="Variant class")
    footnote = f"{len(mapped)} variant(s) shown"
    if n_unmapped:
        footnote += f", {n_unmapped} unmapped (outside aligned/resolved structure region)"

    full_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
</head>
<body style="margin:0;padding:0;font-family:sans-serif;">
<h2 style="margin:8px;">{title}</h2>
<p style="margin:8px;color:#555;font-size:13px;">
  Structure: {struct.chain_id} | alignment identity {alignment.identity:.1%},
  coverage {alignment.coverage:.1%} | {footnote}
</p>
<div style="display:flex;flex-wrap:wrap;gap:24px;">
{domain_legend_html}
{class_legend_html}
</div>
<script>{js_text}</script>
<script>var $3Dmolpromise = Promise.resolve();</script>
{viewer_html}
</body>
</html>"""

    if "cdn.jsdelivr" in full_html or "3dmol.org" in full_html.lower():
        raise RenderError(
            "generated HTML unexpectedly references an external CDN -- "
            "self-contained/offline guarantee violated"
        )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(full_html)
    return out_path


def _group_chains_by_sequence(struct: StructureData) -> dict[str, list[str]]:
    """Chain ids grouped by identical sequence, e.g. the 3 PKD2 copies in a
    PKD1-PKD2 complex collapse to one group/color/legend entry instead of
    three -- "colored by chain" should mean by distinct molecule, not by
    every individual chain letter."""
    groups: dict[str, list[str]] = {}
    for chain_id in struct.all_chain_ca_coords:
        seq = struct.all_chain_sequences.get(chain_id, "")
        groups.setdefault(seq, []).append(chain_id)
    return {", ".join(chain_ids): chain_ids for chain_ids in groups.values()}


def render_chain_overview_png(
    struct: StructureData,
    variants_df: pd.DataFrame,
    alignment: AlignmentResult,
    class_colors: ColorMap,
    chain_colors: ColorMap,
    out_path: str | Path,
    *,
    title: str,
) -> Path:
    """Whole-structure static render with the backbone colored by chain
    (e.g. PKD1 vs. PKD2 in a multimeric complex) instead of by domain.
    Meant to be generated only when the structure file has more than one
    chain. Variants still overlaid and colored by class as everywhere else.
    """
    mapped, n_unmapped = _variant_positions_with_coords(variants_df, struct, alignment)
    groups = _group_chains_by_sequence(struct)

    fig = plt.figure(figsize=(9, 7.5))
    ax = fig.add_subplot(111, projection="3d")

    for label, chain_ids in groups.items():
        color = chain_colors.get(label)
        for chain_id in chain_ids:
            coords_by_resnum = struct.all_chain_ca_coords[chain_id]
            coords = np.array([coords_by_resnum[r] for r in sorted(coords_by_resnum)])
            ax.plot(*coords.T, "-", lw=1.4, color=color, alpha=0.85)

    by_class: dict[str, list[np.ndarray]] = {}
    for item in mapped:
        by_class.setdefault(item["class_name"], []).append(item["coord"])
    for class_name, coords in by_class.items():
        arr = np.array(coords)
        ax.scatter(
            *arr.T, color=class_colors.get(class_name), s=40, edgecolors="white", label=class_name
        )

    footnote = f"{len(mapped)} variant(s) plotted"
    if n_unmapped:
        footnote += f", {n_unmapped} unmapped (outside aligned/resolved structure region)"
    ax.set_title(f"{title}\n{footnote}", fontsize=10)

    class_legend = None
    if by_class:
        class_legend = ax.legend(
            loc="upper left", fontsize=8, title="Variant class", title_fontsize=8
        )
    chain_handles = [Line2D([0], [0], color=chain_colors.get(label), lw=3) for label in groups]
    ax.legend(chain_handles, list(groups), loc="upper right", fontsize=8, title="Chain", title_fontsize=8)
    if class_legend is not None:
        ax.add_artist(class_legend)
    ax.set_axis_off()

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_chain_overview_html(
    struct: StructureData,
    variants_df: pd.DataFrame,
    alignment: AlignmentResult,
    class_colors: ColorMap,
    chain_colors: ColorMap,
    out_path: str | Path,
    *,
    title: str,
    cache_dir: str | Path,
) -> Path:
    """Whole-structure interactive render with the backbone colored by
    chain. Meant to be generated only when the structure file has more than
    one chain. Variants still colored by class as everywhere else."""
    js_path = Path(cache_dir) / "js" / "3Dmol.min.js"
    if not js_path.exists():
        raise RenderError(
            f"no cached 3Dmol.min.js at {js_path} -- run "
            f"`protein-vis fetch --bootstrap-js` on the login node first"
        )
    js_text = js_path.read_text()

    mapped, n_unmapped = _variant_positions_with_coords(variants_df, struct, alignment)
    groups = _group_chains_by_sequence(struct)

    view = py3Dmol.view(width=900, height=650, js="")
    view.addModel(struct.raw_text, struct.fmt)
    view.setStyle({}, {"cartoon": {"color": "lightgray"}})

    chain_legend_items: list[tuple[str, str]] = []
    for label, chain_ids in groups.items():
        color = chain_colors.get(label)
        for chain_id in chain_ids:
            view.setStyle({"chain": chain_id}, {"cartoon": {"color": color}})
        chain_legend_items.append((label, color))

    for item in mapped:
        x, y, z = (float(c) for c in item["coord"])
        view.addSphere(
            {
                "center": {"x": x, "y": y, "z": z},
                "radius": 1.2,
                "color": class_colors.get(item["class_name"]),
            }
        )
    view.zoomTo()
    viewer_html = view.write_html()

    chain_legend_html = _build_legend_html(chain_legend_items, heading="Chain")
    class_legend_html = _build_legend_html(class_colors.legend_items(), heading="Variant class")
    footnote = f"{len(mapped)} variant(s) shown"
    if n_unmapped:
        footnote += f", {n_unmapped} unmapped (outside aligned/resolved structure region)"

    full_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
</head>
<body style="margin:0;padding:0;font-family:sans-serif;">
<h2 style="margin:8px;">{title}</h2>
<p style="margin:8px;color:#555;font-size:13px;">
  Structure: {struct.chain_id} | alignment identity {alignment.identity:.1%},
  coverage {alignment.coverage:.1%} | {footnote}
</p>
<div style="display:flex;flex-wrap:wrap;gap:24px;">
{chain_legend_html}
{class_legend_html}
</div>
<script>{js_text}</script>
<script>var $3Dmolpromise = Promise.resolve();</script>
{viewer_html}
</body>
</html>"""

    if "cdn.jsdelivr" in full_html or "3dmol.org" in full_html.lower():
        raise RenderError(
            "generated HTML unexpectedly references an external CDN -- "
            "self-contained/offline guarantee violated"
        )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(full_html)
    return out_path
