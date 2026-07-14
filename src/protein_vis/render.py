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

import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
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


def _add_variant_label(view, item: dict) -> None:
    """Float the variant's name (e.g. 'R2215W') next to its sphere."""
    x, y, z = (float(c) for c in item["coord"])
    view.addLabel(
        item["raw"],
        {
            "position": {"x": x, "y": y, "z": z},
            "backgroundColor": "black",
            "backgroundOpacity": 0.6,
            "fontColor": "white",
            "fontSize": 11,
            "showBackground": True,
            "inFront": True,
        },
    )


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
    show_variant_labels: bool = False,
) -> Path:
    """Render an interactive, self-contained HTML viewer.

    highlight_resnums -- structure resnums belonging to the domain being
    visualized (if any). When given, that residue range is colored
    distinctly and the initial camera zooms to it instead of the whole
    structure -- without this, every domain's HTML shows the same
    whole-structure view and is indistinguishable from the others except
    for which spheres happen to be colored.

    show_variant_labels -- also float each variant's name (e.g. 'R2215W')
    next to its sphere. Off by default -- callers render a second,
    "_labeled" copy with this on, so both a clean and a labeled view exist
    side by side.
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
        if show_variant_labels:
            _add_variant_label(view, item)
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
    domains: list[Domain], alignment: AlignmentResult, resolved_resnums: set[int]
) -> list[tuple[Domain, list[int]]]:
    """(domain, sorted structure resnums) for every domain with >=1 resolved
    residue, in the given domain order -- shared by every caller that needs
    a domain's resnum set (per-domain zoom renders, and pipeline.py's
    Domain-mode region construction for both the primary chain-group
    (`resolved_resnums=set(struct.ca_coords)`) and any secondary chain-group
    aligned to its own UniProt accession
    (`resolved_resnums=set(struct.all_chain_ca_coords[chain_id])`))."""
    out = []
    for domain in domains:
        resnums = sorted(resnums_for_domain(domain, alignment) & resolved_resnums)
        if resnums:
            out.append((domain, resnums))
    return out


def _group_chains_by_sequence(
    struct: StructureData, chain_labels: dict[str, str] | None = None
) -> dict[str, list[str]]:
    """Chain ids grouped by identical sequence, e.g. the 3 PKD2 copies in a
    PKD1-PKD2 complex collapse to one group/color/legend entry instead of
    three -- "colored by chain" should mean by distinct molecule, not by
    every individual chain letter.

    `chain_labels` (chain id -> display name, e.g. {"D": "PKD1", "A": "PKD2"})
    lets the legend show real protein names instead of raw chain letters --
    structure files carry no such names themselves (an AlphaFold Server
    mmCIF's `_entity.pdbx_description` is empty), so this always has to come
    from the caller. Naming any one member of a group is enough to label the
    whole group.
    """
    chain_labels = chain_labels or {}
    groups: dict[str, list[str]] = {}
    for chain_id in struct.all_chain_ca_coords:
        seq = struct.all_chain_sequences.get(chain_id, "")
        groups.setdefault(seq, []).append(chain_id)
    labeled: dict[str, list[str]] = {}
    for chain_ids in groups.values():
        label = next((chain_labels[c] for c in chain_ids if c in chain_labels), ", ".join(chain_ids))
        labeled[label] = chain_ids
    return labeled


@dataclass
class ModeScheme:
    """One backbone-coloring option in the multi-mode overview's toggle --
    e.g. "Chain", "Domain". `regions` is applied on top of a shared plain
    gray/0.55-opacity base style (so anything not covered by any region
    stays that same gray, at the same opacity -- there's no separate
    "uncategorized" dimming). `legend_items` must already be deduped by name
    (a category like "Cytoplasmic" or a domain that spans several
    disjoint/discontiguous ranges legitimately produces multiple `regions`
    entries sharing one name)."""

    label: str
    regions: list[tuple[str, list[int], str]]  # (chain_id, resi_list, color_hex)
    legend_items: list[tuple[str, str]]  # (name, color), deduped, in display order


def render_multi_mode_overview_html(
    struct: StructureData,
    variants_df: pd.DataFrame,
    alignment: AlignmentResult,
    class_colors: ColorMap,
    modes: dict[str, ModeScheme],
    out_path: str | Path,
    *,
    title: str,
    cache_dir: str | Path,
    default_mode: str = "Chain",
    show_variant_labels: bool = False,
) -> Path:
    """Whole-structure interactive render with a <select> dropdown that
    switches the backbone coloring between several precomputed schemes
    (Chain / EM-AF / Topology / Domain -- see pipeline.run_render) without
    reloading the page. Replaces the old single-purpose
    render_domain_overview_html/render_chain_overview_html: those two modes
    are now just two entries in `modes`, and any protein/complex can add
    more (a new mode is just one more ModeScheme, no new render function).

    Variant spheres/labels are added once, independent of the mode switch --
    the toggle only ever repaints the cartoon backbone, never variant
    coloring (unchanged from every other render in this pipeline: colored by
    class via `class_colors`).

    Unlike every other renderer here, this one does NOT go through
    py3Dmol.view()/write_html() -- switching styles interactively needs a
    real, named JS handle on the live 3Dmol.js viewer object, which
    py3Dmol's generated glue doesn't expose. Instead this hand-writes a
    <script> block that calls the same underlying $3Dmol.createViewer/
    addModel/render API directly (3Dmol.min.js is inlined exactly as
    everywhere else in this module, so this stays just as self-contained/
    offline as every other render).
    """
    js_path = Path(cache_dir) / "js" / "3Dmol.min.js"
    if not js_path.exists():
        raise RenderError(
            f"no cached 3Dmol.min.js at {js_path} -- run "
            f"`protein-vis fetch --bootstrap-js` on the login node first"
        )
    js_text = js_path.read_text()

    if default_mode not in modes:
        raise RenderError(f"default_mode {default_mode!r} not among modes {list(modes)}")

    mapped, n_unmapped = _variant_positions_with_coords(variants_df, struct, alignment)

    mode_regions_json = {
        mode: [{"chain": chain_id, "resi": resi, "color": color} for chain_id, resi, color in scheme.regions]
        for mode, scheme in modes.items()
    }
    mode_legends_json = {
        mode: _build_legend_html(scheme.legend_items, heading=scheme.label) for mode, scheme in modes.items()
    }

    # Populates class_colors via .get() for every variant actually shown --
    # must run before class_colors.legend_items() below, which only reports
    # classes seen so far.
    variant_spheres = [
        {
            "x": float(item["coord"][0]),
            "y": float(item["coord"][1]),
            "z": float(item["coord"][2]),
            "color": class_colors.get(item["class_name"]),
            "label": item["raw"] if show_variant_labels else None,
        }
        for item in mapped
    ]
    class_legend_html = _build_legend_html(class_colors.legend_items(), heading="Variant class")

    footnote = f"{len(mapped)} variant(s) shown"
    if n_unmapped:
        footnote += f", {n_unmapped} unmapped (outside aligned/resolved structure region)"

    mode_options_html = "\n".join(
        f'<option value="{mode}"{" selected" if mode == default_mode else ""}>{mode}</option>'
        for mode in modes
    )

    # Only emitted at all when requested -- so, like every other renderer in
    # this module, the literal "addLabel"/variant-name text is simply absent
    # from a plain (non-labeled) render rather than present-but-inert.
    add_label_js = (
        """viewer.addLabel(v.label, {
      position: {x: v.x, y: v.y, z: v.z},
      backgroundColor: 'black',
      backgroundOpacity: 0.6,
      fontColor: 'white',
      fontSize: 11,
      showBackground: true,
      inFront: true,
    });"""
        if show_variant_labels
        else ""
    )

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
<div style="margin:8px;">
  <label style="font-size:13px;color:#333;">Color by:
    <select id="modeSelect">
{mode_options_html}
    </select>
  </label>
</div>
<div id="viewerContainer" style="width:1400px;height:750px;position:relative;"></div>
<div id="legendContainer" style="display:flex;flex-wrap:wrap;gap:24px;"></div>
<script>{js_text}</script>
<script>
(function() {{
  var element = document.getElementById('viewerContainer');
  var viewer = $3Dmol.createViewer(element, {{backgroundColor: 'white'}});
  viewer.addModel({json.dumps(struct.raw_text)}, {json.dumps(struct.fmt)});

  var modeRegions = {json.dumps(mode_regions_json)};
  var modeLegends = {json.dumps(mode_legends_json)};
  var classLegendHtml = {json.dumps(class_legend_html)};
  var variantSpheres = {json.dumps(variant_spheres)};

  function applyMode(mode) {{
    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray', opacity: 0.55}}}});
    (modeRegions[mode] || []).forEach(function(r) {{
      viewer.setStyle({{chain: r.chain, resi: r.resi}}, {{cartoon: {{color: r.color, opacity: 0.55}}}});
    }});
    document.getElementById('legendContainer').innerHTML = (modeLegends[mode] || '') + classLegendHtml;
    viewer.render();
  }}

  variantSpheres.forEach(function(v) {{
    viewer.addSphere({{
      center: {{x: v.x, y: v.y, z: v.z}},
      radius: 1.9,
      color: v.color,
    }});
    {add_label_js}
  }});

  applyMode({json.dumps(default_mode)});
  viewer.zoomTo();
  viewer.render();

  document.getElementById('modeSelect').addEventListener('change', function(e) {{
    applyMode(e.target.value);
  }});
}})();
</script>
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
