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
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import py3Dmol  # noqa: E402

from .colors import ColorMap
from .structure import AlignmentResult, StructureData


class RenderError(RuntimeError):
    pass


def _build_legend_html(items: list[tuple[str, str]]) -> str:
    rows = "\n".join(
        f'<li><span style="display:inline-block;width:12px;height:12px;'
        f'background:{color};margin-right:6px;border-radius:50%;"></span>{name}</li>'
        for name, color in items
    )
    return f'<ul style="list-style:none;padding:0;margin:8px;font-family:sans-serif;">{rows}</ul>'


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
) -> Path:
    mapped, n_unmapped = _variant_positions_with_coords(variants_df, struct, alignment)

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    backbone = np.array(
        [struct.ca_coords[r] for r in sorted(struct.resnums) if r in struct.ca_coords]
    )
    if len(backbone) > 0:
        ax.plot(*backbone.T, "-", lw=0.6, color="#B0BEC5", alpha=0.6)

    by_class: dict[str, list[np.ndarray]] = {}
    for item in mapped:
        by_class.setdefault(item["class_name"], []).append(item["coord"])

    for class_name, coords in by_class.items():
        arr = np.array(coords)
        ax.scatter(*arr.T, color=colors.get(class_name), s=40, edgecolors="white", label=class_name)

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
) -> Path:
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
    for item in mapped:
        x, y, z = (float(c) for c in item["coord"])
        view.addSphere(
            {
                "center": {"x": x, "y": y, "z": z},
                "radius": 1.2,
                "color": colors.get(item["class_name"]),
            }
        )
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
