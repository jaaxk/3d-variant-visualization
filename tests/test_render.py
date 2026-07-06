from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from protein_vis.colors import ColorMap
from protein_vis.render import RenderError, render_interactive_html, render_static_png
from protein_vis.structure import AlignmentResult, StructureData

FIXTURES = Path(__file__).parent / "fixtures"


def _make_tiny_struct() -> StructureData:
    resnums = list(range(101, 111))
    ca_coords = {r: np.array([i * 3.8, 0.0, 0.0]) for i, r in enumerate(resnums)}
    return StructureData(
        chain_id="A",
        resnums=resnums,
        sequence="MKTAYIAKQR",
        ca_coords=ca_coords,
        raw_text=(FIXTURES / "tiny.pdb").read_text(),
        fmt="pdb",
    )


def _make_alignment() -> AlignmentResult:
    # reference positions 5-14 -> structure resnums 101-110
    return AlignmentResult(
        pos_to_resnum={i + 5: 101 + i for i in range(10)}, identity=1.0, coverage=0.55
    )


def _make_variants_df() -> pd.DataFrame:
    return pd.DataFrame(
        {"class_name": ["classA", "classB"], "raw": ["M5V", "A8G"], "pos": [5, 8]}
    )


def test_render_static_png_writes_file(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    colors = ColorMap()
    out = render_static_png(
        struct, _make_variants_df(), alignment, colors, tmp_path / "out.png", title="test"
    )
    assert out.exists()
    assert out.stat().st_size > 0


def test_render_interactive_html_is_self_contained(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    out = render_interactive_html(
        struct,
        _make_variants_df(),
        alignment,
        colors,
        tmp_path / "out.html",
        title="test",
        cache_dir=cache_dir,
    )
    html = out.read_text()
    assert "THIS_IS_THE_STUB_3DMOL_JS" in html  # js was inlined
    assert "cdn.jsdelivr" not in html
    assert "3dmol.org" not in html.lower()
    assert "$3Dmolpromise = Promise.resolve()" in html


def test_render_static_png_zooms_and_highlights_domain(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    colors = ColorMap()

    whole = render_static_png(
        struct, _make_variants_df(), alignment, colors, tmp_path / "whole.png", title="whole"
    )
    # Highlight just resnums 101-105 (the first half of the fixture backbone).
    zoomed = render_static_png(
        struct, _make_variants_df(), alignment, colors, tmp_path / "zoomed.png", title="zoomed",
        highlight_resnums={101, 102, 103, 104, 105},
    )
    assert whole.exists() and zoomed.exists()
    # A cropped render should differ in size from the uncropped one -- if
    # this ever fails, the highlight/zoom path silently stopped doing anything.
    assert whole.stat().st_size != zoomed.stat().st_size


def test_render_interactive_html_highlights_and_zooms_domain(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    out = render_interactive_html(
        struct,
        _make_variants_df(),
        alignment,
        colors,
        tmp_path / "out.html",
        title="test",
        cache_dir=cache_dir,
        highlight_resnums={101, 102, 103, 104, 105},
    )
    html = out.read_text()
    assert "steelblue" in html
    assert '"resi":[101,102,103,104,105]' in html.replace(" ", "")


def test_render_interactive_html_missing_js_raises(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    colors = ColorMap()
    with pytest.raises(RenderError):
        render_interactive_html(
            struct,
            _make_variants_df(),
            alignment,
            colors,
            tmp_path / "out.html",
            title="test",
            cache_dir=tmp_path / "empty_cache",
        )
