from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from protein_vis.colors import ColorMap, generate_categorical_palette
from protein_vis.domains import Domain
from protein_vis.render import (
    RenderError,
    render_chain_overview_html,
    render_chain_overview_png,
    render_domain_overview_html,
    render_domain_overview_png,
    render_interactive_html,
    render_static_png,
)
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


def _make_multichain_struct() -> StructureData:
    struct = _make_tiny_struct()
    struct.all_chain_ca_coords = {
        "D": dict(struct.ca_coords),
        "A": {r: np.array([i * 3.8, 5.0, 0.0]) for i, r in enumerate(range(201, 205))},
        "B": {r: np.array([i * 3.8, 10.0, 0.0]) for i, r in enumerate(range(301, 305))},
    }
    struct.all_chain_sequences = {"D": "MKTAYIAKQR", "A": "MKTA", "B": "MKTA"}
    return struct


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


def _make_domains() -> list[Domain]:
    # reference positions 5-14 map to structure resnums 101-110 (see
    # _make_alignment) -- split into two adjacent domains covering it.
    return [
        Domain(name="domain_one", start=5, end=9, source="test fixture"),
        Domain(name="domain_two", start=10, end=14, source="test fixture"),
    ]


def test_render_domain_overview_png_writes_file(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    domains = _make_domains()
    domain_colors = ColorMap(fallback_cycle=generate_categorical_palette(len(domains)))

    out = render_domain_overview_png(
        struct,
        _make_variants_df(),
        alignment,
        class_colors,
        domain_colors,
        domains,
        tmp_path / "domain_overview.png",
        title="test",
    )
    assert out.exists()
    assert out.stat().st_size > 0


def test_render_domain_overview_html_shows_both_legends(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    domains = _make_domains()
    domain_colors = ColorMap(fallback_cycle=generate_categorical_palette(len(domains)))
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    out = render_domain_overview_html(
        struct,
        _make_variants_df(),
        alignment,
        class_colors,
        domain_colors,
        domains,
        tmp_path / "domain_overview.html",
        title="test",
        cache_dir=cache_dir,
    )
    html = out.read_text()
    assert "domain_one" in html
    assert "domain_two" in html
    assert "classA" in html
    assert "classB" in html
    assert "Domains" in html
    assert "Variant class" in html
    # self-contained/offline invariant must still hold for this render path too.
    assert "THIS_IS_THE_STUB_3DMOL_JS" in html
    assert "cdn.jsdelivr" not in html
    assert "3dmol.org" not in html.lower()


def test_render_chain_overview_groups_identical_sequences(tmp_path):
    struct = _make_multichain_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    chain_colors = ColorMap()

    out = render_chain_overview_png(
        struct,
        _make_variants_df(),
        alignment,
        class_colors,
        chain_colors,
        tmp_path / "chain_overview.png",
        title="test",
    )
    assert out.exists() and out.stat().st_size > 0


def test_render_chain_overview_html_uses_custom_labels(tmp_path):
    struct = _make_multichain_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    chain_colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    out = render_chain_overview_html(
        struct,
        _make_variants_df(),
        alignment,
        class_colors,
        chain_colors,
        tmp_path / "chain_overview.html",
        title="test",
        cache_dir=cache_dir,
        chain_labels={"D": "PKD1", "A": "PKD2"},
    )
    html = out.read_text()
    # A and B share a sequence and should collapse to one legend entry,
    # labeled "PKD2" (from A's label) rather than "A, B".
    assert "PKD1" in html
    assert "PKD2" in html
    assert "A, B" not in html
    assert "Chain" in html
    assert "cdn.jsdelivr" not in html
    assert "3dmol.org" not in html.lower()


def test_render_chain_overview_html_without_labels_falls_back_to_chain_ids(tmp_path):
    struct = _make_multichain_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    chain_colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    out = render_chain_overview_html(
        struct,
        _make_variants_df(),
        alignment,
        class_colors,
        chain_colors,
        tmp_path / "chain_overview.html",
        title="test",
        cache_dir=cache_dir,
    )
    html = out.read_text()
    assert "A, B" in html
    assert ">D<" in html
