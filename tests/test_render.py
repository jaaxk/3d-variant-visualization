from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from protein_vis.colors import ColorMap, generate_categorical_palette
from protein_vis.domains import Domain
from protein_vis.render import (
    ModeScheme,
    RenderError,
    _build_checkbox_legend_html,
    render_interactive_html,
    render_multi_mode_overview_html,
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


def _make_domain_mode_scheme() -> ModeScheme:
    domains = _make_domains()
    domain_colors = ColorMap(fallback_cycle=generate_categorical_palette(len(domains)))
    alignment = _make_alignment()
    struct = _make_tiny_struct()
    regions = []
    legend = []
    for domain in domains:
        resnums = sorted(
            {alignment.pos_to_resnum[p] for p in range(domain.start, domain.end + 1)
             if p in alignment.pos_to_resnum} & set(struct.ca_coords)
        )
        color = domain_colors.get(domain.name)
        regions.append((struct.chain_id, resnums, color, domain.name))
        legend.append((domain.name, color))
    return ModeScheme(label="Domain", regions=regions, legend_items=legend)


def _make_chain_mode_scheme(struct: StructureData, chain_labels=None) -> ModeScheme:
    from protein_vis.render import _group_chains_by_sequence

    chain_colors = ColorMap()
    groups = _group_chains_by_sequence(struct, chain_labels)
    regions, legend = [], []
    for label, chain_ids in groups.items():
        color = chain_colors.get(label)
        legend.append((label, color))
        for chain_id in chain_ids:
            regions.append((chain_id, sorted(struct.all_chain_ca_coords[chain_id]), color, label))
    return ModeScheme(label="Chain", regions=regions, legend_items=legend)


def test_render_multi_mode_overview_writes_all_modes(tmp_path):
    struct = _make_multichain_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    modes = {
        "Chain": _make_chain_mode_scheme(struct, chain_labels={"D": "PKD1", "A": "PKD2"}),
        "Domain": _make_domain_mode_scheme(),
    }

    out = render_multi_mode_overview_html(
        struct, _make_variants_df(), alignment, class_colors, modes,
        tmp_path / "overview.html", title="test", cache_dir=cache_dir,
    )
    html = out.read_text()
    # Every mode's regions/legend are embedded up front (so the dropdown can
    # switch client-side without reloading) -- not just the default mode's.
    assert "domain_one" in html
    assert "domain_two" in html
    assert "PKD1" in html
    assert "PKD2" in html
    assert '<option value="Chain" selected>' in html
    assert '<option value="Domain">' in html
    assert "modeSelect" in html
    assert "classA" in html and "classB" in html
    # The legend must be visible in the raw HTML itself (statically
    # pre-populated for the default mode), not only injected by JS at
    # runtime -- so it's never blank regardless of how/whether the 3Dmol.js
    # viewer script actually executes in the viewer.
    legend_div = html[html.index('id="legendContainer"'):html.index("</div>", html.index('id="legendContainer"'))]
    assert "PKD1" in legend_div  # default mode is Chain
    assert "classA" in legend_div
    # self-contained/offline invariant must still hold for this render path too.
    assert "THIS_IS_THE_STUB_3DMOL_JS" in html
    assert "cdn.jsdelivr" not in html
    assert "3dmol.org" not in html.lower()


def test_render_multi_mode_overview_legend_includes_unmapped_classes(tmp_path):
    """A variant class with zero markers actually shown on THIS structure
    (e.g. it only covers a region the structure doesn't resolve) must still
    appear in the legend, flagged "(0 shown)" -- never silently vanish, since
    that's indistinguishable from the class not existing at all."""
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    variants_df = pd.DataFrame(
        {
            "class_name": ["classA", "classB", "unmapped_class"],
            "raw": ["M5V", "A8G", "M5X"],
            "pos": [5, 8, 999],  # 999 is outside the aligned/resolved region
        }
    )
    modes = {"Domain": _make_domain_mode_scheme()}

    out = render_multi_mode_overview_html(
        struct, variants_df, alignment, class_colors, modes,
        tmp_path / "overview.html", title="test", cache_dir=cache_dir, default_mode="Domain",
    )
    html = out.read_text()
    assert "classA" in html
    assert "unmapped_class (0 shown)" in html


def test_render_multi_mode_overview_respects_default_mode(tmp_path):
    struct = _make_multichain_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    modes = {
        "Chain": _make_chain_mode_scheme(struct),
        "Domain": _make_domain_mode_scheme(),
    }
    out = render_multi_mode_overview_html(
        struct, _make_variants_df(), alignment, class_colors, modes,
        tmp_path / "overview.html", title="test", cache_dir=cache_dir, default_mode="Domain",
    )
    html = out.read_text()
    assert '<option value="Domain" selected>' in html
    assert '<option value="Chain">' in html


def test_render_multi_mode_overview_rejects_unknown_default_mode(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    modes = {"Chain": _make_chain_mode_scheme(_make_multichain_struct())}
    with pytest.raises(RenderError):
        render_multi_mode_overview_html(
            struct, _make_variants_df(), alignment, class_colors, modes,
            tmp_path / "overview.html", title="test", cache_dir=cache_dir, default_mode="Nope",
        )


def test_render_interactive_html_shows_variant_labels_when_requested(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    plain = render_interactive_html(
        struct, _make_variants_df(), alignment, colors, tmp_path / "plain.html",
        title="test", cache_dir=cache_dir,
    )
    labeled = render_interactive_html(
        struct, _make_variants_df(), alignment, colors, tmp_path / "labeled.html",
        title="test", cache_dir=cache_dir, show_variant_labels=True,
    )
    # addLabel calls (and hence the variant names) should only appear in the
    # labeled render -- the plain one stays exactly as before.
    assert "addLabel" not in plain.read_text()
    labeled_html = labeled.read_text()
    assert "addLabel" in labeled_html
    assert "M5V" in labeled_html
    assert "A8G" in labeled_html


def test_render_multi_mode_overview_shows_variant_labels_when_requested(tmp_path):
    struct = _make_tiny_struct()
    alignment = _make_alignment()
    class_colors = ColorMap()
    cache_dir = tmp_path / "cache"
    (cache_dir / "js").mkdir(parents=True)
    (cache_dir / "js" / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())
    modes = {"Domain": _make_domain_mode_scheme()}

    plain = render_multi_mode_overview_html(
        struct, _make_variants_df(), alignment, class_colors, modes,
        tmp_path / "plain.html", title="test", cache_dir=cache_dir, default_mode="Domain",
    )
    labeled = render_multi_mode_overview_html(
        struct, _make_variants_df(), alignment, class_colors, modes,
        tmp_path / "labeled.html", title="test", cache_dir=cache_dir, default_mode="Domain",
        show_variant_labels=True,
    )
    assert "addLabel" not in plain.read_text()
    labeled_html = labeled.read_text()
    assert "addLabel" in labeled_html
    assert "M5V" in labeled_html
    assert "A8G" in labeled_html


def test_build_checkbox_legend_html_has_checkboxes_and_uncheck_all():
    html = _build_checkbox_legend_html([("domain_one", "#111111"), ("domain_two", "#222222")], heading="Domain")
    assert html.count('class="modeCategoryCheckbox"') == 2
    assert 'data-name="domain_one"' in html
    assert 'data-name="domain_two"' in html
    assert "checked" in html
    assert 'class="modeUncheckAllBtn"' in html
    assert "Uncheck all" in html


def test_build_checkbox_legend_html_no_button_when_empty():
    html = _build_checkbox_legend_html([], heading="Domain")
    assert "modeUncheckAllBtn" not in html
