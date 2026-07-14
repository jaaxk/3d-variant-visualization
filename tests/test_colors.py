from protein_vis.colors import (
    CLASS_COLORS,
    CONFIDENCE_COLORS,
    FALLBACK_CYCLE,
    ColorMap,
    confidence_bucket,
    generate_categorical_palette,
)


def test_known_class_uses_exact_palette():
    cm = ColorMap()
    assert cm.get("Pathogenic") == CLASS_COLORS["Pathogenic"]
    assert cm.get("hypomorphic") == CLASS_COLORS["hypomorphic"]


def test_unknown_class_gets_stable_fallback():
    cm = ColorMap()
    first = cm.get("mystery_class")
    second = cm.get("mystery_class")
    assert first == second
    assert first == FALLBACK_CYCLE[0]


def test_fallback_cycle_advances_per_distinct_class():
    cm = ColorMap()
    colors = [cm.get(f"class_{i}") for i in range(3)]
    assert colors == FALLBACK_CYCLE[:3]


def test_legend_items_preserve_first_seen_order():
    cm = ColorMap()
    cm.get("Benign")
    cm.get("Pathogenic")
    cm.get("Benign")  # repeat should not duplicate or reorder
    assert cm.legend_items() == [
        ("Benign", CLASS_COLORS["Benign"]),
        ("Pathogenic", CLASS_COLORS["Pathogenic"]),
    ]


def test_generate_categorical_palette_returns_n_unique_colors():
    for n in (0, 1, 3, len(FALLBACK_CYCLE), len(FALLBACK_CYCLE) + 1, 31):
        palette = generate_categorical_palette(n)
        assert len(palette) == n
        assert len(set(palette)) == n
        for color in palette:
            assert color.startswith("#") and len(color) == 7


def test_colormap_uses_custom_fallback_cycle():
    palette = generate_categorical_palette(5)
    cm = ColorMap(fallback_cycle=palette)
    assigned = [cm.get(f"domain_{i}") for i in range(5)]
    assert assigned == palette


def test_confidence_bucket_thresholds():
    assert confidence_bucket(95) == "Very high (pLDDT > 90)"
    assert confidence_bucket(90) == "Confident (70-90)"  # boundary: not strictly > 90
    assert confidence_bucket(80) == "Confident (70-90)"
    assert confidence_bucket(70) == "Low (50-70)"
    assert confidence_bucket(60) == "Low (50-70)"
    assert confidence_bucket(50) == "Very low (< 50)"
    assert confidence_bucket(10) == "Very low (< 50)"


def test_confidence_colors_has_experimentally_resolved_and_four_bands():
    assert set(CONFIDENCE_COLORS) == {
        "Experimentally resolved",
        "Very high (pLDDT > 90)",
        "Confident (70-90)",
        "Low (50-70)",
        "Very low (< 50)",
    }
