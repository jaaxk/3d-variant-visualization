from protein_vis.colors import ColorMap, CLASS_COLORS, FALLBACK_CYCLE


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
