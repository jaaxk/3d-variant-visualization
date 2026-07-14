"""Class -> color assignment, shared across every render in a run.

Colors for the classification labels already used elsewhere in this user's
projects (dms_contrastive/inference/visualize_brca2.py) are reused verbatim so
that visualizations stay visually consistent across repos.
"""

from __future__ import annotations

import colorsys

# Reused verbatim from dms_contrastive/inference/visualize_brca2.py for
# cross-project visual consistency.
CLASS_COLORS: dict[str, str] = {
    "Benign": "#2196F3",
    "Pathogenic": "#F44336",
    "Uncertain": "#BDBDBD",
    "Functional": "#2196F3",
    "Non-functional": "#F44336",
    "Intermediate": "#BDBDBD",
    "hypomorphic": "#FB8C00",
}

# Fixed (not generated) -- topology is always exactly these 3 categories,
# and provenance always exactly these 3, regardless of protein.
TOPOLOGY_COLORS: dict[str, str] = {
    "Cytoplasmic": "#66BB6A",
    "Transmembrane": "#FFA726",
    "Extracellular": "#42A5F5",
}

PROVENANCE_COLORS: dict[str, str] = {
    "6A70": "#8E24AA",
    "AlphaFold2: Complex": "#00897B",
    "AlphaFold2: PKD1 monomer": "#EF6C00",
}

# AlphaFold's own official per-residue pLDDT confidence palette (same 4
# bands/colors AlphaFold DB and the AlphaFold Server use in their own
# confidence-colored structure views) -- reusing this rather than inventing
# a new scheme keeps the coloring immediately recognizable to anyone who's
# used AlphaFold's own tools. "Experimentally resolved" is not a confidence
# band at all (a real deposited structure -- e.g. 6A70 -- has no pLDDT), so
# it gets its own fixed color rather than being folded into the bands;
# reuses PROVENANCE_COLORS' "6A70" purple for visual consistency with the
# EM/AF mode's own "this part is real data" color.
CONFIDENCE_COLORS: dict[str, str] = {
    "Experimentally resolved": PROVENANCE_COLORS["6A70"],
    "Very high (pLDDT > 90)": "#0053D6",
    "Confident (70-90)": "#65CBF3",
    "Low (50-70)": "#FFDB13",
    "Very low (< 50)": "#FF7D45",
}

DEFAULT_COLOR = "#9E9E9E"


def confidence_bucket(plddt: float) -> str:
    """Bucket a per-residue pLDDT value into AlphaFold's own 4 confidence
    bands (see CONFIDENCE_COLORS) -- same thresholds AlphaFold DB/Server use."""
    if plddt > 90:
        return "Very high (pLDDT > 90)"
    if plddt > 70:
        return "Confident (70-90)"
    if plddt > 50:
        return "Low (50-70)"
    return "Very low (< 50)"

# Assigned in order to any class name not present in CLASS_COLORS.
FALLBACK_CYCLE: list[str] = [
    "#4CAF50",
    "#9C27B0",
    "#00BCD4",
    "#795548",
    "#3F51B5",
    "#8BC34A",
    "#E91E63",
]


def generate_categorical_palette(n: int) -> list[str]:
    """Generate `n` visually-distinct hex colors for cases (like a protein's
    full domain architecture) where the count of categories genuinely
    exceeds what any fixed qualitative palette can distinguish.

    NOTE: this deliberately departs from the `dataviz` skill's categorical
    color formula, whose "8 fixed hue anchors, never generated, fold a 9th
    series into Other/small-multiples/composite-encoding" rule assumes a
    business-chart series count. A protein's domain architecture is a
    different kind of artifact -- an established structural-biology
    convention (PyMOL/ChimeraX "spectrum" domain coloring, UniProt's own
    feature viewer) where every domain needs its own identifiable color
    regardless of count, precisely because the whole point of this view is
    seeing all of them at once (folding extras into "Other" or splitting
    into small multiples would defeat that -- the per-domain zoomed renders
    already are the small-multiples view). Beyond ~12 categories no palette
    stays reliably colorblind-distinguishable (the dataviz skill's own CVD
    check tops out there); this trades that off deliberately, and the
    legend always prints each domain's name as text directly next to its
    swatch so identity is never color-alone.

    Colors are evenly spaced around the hue wheel at fixed
    saturation/lightness chosen to keep every hue legible against both a
    light and a dark HTML background (mid-lightness, high-ish saturation --
    approximating the dataviz skill's OKLCH lightness-band/chroma-floor
    intent without pulling in an OKLCH-capable dependency for this one
    secondary, non-dashboard visualization).
    """
    if n <= 0:
        return []
    saturation = 0.60
    lightness = 0.55
    colors = []
    for i in range(n):
        hue = i / n
        r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        colors.append("#{:02X}{:02X}{:02X}".format(round(r * 255), round(g * 255), round(b * 255)))
    return colors


class ColorMap:
    """Stateful class -> color lookup.

    Kept stateful (rather than a bare dict) so that a class name not in
    CLASS_COLORS gets the *same* fallback color, and the *same* legend
    position, across every domain sub-plot rendered in a single pipeline run.
    """

    def __init__(
        self,
        overrides: dict[str, str] | None = None,
        fallback_cycle: list[str] | None = None,
    ) -> None:
        self._known = dict(CLASS_COLORS)
        if overrides:
            self._known.update(overrides)
        self._fallback_cycle = fallback_cycle if fallback_cycle is not None else FALLBACK_CYCLE
        self._assigned: dict[str, str] = {}
        self._order: list[str] = []
        self._next_fallback = 0

    def get(self, class_name: str) -> str:
        if class_name in self._assigned:
            return self._assigned[class_name]

        color = self._known.get(class_name)
        if color is None:
            # Case-insensitive fallback lookup before assigning a new color.
            lowered = class_name.lower()
            for known_name, known_color in self._known.items():
                if known_name.lower() == lowered:
                    color = known_color
                    break

        if color is None:
            if self._next_fallback < len(self._fallback_cycle):
                color = self._fallback_cycle[self._next_fallback]
                self._next_fallback += 1
            else:
                color = DEFAULT_COLOR

        self._assigned[class_name] = color
        self._order.append(class_name)
        return color

    def legend_items(self) -> list[tuple[str, str]]:
        """Classes seen so far via get(), in first-seen order."""
        return [(name, self._assigned[name]) for name in self._order]
