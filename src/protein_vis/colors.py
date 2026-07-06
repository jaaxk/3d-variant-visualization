"""Class -> color assignment, shared across every render in a run.

Colors for the classification labels already used elsewhere in this user's
projects (dms_contrastive/inference/visualize_brca2.py) are reused verbatim so
that visualizations stay visually consistent across repos.
"""

from __future__ import annotations

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

DEFAULT_COLOR = "#9E9E9E"

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


class ColorMap:
    """Stateful class -> color lookup.

    Kept stateful (rather than a bare dict) so that a class name not in
    CLASS_COLORS gets the *same* fallback color, and the *same* legend
    position, across every domain sub-plot rendered in a single pipeline run.
    """

    def __init__(self, overrides: dict[str, str] | None = None) -> None:
        self._known = dict(CLASS_COLORS)
        if overrides:
            self._known.update(overrides)
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
            if self._next_fallback < len(FALLBACK_CYCLE):
                color = FALLBACK_CYCLE[self._next_fallback]
                self._next_fallback += 1
            else:
                color = DEFAULT_COLOR

        self._assigned[class_name] = color
        self._order.append(class_name)
        return color

    def legend_items(self) -> list[tuple[str, str]]:
        """Classes seen so far via get(), in first-seen order."""
        return [(name, self._assigned[name]) for name in self._order]
