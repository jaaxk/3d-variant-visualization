"""Structural domain configuration and variant-to-domain assignment.

Domains are expressed as residue ranges in the reference protein's canonical
(e.g. UniProt) numbering -- the same numbering used by the variant CSV -- and
are mapped onto whichever structure is actually loaded via
structure.align_to_reference's position map, never via raw residue-number
equality against the structure file itself.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

from .structure import AlignmentResult


@dataclass(frozen=True)
class Domain:
    """A named residue set in reference-sequence numbering.

    Either a contiguous `start..end` range (the common case) or an explicit
    `positions` set (for genuinely discontiguous categories, e.g. a
    structural contact interface -- see scripts/compute_pkd1_pkd2_interface.py)
    -- exactly one of the two must be given.
    """

    name: str
    start: int | None = None
    end: int | None = None
    source: str | None = None
    positions: frozenset[int] | None = None

    def __post_init__(self) -> None:
        if self.positions is None and (self.start is None or self.end is None):
            raise ValueError(
                f"Domain {self.name!r} needs either start/end or positions"
            )

    def contains(self, pos: int) -> bool:
        if self.positions is not None:
            return pos in self.positions
        return self.start <= pos <= self.end


def load_domain_config(yaml_path: str | Path) -> list[Domain]:
    data = yaml.safe_load(Path(yaml_path).read_text())
    domains = []
    for entry in data.get("domains", []):
        positions = entry.get("positions")
        domains.append(
            Domain(
                name=entry["name"],
                start=int(entry["start"]) if "start" in entry else None,
                end=int(entry["end"]) if "end" in entry else None,
                source=entry.get("source"),
                positions=frozenset(int(p) for p in positions) if positions is not None else None,
            )
        )
    return domains


def auto_domains_from_uniprot(uniprot_json_path: str | Path) -> list[Domain]:
    """Fallback domain derivation for proteins without a curated YAML.

    Uses UniProt "Domain", "Region", and "Repeat" feature types. Note: for
    proteins like BRCA2 whose UniProt record only has functional interaction
    Regions and Repeats (no true structural "Domain" features), this yields
    a coarser split than a hand-curated config -- callers should prefer a
    curated YAML when one exists.
    """
    data = json.loads(Path(uniprot_json_path).read_text())
    wanted_types = {"Domain", "Region", "Repeat"}
    domains: list[Domain] = []
    name_counts: dict[str, int] = {}

    for feature in data.get("features", []):
        if feature.get("type") not in wanted_types:
            continue
        location = feature.get("location", {})
        start = location.get("start", {}).get("value")
        end = location.get("end", {}).get("value")
        if start is None or end is None:
            continue
        raw_name = feature.get("description") or feature.get("type")
        slug = "".join(c if c.isalnum() else "_" for c in raw_name).strip("_") or "region"
        name_counts[slug] = name_counts.get(slug, 0) + 1
        name = slug if name_counts[slug] == 1 else f"{slug}_{name_counts[slug]}"
        domains.append(Domain(name=name, start=int(start), end=int(end), source="UniProt (auto)"))

    return domains


def auto_topology_from_uniprot(uniprot_json_path: str | Path) -> list[Domain]:
    """Membrane-topology categories derived from UniProt's own
    "Transmembrane"/"Topological domain"/"Intramembrane" feature types --
    generalizes to any membrane protein with these UniProt annotations, no
    curated config needed. "Transmembrane" and "Intramembrane" features
    always bucket to "Transmembrane" regardless of their own description
    (helix name, "Pore-forming", etc.); "Topological domain" features bucket
    by their own description text ("Cytoplasmic"/"Extracellular"/"Lumenal"
    etc., taken verbatim from UniProt rather than guessed). Multiple ranges
    commonly share the same bucket name (e.g. every cytoplasmic loop) --
    that's fine for coloring (same category, same color) but callers that
    build a legend must dedupe by name first.
    """
    data = json.loads(Path(uniprot_json_path).read_text())
    domains: list[Domain] = []

    for feature in data.get("features", []):
        ftype = feature.get("type")
        location = feature.get("location", {})
        start = location.get("start", {}).get("value")
        end = location.get("end", {}).get("value")
        if start is None or end is None:
            continue

        if ftype in ("Transmembrane", "Intramembrane"):
            bucket = "Transmembrane"
        elif ftype == "Topological domain":
            bucket = (feature.get("description") or "").split(";")[0].strip() or "Extracellular"
        else:
            continue

        domains.append(
            Domain(name=bucket, start=int(start), end=int(end), source="UniProt (auto topology)")
        )

    return domains


def _domain_to_yaml_entry(d: Domain) -> dict:
    if d.positions is not None:
        return {"name": d.name, "positions": sorted(d.positions), "source": d.source}
    return {"name": d.name, "start": d.start, "end": d.end, "source": d.source}


def write_domain_config(
    domains: list[Domain], out_path: str | Path, *, accession: str, note: str
) -> None:
    payload = {
        "uniprot_accession": accession,
        "note": note,
        "domains": [_domain_to_yaml_entry(d) for d in domains],
    }
    Path(out_path).write_text(yaml.safe_dump(payload, sort_keys=False))


def resnums_for_domain(domain: Domain, alignment: AlignmentResult) -> set[int]:
    """Structure resnums covered by a domain's reference-sequence range.

    Used to highlight/zoom to a domain's own backbone segment (per-domain
    renders) or to color a domain's segment distinctly (whole-structure
    domain overview) -- without this, every domain's visualization would
    show the same whole-structure view apart from which variants happen to
    be colored.
    """
    positions = domain.positions if domain.positions is not None else range(domain.start, domain.end + 1)
    return {
        alignment.pos_to_resnum[pos]
        for pos in positions
        if pos in alignment.pos_to_resnum
    }


def assign_domains(pos: int, domains: list[Domain]) -> list[str]:
    """All domains containing pos (inclusive range, or explicit position set
    for discontiguous domains). Overlaps are allowed by design -- e.g. a
    sub-domain nested within a larger domain's span -- so a variant may
    legitimately belong to more than one domain."""
    return [d.name for d in domains if d.contains(pos)]


def group_variants_by_domain(
    variants_df: pd.DataFrame, domains: list[Domain]
) -> dict[str, pd.DataFrame]:
    """Explode variants into per-domain groups (a multi-domain variant appears
    in every matching group). Only domains with >=1 assigned variant are
    returned; callers should always additionally render the full, unfiltered
    overview separately."""
    rows_by_domain: dict[str, list[int]] = {}
    for idx, pos in variants_df["pos"].items():
        for name in assign_domains(pos, domains):
            rows_by_domain.setdefault(name, []).append(idx)

    return {
        name: variants_df.loc[indices].reset_index(drop=True)
        for name, indices in rows_by_domain.items()
    }
