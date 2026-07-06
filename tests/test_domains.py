from pathlib import Path

import pandas as pd

from protein_vis.domains import (
    Domain,
    assign_domains,
    auto_domains_from_uniprot,
    group_variants_by_domain,
    load_domain_config,
)

FIXTURES = Path(__file__).parent / "fixtures"


def test_load_domain_config():
    domains = load_domain_config(FIXTURES / "tiny_domains.yaml")
    names = {d.name for d in domains}
    assert names == {"core_domain", "unused_domain"}


def test_assign_domains_overlap():
    domains = [
        Domain(name="A", start=1, end=10),
        Domain(name="B", start=5, end=15),
    ]
    assert assign_domains(7, domains) == ["A", "B"]
    assert assign_domains(12, domains) == ["B"]
    assert assign_domains(20, domains) == []


def test_group_variants_by_domain_excludes_zero_hit_domains():
    domains = load_domain_config(FIXTURES / "tiny_domains.yaml")
    variants_df = pd.DataFrame(
        {"class_name": ["classA", "classA"], "raw": ["M5V", "K6L"], "pos": [5, 6]}
    )
    groups = group_variants_by_domain(variants_df, domains)
    assert "core_domain" in groups
    assert "unused_domain" not in groups
    assert len(groups["core_domain"]) == 2


def test_group_variants_by_domain_explodes_overlapping_hits():
    domains = [Domain(name="A", start=1, end=10), Domain(name="B", start=5, end=15)]
    variants_df = pd.DataFrame({"class_name": ["x"], "raw": ["M7V"], "pos": [7]})
    groups = group_variants_by_domain(variants_df, domains)
    assert set(groups) == {"A", "B"}
    assert len(groups["A"]) == 1
    assert len(groups["B"]) == 1


def test_auto_domains_from_uniprot_uses_wanted_feature_types():
    domains = auto_domains_from_uniprot(FIXTURES / "tiny_uniprot_features.json")
    names = {d.name for d in domains}
    assert "Core_domain" in names
    assert "BRCA2_1" in names
    # "Chain" feature type is not in the wanted set and must be excluded.
    assert not any("Ignored" in n for n in names)
    assert len(domains) == 2
