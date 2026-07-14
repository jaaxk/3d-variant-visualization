"""Full offline smoke test of the `render` path -- no network, no cluster.

Exercises variants -> structure -> alignment -> domain grouping -> rendering
end-to-end against tiny synthetic fixtures.
"""

import json
from pathlib import Path

from protein_vis import pipeline
from protein_vis import structure as structure_mod

FIXTURES = Path(__file__).parent / "fixtures"


def _setup_cache(tmp_path):
    cache_dir = tmp_path / "cache"
    structure_spec = f"file:{FIXTURES / 'tiny.pdb'}"
    structure_mod.fetch_structure(structure_spec, cache_dir)

    uniprot_dir = cache_dir / "uniprot"
    uniprot_dir.mkdir(parents=True)
    (uniprot_dir / "TEST0001.fasta").write_text((FIXTURES / "tiny_uniprot.fasta").read_text())
    (uniprot_dir / "TEST0001.json").write_text(
        (FIXTURES / "tiny_uniprot_features.json").read_text()
    )
    structure_mod._write_manifest_entry(
        cache_dir,
        "uniprot:TEST0001",
        {
            "url": None,
            "local_path": "uniprot/TEST0001.json",
            "fetched_at": "test",
            "sha256": "test",
        },
    )

    js_dir = cache_dir / "js"
    js_dir.mkdir(parents=True)
    (js_dir / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    return cache_dir, structure_spec


def test_pipeline_smoke_end_to_end(tmp_path):
    cache_dir, structure_spec = _setup_cache(tmp_path)
    output_dir = tmp_path / "output"

    pipeline.run_render(
        variants_csv=FIXTURES / "tiny_variants.csv",
        structure_spec=structure_spec,
        uniprot_accession="TEST0001",
        cache_dir=cache_dir,
        domains_arg=str(FIXTURES / "tiny_domains.yaml"),
        output_dir=output_dir,
        strict_wt=True,
        min_identity=0.5,
    )

    assert (output_dir / "overview.html").exists()
    assert (output_dir / "overview.png").exists()
    # domain_overview.html/chain_overview.html no longer exist -- both are
    # now just "Domain"/"Chain" toggle options inside overview.html.
    assert not (output_dir / "domain_overview.html").exists()
    assert not (output_dir / "chain_overview.html").exists()
    assert (output_dir / "core_domain.html").exists()
    assert (output_dir / "core_domain.png").exists()
    # unused_domain has zero assigned variants -- must not be rendered.
    assert not (output_dir / "unused_domain.html").exists()
    # every interactive HTML also gets a "_labeled" twin (variant names
    # floated next to each point) -- doubles the HTML output count.
    assert (output_dir / "overview_labeled.html").exists()
    assert (output_dir / "core_domain_labeled.html").exists()
    assert "addLabel" not in (output_dir / "overview.html").read_text()
    assert "addLabel" in (output_dir / "overview_labeled.html").read_text()

    # overview.html's Domain mode colors every configured domain that has
    # >=1 resolved residue -- core_domain (5-14) overlaps the structure's
    # resolved region and appears in the legend. (unused_domain, 1-2, falls
    # entirely outside the aligned/resolved region for this fixture
    # structure, so it correctly has nothing to color -- same filtering
    # every domain-coloring path in this pipeline applies.)
    overview_html = (output_dir / "overview.html").read_text()
    assert "core_domain" in overview_html
    assert '<option value="Chain"' in overview_html
    assert '<option value="Domain"' in overview_html
    assert '<option value="Topology"' in overview_html
    assert "EM/AF" not in overview_html  # no --provenance was passed

    report = json.loads((output_dir / "run_report.json").read_text())
    assert report["total_variants"] == 5
    assert "core_domain" in report["domains_rendered"]
    assert "unused_domain" not in report["domains_rendered"]
    assert report["domains_rendered"]["core_domain"]["n_unmapped"] == 0
    # Per-domain zoom renders record which structure residues were
    # highlighted/zoomed to.
    assert report["domains_rendered"]["core_domain"]["n_structure_residues_highlighted"] == 10
    assert report["overview"]["n_variants"] == 5
    assert set(report["overview"]["modes"]) == {"Chain", "Domain", "Topology"}


def test_pipeline_smoke_auto_domains(tmp_path):
    cache_dir, structure_spec = _setup_cache(tmp_path)
    output_dir = tmp_path / "output_auto"

    pipeline.run_render(
        variants_csv=FIXTURES / "tiny_variants.csv",
        structure_spec=structure_spec,
        uniprot_accession="TEST0001",
        cache_dir=cache_dir,
        domains_arg="auto",
        output_dir=output_dir,
        strict_wt=True,
        min_identity=0.5,
    )

    overview_html = (output_dir / "overview.html").read_text()
    # "Core_domain" (Domain feature, 5-14) overlaps the structure's resolved
    # region and shows up in the Domain mode's legend. ("BRCA2_1", a Repeat
    # feature at 1-2, falls entirely outside the aligned/resolved region for
    # this tiny fixture structure, so -- correctly -- it has no resolved
    # residues to color and is excluded, same as every other domain-coloring
    # path in this pipeline.)
    assert "Core_domain" in overview_html

    report = json.loads((output_dir / "run_report.json").read_text())
    # tiny_uniprot_features.json's "Core_domain" feature (5-14) should catch
    # all 5 fixture variants (positions 5-9).
    assert "Core_domain" in report["domains_rendered"]


def test_pipeline_smoke_class_color_and_variant_overrides(tmp_path):
    cache_dir, structure_spec = _setup_cache(tmp_path)
    output_dir = tmp_path / "output_overrides"

    pipeline.run_render(
        variants_csv=FIXTURES / "tiny_variants.csv",
        structure_spec=structure_spec,
        uniprot_accession="TEST0001",
        cache_dir=cache_dir,
        domains_arg=str(FIXTURES / "tiny_domains.yaml"),
        output_dir=output_dir,
        strict_wt=True,
        min_identity=0.5,
        class_color_overrides={"classA": "#E53935", "special": "#1E88E5"},
        variant_class_overrides={"M5V": "special"},
    )

    report = json.loads((output_dir / "run_report.json").read_text())
    # M5V was reassigned out of classA into the synthetic "special" class --
    # this takes precedence over its original class for both counting and
    # coloring, since the reassignment happens before anything else runs.
    assert report["class_counts"]["special"] == 1
    assert report["class_counts"]["classA"] == 2  # K6L, T7S remain

    html = (output_dir / "overview.html").read_text()
    assert "#1E88E5" in html  # special's overridden color
    assert "#E53935" in html  # classA's overridden color
    assert "special" in html  # legend entry for the synthetic class


def _setup_multichain_cache(tmp_path):
    """A 2-protein cache (chain A = TEST0001, chain B = TEST0002) for
    exercising chain_uniprot/domains_for/provenance/interface_json --
    mirrors the real PKD1(D)/PKD2(A) grafted-structure run's shape at
    fixture scale."""
    cache_dir = tmp_path / "cache"
    structure_spec = f"file:{FIXTURES / 'tiny_multichain.pdb'}"
    structure_mod.fetch_structure(structure_spec, cache_dir)

    uniprot_dir = cache_dir / "uniprot"
    uniprot_dir.mkdir(parents=True)
    for accession, fasta_name, json_name in [
        ("TEST0001", "tiny_uniprot.fasta", "tiny_uniprot_features.json"),
        ("TEST0002", "tiny_uniprot2.fasta", "tiny_uniprot2_features.json"),
    ]:
        (uniprot_dir / f"{accession}.fasta").write_text((FIXTURES / fasta_name).read_text())
        (uniprot_dir / f"{accession}.json").write_text((FIXTURES / json_name).read_text())
        structure_mod._write_manifest_entry(
            cache_dir, f"uniprot:{accession}",
            {"url": None, "local_path": f"uniprot/{accession}.json", "fetched_at": "test", "sha256": "test"},
        )

    js_dir = cache_dir / "js"
    js_dir.mkdir(parents=True)
    (js_dir / "3Dmol.min.js").write_text((FIXTURES / "3Dmol.min.js").read_text())

    return cache_dir, structure_spec


def test_pipeline_smoke_multichain_domains_topology_provenance_interface(tmp_path):
    cache_dir, structure_spec = _setup_multichain_cache(tmp_path)
    output_dir = tmp_path / "output_multichain"

    provenance_path = tmp_path / "provenance.json"
    provenance_path.write_text(json.dumps({"A": [101, 102, 103], "B": [201]}))

    interface_json = tmp_path / "interface.json"
    interface_json.write_text(json.dumps({"TEST0001": [6], "TEST0002": [1]}))

    pipeline.run_render(
        variants_csv=FIXTURES / "tiny_variants.csv",
        structure_spec=structure_spec,
        uniprot_accession="TEST0001",
        cache_dir=cache_dir,
        domains_arg=str(FIXTURES / "tiny_domains.yaml"),
        output_dir=output_dir,
        strict_wt=True,
        min_identity=0.5,
        chain_uniprot={"B": "TEST0002"},
        domains_for={"TEST0002": "auto"},
        provenance_path=provenance_path,
        interface_json=interface_json,
    )

    overview_html = (output_dir / "overview.html").read_text()
    # EM/AF mode only appears because --provenance was passed this time.
    assert '<option value="EM/AF"' in overview_html
    assert "6A70" in overview_html
    assert "AlphaFold2" in overview_html
    # Topology mode picks up TEST0002's Transmembrane feature via its own
    # chain-B alignment (auto_topology_from_uniprot), not just the primary
    # chain's.
    assert "Transmembrane" in overview_html
    # The interface positions (TEST0001 pos 6, TEST0002 pos 1) merge into
    # both accessions' Domain mode as one shared "Interface" legend entry.
    assert "Interface" in overview_html

    report = json.loads((output_dir / "run_report.json").read_text())
    assert set(report["overview"]["modes"]) == {"Chain", "Domain", "Topology", "EM/AF"}
