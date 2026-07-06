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
    assert (output_dir / "core_domain.html").exists()
    assert (output_dir / "core_domain.png").exists()
    # unused_domain has zero assigned variants -- must not be rendered.
    assert not (output_dir / "unused_domain.html").exists()

    report = json.loads((output_dir / "run_report.json").read_text())
    assert report["total_variants"] == 5
    assert "core_domain" in report["domains_rendered"]
    assert "unused_domain" not in report["domains_rendered"]
    assert report["domains_rendered"]["core_domain"]["n_unmapped"] == 0
    # Domain renders should record which structure residues were
    # highlighted/zoomed to; the overview has no single domain to highlight.
    assert report["domains_rendered"]["core_domain"]["n_structure_residues_highlighted"] == 10
    assert report["domains_rendered"]["overview"]["n_structure_residues_highlighted"] is None


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

    report = json.loads((output_dir / "run_report.json").read_text())
    # tiny_uniprot_features.json's "Core_domain" feature (5-14) should catch
    # all 5 fixture variants (positions 5-9).
    assert "Core_domain" in report["domains_rendered"]
