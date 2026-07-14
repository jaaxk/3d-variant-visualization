from pathlib import Path

import numpy as np
import pytest

from protein_vis.structure import (
    AlignmentError,
    StructureData,
    align_to_reference,
    fetch_structure,
    load_structure,
    parse_source_spec,
)

FIXTURES = Path(__file__).parent / "fixtures"
REFERENCE_SEQ = "EEEEMKTAYIAKQREEEE"  # positions 5-14 = MKTAYIAKQR (matches tiny.pdb)


@pytest.mark.parametrize(
    "spec,kind,acc,chain",
    [
        ("alphafold:P51587", "alphafold", "P51587", None),
        ("pdb:1IYJ", "pdb", "1IYJ", None),
        ("pdb:1IYJ:B", "pdb", "1IYJ", "B"),
        ("file:/tmp/foo.pdb", "file", "/tmp/foo.pdb", None),
    ],
)
def test_parse_source_spec(spec, kind, acc, chain):
    source = parse_source_spec(spec)
    assert source.kind == kind
    assert source.accession_or_path == acc
    assert source.chain == chain


def test_parse_source_spec_rejects_unknown_kind():
    with pytest.raises(ValueError):
        parse_source_spec("bogus:1234")


def _make_tiny_struct() -> StructureData:
    core = "MKTAYIAKQR"
    resnums = list(range(101, 111))
    ca_coords = {r: np.array([i * 3.8, 0.0, 0.0]) for i, r in enumerate(resnums)}
    return StructureData(
        chain_id="A",
        resnums=resnums,
        sequence=core,
        ca_coords=ca_coords,
        raw_text="",
        fmt="pdb",
    )


def test_align_to_reference_maps_positions_correctly():
    struct = _make_tiny_struct()
    result = align_to_reference(struct, REFERENCE_SEQ, min_identity=0.5)
    assert result.identity > 0.99
    # reference position 5 ("M") should map to structure resnum 101
    assert result.pos_to_resnum[5] == 101
    # reference position 14 ("R") should map to structure resnum 110
    assert result.pos_to_resnum[14] == 110
    assert 0.0 < result.coverage <= 1.0


def test_align_to_reference_raises_below_min_identity():
    struct = _make_tiny_struct()
    unrelated_reference = "QWERTYQWERTYQWERTYQWERTY"
    with pytest.raises(AlignmentError):
        align_to_reference(struct, unrelated_reference, min_identity=0.9)


def test_load_structure_extracts_ca_bfactor(tmp_path):
    """CA atom B-factor is a real crystallographic/EM B-factor for
    experimental structures, but AlphaFold (DB and Server) writes per-atom
    pLDDT into that same column -- pipeline.py's --confidence mode reads it
    from here, so load_structure must always carry it through regardless of
    which interpretation actually applies for a given file."""
    pdb_text = (
        "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 95.00           C\n"
        "ATOM      2  CA  LYS A   2       3.800   0.000   0.000  1.00 42.50           C\n"
        "TER\nEND\n"
    )
    pdb_path = tmp_path / "bfactor.pdb"
    pdb_path.write_text(pdb_text)
    cache_dir = tmp_path / "cache"
    spec = f"file:{pdb_path}"
    fetch_structure(spec, cache_dir)
    struct = load_structure(spec, cache_dir)
    assert struct.all_chain_ca_bfactor["A"][1] == pytest.approx(95.00)
    assert struct.all_chain_ca_bfactor["A"][2] == pytest.approx(42.50)
