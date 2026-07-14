#!/usr/bin/env python3
"""Compute the real PKD1-PKD2 contact interface from PDB 6A70's deposited
coordinates, for pkd1 clustering / DMS-variant interface classification.

Why: the user wants an "interface"/"binding" domain for the PKD1-PKD2
overview coloring, but no UniProt feature is literally annotated as such.
Rather than guess, this computes it directly from the real cryo-EM
coordinates: any residue with a heavy atom within `--cutoff` (default 5 A)
of the other protein's chain in 6A70 is an interface residue. Run once
against the structure already cached by `protein-vis fetch` -- this is a
property of the real experimental complex, not of whichever predicted
structure gets visualized later, so it never needs to be re-run per-run.

Output: configs/domains/pkd1_pkd2_interface.json --
    {"P98161": [uniprot_pos, ...], "Q13563": [uniprot_pos, ...]}
a plain position-list artifact, independent of protein_vis's own domain
YAML format, meant to be joined against any variant table later (a DMS CSV
with a `pos` column gets an `is_interface` column via
`df["pos"].isin(interface["P98161"])`) as well as consumed by
protein_vis's Domain coloring mode.

Usage (cheap, CPU-only -- run via srun, not the login node):
    python scripts/compute_pkd1_pkd2_interface.py \\
        --cache-dir /scratch/jv2807/pkd1/structure_cache \\
        --output configs/domains/pkd1_pkd2_interface.json
"""

from __future__ import annotations

import argparse
import json
from io import StringIO
from pathlib import Path

from Bio.PDB import NeighborSearch, PDBParser

from protein_vis import structure as structure_mod

PKD1_ACCESSION = "P98161"
PKD2_ACCESSION = "Q13563"
PKD1_CHAIN = "B"
PKD2_CHAINS = ("A", "F", "G")  # 3 copies of PKD2 in 6A70; only the truly
                                # adjacent copy will actually contribute contacts.


def _heavy_atoms(chain) -> list:
    return [atom for res in chain for atom in res if res.id[0] == " "]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", required=True,
                     help="protein_vis structure cache (must already have 'pdb:6A70' and "
                          "the P98161/Q13563 UniProt references fetched).")
    ap.add_argument("--output", required=True, help="Output JSON path.")
    ap.add_argument("--cutoff", type=float, default=5.0, help="Contact distance cutoff, Angstroms.")
    args = ap.parse_args()

    cache_dir = Path(args.cache_dir)

    print(f"[1/3] loading 6A70 chain {PKD1_CHAIN} (PKD1) and aligning to {PKD1_ACCESSION}")
    pkd1_struct = structure_mod.load_structure(f"pdb:6A70:{PKD1_CHAIN}", cache_dir)
    pkd1_ref = structure_mod.load_uniprot_sequence(cache_dir, PKD1_ACCESSION)
    pkd1_alignment = structure_mod.align_to_reference(pkd1_struct, pkd1_ref)
    print(f"      identity={pkd1_alignment.identity:.1%} coverage={pkd1_alignment.coverage:.1%}")
    pkd1_resnum_to_pos = {resnum: pos for pos, resnum in pkd1_alignment.pos_to_resnum.items()}

    print(f"[2/3] loading 6A70 chains {PKD2_CHAINS} (PKD2 copies) and aligning to {PKD2_ACCESSION}")
    pkd2_ref = structure_mod.load_uniprot_sequence(cache_dir, PKD2_ACCESSION)
    pkd2_resnum_to_pos: dict[str, dict[int, int]] = {}
    for chain_id in PKD2_CHAINS:
        struct = structure_mod.load_structure(f"pdb:6A70:{chain_id}", cache_dir)
        alignment = structure_mod.align_to_reference(struct, pkd2_ref)
        print(f"      chain {chain_id}: identity={alignment.identity:.1%} coverage={alignment.coverage:.1%}")
        pkd2_resnum_to_pos[chain_id] = {
            resnum: pos for pos, resnum in alignment.pos_to_resnum.items()
        }

    print(f"[3/3] computing heavy-atom contacts (<= {args.cutoff} A) directly from 6A70's coordinates")
    # pkd1_struct.raw_text is the *whole* 6A70 file (load_structure reads the full file
    # text regardless of which single chain it auto-selects), so re-parsing it directly
    # gives every chain's full atom list without touching the cache manifest again.
    model = next(iter(PDBParser(QUIET=True).get_structure("6A70", StringIO(pkd1_struct.raw_text))))

    pkd1_atoms = _heavy_atoms(model[PKD1_CHAIN])
    ns = NeighborSearch(pkd1_atoms)

    pkd1_positions: set[int] = set()
    pkd2_positions: set[int] = set()
    for chain_id in PKD2_CHAINS:
        for atom in _heavy_atoms(model[chain_id]):
            neighbors = ns.search(atom.coord, args.cutoff)
            if not neighbors:
                continue
            pkd2_resnum = atom.get_parent().id[1]
            pos = pkd2_resnum_to_pos[chain_id].get(pkd2_resnum)
            if pos is not None:
                pkd2_positions.add(pos)
            for neighbor_atom in neighbors:
                pkd1_resnum = neighbor_atom.get_parent().id[1]
                pkd1_pos = pkd1_resnum_to_pos.get(pkd1_resnum)
                if pkd1_pos is not None:
                    pkd1_positions.add(pkd1_pos)

    print(f"      {len(pkd1_positions)} PKD1 ({PKD1_ACCESSION}) interface residue(s)")
    print(f"      {len(pkd2_positions)} PKD2 ({PKD2_ACCESSION}) interface residue(s)")

    # Bio.Align's alignment.aligned blocks hand back numpy ints, which
    # propagate into these sets via align_to_reference's pos_to_resnum --
    # cast to plain int so json.dumps doesn't choke on numpy.int64.
    payload = {
        PKD1_ACCESSION: sorted(int(p) for p in pkd1_positions),
        PKD2_ACCESSION: sorted(int(p) for p in pkd2_positions),
    }
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
