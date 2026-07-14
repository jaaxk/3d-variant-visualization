#!/usr/bin/env python3
"""Graft 6A70's real deposited coordinates onto the overlapping region of a
new AlphaFold2 PKD1/PKD2 complex prediction.

Why: AlphaFold's template mechanism feeds a template's structure into the
network as a learned input feature -- it never copies template coordinates
directly. So even with 6A70 found as a template during modeling, the
predicted region overlapping 6A70 will be *influenced by* but not identical
to 6A70's real coordinates. This script makes it identical, as a
post-processing step that doesn't depend on which templates AlphaFold's own
search happened to use:

1. Match each chain in the prediction to its corresponding chain in 6A70 by
   sequence (never assume chain-letter conventions match between the two
   files -- same principle protein_vis itself follows throughout).
2. Align each matched pair's sequences (BLOSUM62, same method as
   protein_vis.structure.align_to_reference) to find the shared residue
   range, independent of either file's internal numbering.
3. Compute ONE rigid-body superposition (Bio.PDB.Superimposer) from ALL
   matched/shared CA atoms across ALL chains combined, and apply it to the
   entire predicted structure -- so the whole assembly moves into 6A70's
   reference frame together, preserving whatever inter-chain packing the
   prediction produced for the new region.
4. Splice: keep 6A70's real atoms for shared residues, keep the
   (transformed) prediction's atoms for new residues. This is a hard splice
   -- it may leave a small bond-length/angle discontinuity right at the
   junction, since the prediction's coordinates there won't land exactly on
   6A70's. A rigorous fix (e.g. a short OpenMM energy minimization localized
   to the junction residues) is a documented follow-up in the README, not
   built here -- not worth introducing a whole MD dependency for what's
   likely a sub-angstrom kink at one or two residues, until/unless it turns
   out to actually matter for downstream use.

Usage:
    python graft_6a70_onto_prediction.py \\
        --predicted /scratch/.../ranked_0.pdb \\
        --template /scratch/jv2807/pkd1/structure_cache/structures/PDB-6A70.pdb \\
        --output /scratch/.../pkd1_pkd2_grafted.pdb
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Superimposer
from Bio.PDB.Chain import Chain

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

@dataclass
class ChainInfo:
    chain: Chain
    resnums: list
    sequence: str


def parse_structure(path: str):
    fmt = "cif" if path.lower().endswith((".cif", ".mmcif")) else "pdb"
    parser = MMCIFParser(QUIET=True) if fmt == "cif" else PDBParser(QUIET=True)
    structure = parser.get_structure("s", path)
    return next(iter(structure))


def chain_info(chain: Chain) -> ChainInfo:
    resnums, seq_chars = [], []
    for res in chain:
        if res.id[0] != " " or res.resname not in THREE_TO_ONE:
            continue
        resnums.append(res.id[1])
        seq_chars.append(THREE_TO_ONE[res.resname])
    return ChainInfo(chain=chain, resnums=resnums, sequence="".join(seq_chars))


def align_sequences(seq_a: str, seq_b: str):
    """BLOSUM62 local alignment -- same method as
    protein_vis.structure.align_to_reference. Returns (identity, list of
    (index_in_a, index_in_b) for exactly-matching aligned positions)."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(seq_a, seq_b)[0]
    a_blocks, b_blocks = alignment.aligned
    matches = []
    n_aligned = n_matches = 0
    for (a_start, a_end), (b_start, b_end) in zip(a_blocks, b_blocks):
        for i in range(a_end - a_start):
            ai, bi = a_start + i, b_start + i
            n_aligned += 1
            if seq_a[ai] == seq_b[bi]:
                n_matches += 1
                matches.append((ai, bi))
    identity = n_matches / n_aligned if n_aligned else 0.0
    return identity, matches


def matched_ca_pairs(pred_info: ChainInfo, tmpl_info: ChainInfo, matches) -> list[tuple]:
    pairs = []
    for pred_idx, tmpl_idx in matches:
        pred_res = pred_info.chain[pred_info.resnums[pred_idx]]
        tmpl_res = tmpl_info.chain[tmpl_info.resnums[tmpl_idx]]
        if "CA" in pred_res and "CA" in tmpl_res:
            pairs.append((pred_res["CA"], tmpl_res["CA"], pred_info.resnums[pred_idx], tmpl_info.resnums[tmpl_idx]))
    return pairs


def assign_chains(predicted_chains: list[ChainInfo], template_chains: list[ChainInfo], min_identity: float):
    """One-to-one assignment of predicted chains to template chains.

    Sequence identity alone can't disambiguate *identical* repeated chains
    (e.g. the 3 PKD2 copies) -- every copy scores 100% against every other
    copy's template counterpart. So: group predicted chains by identical
    sequence, and within a group of size >1, brute-force every permutation
    of assigning that group's template candidates and keep whichever
    permutation gives the lowest trial superposition RMSD (a geometric
    criterion, since sequence can't break the tie). Groups are small here
    (at most 3 PKD2 copies), so brute force is trivial.
    """
    import itertools

    # All (predicted, template) pairs with acceptable identity, plus their matched atoms.
    candidates: dict[int, list[tuple]] = {}  # id(pred) -> [(tmpl_info, identity, matches), ...]
    for pred_info in predicted_chains:
        scored = []
        for tmpl_info in template_chains:
            identity, matches = align_sequences(pred_info.sequence, tmpl_info.sequence)
            if identity >= min_identity and matches:
                scored.append((tmpl_info, identity, matches))
        candidates[id(pred_info)] = scored

    # Group predicted chains by identical sequence (the repeated-chain case).
    groups: dict[str, list[ChainInfo]] = {}
    for pred_info in predicted_chains:
        groups.setdefault(pred_info.sequence, []).append(pred_info)

    assignment: dict[str, tuple] = {}  # predicted chain id -> (tmpl_info, identity, matches)
    for group in groups.values():
        if len(group) == 1:
            pred_info = group[0]
            scored = candidates[id(pred_info)]
            if scored:
                best = max(scored, key=lambda s: len(s[2]))
                assignment[pred_info.chain.id] = best
            continue

        # Repeated chains: try every permutation of template candidates shared by
        # (at least one of) the group, pick the assignment with lowest combined RMSD.
        shared_candidates = {t.chain.id: t for pred_info in group for t, _, _ in candidates[id(pred_info)]}
        template_options = list(shared_candidates.values())
        best_perm, best_rms = None, float("inf")
        for perm in itertools.permutations(template_options, min(len(group), len(template_options))):
            fixed, moving = [], []
            per_pred = {}
            for pred_info, tmpl_info in zip(group, perm):
                identity, matches = align_sequences(pred_info.sequence, tmpl_info.sequence)
                pairs = matched_ca_pairs(pred_info, tmpl_info, matches)
                if not pairs:
                    break
                per_pred[pred_info.chain.id] = (tmpl_info, identity, matches)
                fixed.extend(p[1] for p in pairs)
                moving.extend(p[0] for p in pairs)
            else:
                sup = Superimposer()
                sup.set_atoms(fixed, moving)
                if sup.rms < best_rms:
                    best_rms, best_perm = sup.rms, per_pred
        if best_perm:
            assignment.update(best_perm)

    return assignment


def compute_chain_matches(predicted_chains, template_chains, min_identity: float):
    """One-to-one chain assignment (see assign_chains) plus, for every
    assigned chain, the per-residue match map -- predicted resnum -> (matched
    template chain id, template resnum) -- for every residue whose sequence
    aligned to a residue *actually present with coordinates* in the template
    (chain_info() only ever extracts residues with modeled coordinates, so
    "sequence-matched" and "covered by the template's resolved density" are
    the same thing here).

    Shared by both the graft path (which additionally pulls CA atoms from
    this map for superposition) and the read-only provenance path (which
    only needs the map itself, to know which residues are EM vs. AF-only --
    see --provenance-only below). Kept separate from the CA-atom/superposition
    step so provenance can be computed without redoing any of that.
    """
    assignment = assign_chains(predicted_chains, template_chains, min_identity)
    per_chain_matches: dict[str, dict[int, tuple]] = {}
    for pred_info in predicted_chains:
        assigned = assignment.get(pred_info.chain.id)
        if assigned is None:
            print(f"  chain {pred_info.chain.id}: no confident template match -- "
                  f"treated as fully new, no grafting")
            per_chain_matches[pred_info.chain.id] = {}
            continue
        match, identity, matches = assigned
        print(f"  chain {pred_info.chain.id} ({len(pred_info.sequence)} aa) -> "
              f"template chain {match.chain.id} ({len(match.sequence)} aa), "
              f"identity={identity:.1%}, {len(matches)} shared residue(s)")

        resnum_map = {}
        for pred_idx, tmpl_idx in matches:
            pred_resnum = pred_info.resnums[pred_idx]
            tmpl_resnum = match.resnums[tmpl_idx]
            resnum_map[pred_resnum] = (match.chain.id, tmpl_resnum)
        per_chain_matches[pred_info.chain.id] = resnum_map
    return assignment, per_chain_matches


def write_provenance(per_chain_matches: dict, out_path: str) -> None:
    """{predicted_chain_id: [resnum, ...]} of every residue covered by the
    template ("EM"); any resnum absent from its chain's list is implicitly
    AF-predicted -- consumers never need a second "AF" list."""
    provenance = {
        chain_id: sorted(resnum_map)
        for chain_id, resnum_map in per_chain_matches.items()
    }
    Path(out_path).write_text(json.dumps(provenance, indent=2, sort_keys=True))
    n_em = sum(len(v) for v in provenance.values())
    print(f"wrote provenance ({n_em} EM residue(s) across {len(provenance)} chain(s)) -> {out_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predicted", required=True,
                     help="AlphaFold2 prediction (PDB/mmCIF) -- or, with --provenance-only, "
                          "an already-grafted structure to compute provenance for.")
    ap.add_argument("--template", required=True, help="6A70 structure (PDB/mmCIF)")
    ap.add_argument("--output", required=False, help="Grafted structure output path.")
    ap.add_argument("--min-identity", type=float, default=0.9,
                     help="Minimum identity to accept a chain match (default 0.9 -- these "
                          "are the same real protein, near-identical sequence expected).")
    ap.add_argument("--provenance-output",
                     help="Also write a JSON sidecar (chain_id -> [EM resnum, ...]) recording "
                          "which residues came directly from --template.")
    ap.add_argument("--provenance-only", action="store_true",
                     help="Skip superposition/splice entirely. Just compute which residues in "
                          "--predicted have a sequence-matched (= template-resolved-density-"
                          "covered) counterpart in --template and write --provenance-output. "
                          "Use this against an ALREADY-GRAFTED structure to get EM-vs-AF "
                          "provenance without re-running the graft -- the same chain-matching "
                          "the graft itself uses already defines 'covered by 6A70' exactly, so "
                          "nothing needs to be re-modeled, only recomputed against the files "
                          "already on disk.")
    args = ap.parse_args()

    predicted_model = parse_structure(args.predicted)
    template_model = parse_structure(args.template)

    predicted_chains = [chain_info(c) for c in predicted_model if len(chain_info(c).sequence) > 0]
    template_chains = [chain_info(c) for c in template_model if len(chain_info(c).sequence) > 0]

    print(f"predicted structure: {len(predicted_chains)} chain(s) "
          f"({[c.chain.id for c in predicted_chains]})")
    print(f"template structure:  {len(template_chains)} chain(s) "
          f"({[c.chain.id for c in template_chains]})")

    assignment, per_chain_matches = compute_chain_matches(
        predicted_chains, template_chains, args.min_identity
    )

    if args.provenance_only:
        if not args.provenance_output:
            raise SystemExit("--provenance-only requires --provenance-output")
        write_provenance(per_chain_matches, args.provenance_output)
        return

    if not args.output:
        raise SystemExit("--output is required unless --provenance-only")

    # 3. One global rigid-body superposition from all matched CA atoms (pulled from
    # per_chain_matches, combined across ALL chains), applied to every atom in the
    # predicted structure.
    fixed_atoms, moving_atoms = [], []
    for pred_info in predicted_chains:
        resnum_map = per_chain_matches[pred_info.chain.id]
        for pred_resnum, (tmpl_chain_id, tmpl_resnum) in resnum_map.items():
            tmpl_chain = next(c for c in template_model if c.id == tmpl_chain_id)
            pred_res, tmpl_res = pred_info.chain[pred_resnum], tmpl_chain[tmpl_resnum]
            if "CA" in pred_res and "CA" in tmpl_res:
                moving_atoms.append(pred_res["CA"])
                fixed_atoms.append(tmpl_res["CA"])

    if not fixed_atoms:
        raise SystemExit("no chain matched the template with sufficient identity -- nothing to graft")

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    print(f"superposition RMSD over {len(fixed_atoms)} shared CA atom(s): {sup.rms:.3f} A")
    all_predicted_atoms = [atom for chain in predicted_model for res in chain for atom in res]
    sup.apply(all_predicted_atoms)

    # 4. Splice: overwrite each shared residue's atoms with the template's real
    # coordinates (matching atom names only -- side chains can differ in atom count/
    # rotamer, backbone N/CA/C/O always match). Residues with no template counterpart
    # keep the (already-transformed) predicted coordinates untouched.
    grafted_residue_count = 0
    for pred_info in predicted_chains:
        resnum_map = per_chain_matches[pred_info.chain.id]
        for pred_resnum, (tmpl_chain_id, tmpl_resnum) in resnum_map.items():
            tmpl_chain = next(c for c in template_model if c.id == tmpl_chain_id)
            if tmpl_resnum not in tmpl_chain or pred_resnum not in pred_info.chain:
                continue
            pred_res, tmpl_res = pred_info.chain[pred_resnum], tmpl_chain[tmpl_resnum]
            for atom in pred_res:
                if atom.get_name() in tmpl_res:
                    atom.set_coord(tmpl_res[atom.get_name()].get_coord())
            grafted_residue_count += 1

    print(f"grafted {grafted_residue_count} residue(s) directly from the template")

    io = PDBIO()
    io.set_structure(predicted_model)
    io.save(args.output)
    print(f"wrote {args.output}")

    if args.provenance_output:
        write_provenance(per_chain_matches, args.provenance_output)


if __name__ == "__main__":
    main()
