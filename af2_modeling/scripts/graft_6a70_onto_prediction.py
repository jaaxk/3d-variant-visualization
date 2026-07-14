#!/usr/bin/env python3
"""Graft a template structure's real coordinates onto the overlapping region
of an AlphaFold2 prediction -- and, via chaining, onto the output of a prior
run of this same script.

Why: AlphaFold's template mechanism feeds a template's structure into the
network as a learned input feature -- it never copies template coordinates
directly. So even with a template found during modeling, the predicted
region overlapping it will be *influenced by* but not identical to its real
coordinates. This script makes it identical, as a post-processing step that
doesn't depend on which templates AlphaFold's own search happened to use:

1. Match each chain in the prediction to its corresponding chain in the
   template by sequence (never assume chain-letter conventions match between
   the two files -- same principle protein_vis itself follows throughout).
2. Align each matched pair's sequences (BLOSUM62, same method as
   protein_vis.structure.align_to_reference) to find the shared residue
   range, independent of either file's internal numbering.
3. Compute ONE rigid-body superposition (Bio.PDB.Superimposer) from ALL
   matched/shared CA atoms across ALL chains combined, and apply it to the
   entire predicted structure -- so the whole assembly moves into the
   template's reference frame together, preserving whatever inter-chain
   packing the prediction produced for the new region.
4. Splice: keep the template's real atoms for shared residues, keep the
   (transformed) prediction's atoms for new residues. This is a hard splice
   -- it may leave a small bond-length/angle discontinuity right at the
   junction, since the prediction's coordinates there won't land exactly on
   the template's. A rigorous fix (e.g. a short OpenMM energy minimization
   localized to the junction residues) is a documented follow-up in the
   README, not built here -- not worth introducing a whole MD dependency for
   what's likely a sub-angstrom kink at one or two residues, until/unless it
   turns out to actually matter for downstream use.
5. Optionally (--include-unmatched-template-chains) carry through, unchanged,
   any template chain that no predicted chain matched -- e.g. when grafting a
   PKD1-only monomer prediction onto an already-grafted PKD1/PKD2 complex,
   the complex's 3 PKD2 chains have no counterpart in the monomer prediction
   at all, but still belong in the final merged structure.

Chaining across multiple graft stages: each run can optionally read the
*previous* stage's provenance sidecar (--template-provenance) to know what
--template's own residues already are (e.g. "6A70" vs "AlphaFold2: Complex"),
so a matched residue's output label is inherited from that prior stage
rather than collapsed back to a single generic "template" bucket. This is
how a chain of grafts accumulates a full N-way per-residue provenance map
across every stage, not just a binary this-run-vs-that-run split.

Usage (stage 1 -- template is a real deposited structure, e.g. 6A70):
    python graft_6a70_onto_prediction.py \\
        --predicted /scratch/.../ranked_0.pdb \\
        --template /scratch/jv2807/pkd1/structure_cache/structures/PDB-6A70.pdb \\
        --output /scratch/.../pkd1_pkd2_grafted.pdb \\
        --provenance-output /scratch/.../pkd1_pkd2_grafted.provenance.json \\
        --new-label "AlphaFold2: Complex"

Usage (stage 2 -- template is stage 1's own output, chaining provenance):
    python graft_6a70_onto_prediction.py \\
        --predicted /scratch/.../fold_pkd1_human_model_0.cif \\
        --template /scratch/.../pkd1_pkd2_grafted.pdb \\
        --template-provenance /scratch/.../pkd1_pkd2_grafted.provenance.json \\
        --new-label "AlphaFold2: PKD1 monomer" \\
        --include-unmatched-template-chains \\
        --output /scratch/.../pkd1_full_pkd2_grafted.pdb \\
        --provenance-output /scratch/.../pkd1_full_pkd2_grafted.provenance.json
"""

from __future__ import annotations

import argparse
import json
import string
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


def assign_passthrough_ids(
    existing_ids: set,
    template_chains: list,
    assigned_tmpl_chain_ids: set,
) -> dict:
    """Decide the final output chain id for every --template chain that no
    --predicted chain matched (the "passthrough" chains -- e.g. the 3 PKD2
    chains when grafting a PKD1-only monomer onto an already-grafted
    complex). Matched --predicted chains are never renamed by this script;
    passthrough chains are renamed only if their original id collides with
    an id already in use, via a single running `used` set updated
    immediately after each assignment (so a later passthrough chain's
    collision check sees earlier renames, not just the original id set).
    Iterates template_chains in their original file order for a
    deterministic result regardless of dict/set iteration order elsewhere.
    """
    used = set(existing_ids)
    pool = list(string.ascii_uppercase) + [str(d) for d in range(10)]
    mapping = {}
    for tmpl_info in template_chains:
        orig_id = tmpl_info.chain.id
        if orig_id in assigned_tmpl_chain_ids:
            continue
        if orig_id not in used:
            new_id = orig_id
        else:
            new_id = next(c for c in pool if c not in used)
            print(f"  passthrough template chain {orig_id!r} renamed -> {new_id!r} (id collision)")
        used.add(new_id)
        mapping[orig_id] = new_id
    return mapping


def build_provenance_labels(
    predicted_chains: list,
    per_chain_matches: dict,
    template_provenance: dict | None,
    default_label: str,
    new_label: str,
    passthrough_chains: dict,
) -> dict:
    """{output_chain_id: {resnum_str: label}} for every residue in the final
    output structure -- both the (possibly partially-spliced) predicted
    chains and any passthrough template chains.

    For a predicted chain's residue matched to the template: label is
    inherited from `template_provenance[tmpl_chain_id][str(tmpl_resnum)]` if
    that sidecar was given and has an entry, else `default_label` (the
    template genuinely IS the original source, e.g. 6A70 itself, when no
    prior-stage provenance exists). For a predicted chain's residue with no
    template match at all: `new_label` (genuinely new this stage).

    For a passthrough chain (no predicted counterpart): every residue's
    label is copied verbatim from `template_provenance` for that chain's
    *original* template id, or `default_label` for all of it if no
    `template_provenance` was given.

    Pure function of already-computed data -- no alignment/matching here.
    """
    template_provenance = template_provenance or {}
    labels: dict = {}

    for pred_info in predicted_chains:
        resnum_map = per_chain_matches.get(pred_info.chain.id, {})
        chain_labels = {}
        for resnum in pred_info.resnums:
            if resnum in resnum_map:
                tmpl_chain_id, tmpl_resnum = resnum_map[resnum]
                label = template_provenance.get(tmpl_chain_id, {}).get(str(tmpl_resnum), default_label)
            else:
                label = new_label
            chain_labels[str(resnum)] = label
        labels[pred_info.chain.id] = chain_labels

    for new_id, (orig_tmpl_id, tmpl_info) in passthrough_chains.items():
        tmpl_prov = template_provenance.get(orig_tmpl_id, {})
        labels[new_id] = {
            str(resnum): tmpl_prov.get(str(resnum), default_label)
            for resnum in tmpl_info.resnums
        }

    return labels


def write_provenance(labels: dict, out_path: str) -> None:
    """Serialize the {chain_id: {resnum_str: label}} provenance map and log
    per-label residue counts."""
    Path(out_path).write_text(json.dumps(labels, indent=2, sort_keys=True))
    counts: dict = {}
    for chain_labels in labels.values():
        for label in chain_labels.values():
            counts[label] = counts.get(label, 0) + 1
    n_total = sum(counts.values())
    counts_str = ", ".join(f"{label}={n}" for label, n in sorted(counts.items()))
    print(f"wrote provenance ({n_total} residue(s) across {len(labels)} chain(s): {counts_str}) -> {out_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predicted", required=True,
                     help="AlphaFold2 prediction (PDB/mmCIF) -- or, with --provenance-only, "
                          "an already-grafted structure to compute provenance for.")
    ap.add_argument("--template", required=True,
                     help="Structure to splice real coordinates from (PDB/mmCIF) -- either a "
                          "deposited structure (e.g. 6A70) or a prior stage's own graft output.")
    ap.add_argument("--output", required=False, help="Grafted structure output path.")
    ap.add_argument("--min-identity", type=float, default=0.9,
                     help="Minimum identity to accept a chain match (default 0.9 -- these "
                          "are the same real protein, near-identical sequence expected).")
    ap.add_argument("--template-provenance",
                     help="Provenance JSON sidecar (see --provenance-output) describing "
                          "--template's OWN per-residue provenance, from a prior stage's graft. "
                          "A matched residue's output label is inherited from this file "
                          "(falling back to --default-label if the specific residue is "
                          "missing); omit this when --template is a real deposited structure "
                          "with no prior grafting stage behind it.")
    ap.add_argument("--default-label", default="6A70",
                     help="Label for residues spliced from --template when --template-provenance "
                          "wasn't given (i.e. --template genuinely IS the original source, e.g. "
                          "6A70 itself). Default: %(default)s")
    ap.add_argument("--new-label", default="AlphaFold2",
                     help="Label for residues in --predicted with no counterpart in --template "
                          "at all (genuinely new this stage). Default: %(default)s")
    ap.add_argument("--include-unmatched-template-chains", action="store_true",
                     help="Carry through, unchanged, any --template chain that no --predicted "
                          "chain matched (e.g. the 3 PKD2 chains when --predicted is a PKD1-only "
                          "monomer). These chains are already in the correct reference frame, "
                          "so no transform is applied. No-ops (with a printed note) under "
                          "--provenance-only, since re-matching an already-merged structure by "
                          "sequence naturally re-discovers these as ordinary matched chains.")
    ap.add_argument("--provenance-output",
                     help="Also write a JSON sidecar ({chain_id: {resnum_str: label}}) recording "
                          "each residue's provenance label.")
    ap.add_argument("--provenance-only", action="store_true",
                     help="Skip superposition/splice entirely. Just compute which residues in "
                          "--predicted have a sequence-matched (= template-resolved-density-"
                          "covered) counterpart in --template and write --provenance-output. "
                          "Use this against an ALREADY-GRAFTED structure to get provenance "
                          "without re-running the graft -- the same chain-matching the graft "
                          "itself uses already defines 'covered by --template' exactly, so "
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

    template_provenance = None
    if args.template_provenance:
        template_provenance = json.loads(Path(args.template_provenance).read_text())
        if any(isinstance(v, list) for v in template_provenance.values()):
            raise SystemExit(
                f"{args.template_provenance}: old-format (flat list) provenance JSON -- "
                "regenerate it with this script's current version first"
            )

    assignment, per_chain_matches = compute_chain_matches(
        predicted_chains, template_chains, args.min_identity
    )

    if args.provenance_only:
        if not args.provenance_output:
            raise SystemExit("--provenance-only requires --provenance-output")
        if args.include_unmatched_template_chains:
            print("  --include-unmatched-template-chains is a no-op under --provenance-only "
                  "(re-matching an already-merged structure by sequence already finds these "
                  "chains as ordinary matches)")
        labels = build_provenance_labels(
            predicted_chains, per_chain_matches, template_provenance,
            args.default_label, args.new_label, passthrough_chains={},
        )
        write_provenance(labels, args.provenance_output)
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

    # 5. Optionally carry through template chains with no predicted counterpart at all,
    # unchanged (already in the correct reference frame -- no transform needed).
    passthrough_chains = {}
    if args.include_unmatched_template_chains:
        assigned_tmpl_ids = {match.chain.id for match, _, _ in assignment.values()}
        existing_ids = {c.chain.id for c in predicted_chains}
        id_mapping = assign_passthrough_ids(existing_ids, template_chains, assigned_tmpl_ids)
        for orig_id, new_id in id_mapping.items():
            tmpl_info = next(c for c in template_chains if c.chain.id == orig_id)
            new_chain = tmpl_info.chain.copy()
            new_chain.id = new_id
            predicted_model.add(new_chain)
            passthrough_chains[new_id] = (orig_id, tmpl_info)
        if passthrough_chains:
            print(f"carried through {len(passthrough_chains)} unmatched template chain(s): "
                  f"{list(passthrough_chains)}")

    io = PDBIO()
    io.set_structure(predicted_model)
    io.save(args.output)
    print(f"wrote {args.output}")

    if args.provenance_output:
        labels = build_provenance_labels(
            predicted_chains, per_chain_matches, template_provenance,
            args.default_label, args.new_label, passthrough_chains,
        )
        write_provenance(labels, args.provenance_output)


if __name__ == "__main__":
    main()
