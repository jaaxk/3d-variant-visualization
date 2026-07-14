"""Structure loading: fetch (network, login node) and load (offline, compute node).

Torch HPC compute nodes may lack outbound internet access, so this module
enforces a hard split:
  - fetch_structure() / fetch_uniprot_reference(): NETWORK. Only ever called
    from `protein-vis fetch`, meant to run on the login node.
  - load_structure(): LOCAL ONLY. Reads from the manifest written by fetch_*;
    raises a clear error if the requested source hasn't been fetched yet.

Structure sources are specified as a short string:
    "alphafold:<UniProt accession>"
    "pdb:<PDB ID>" or "pdb:<PDB ID>:<CHAIN>"
    "file:<local path>"

A structure's own residue numbering is NOT assumed to match the reference
protein's canonical (e.g. UniProt) numbering in general -- experimental
structures routinely have offsets, gaps, engineered constructs, or (as with
BRCA2, whose only PDB structural coverage is the rat ortholog 1IYJ) are a
different organism's ortholog entirely. align_to_reference() performs a real
pairwise sequence alignment to build the correspondence; AlphaFold-sourced
structures are expected to align ~1:1 (verified, not assumed).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

import numpy as np
import requests
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, PDBParser

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"
UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/{accession}.json"
THREEDMOL_JS_URL = "https://3Dmol.org/build/3Dmol-min.js"

# Standard 20 amino acids, 3-letter -> 1-letter. Kept as a local mapping
# rather than relying on Bio.PDB.Polypeptide's three_to_one, whose API has
# changed across biopython versions.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


class StructureFetchError(RuntimeError):
    pass


class StructureLoadError(RuntimeError):
    pass


class AlignmentError(RuntimeError):
    pass


@dataclass(frozen=True)
class StructureSource:
    kind: Literal["alphafold", "pdb", "file"]
    accession_or_path: str
    chain: str | None = None


@dataclass
class StructureData:
    chain_id: str
    resnums: list[int]
    sequence: str
    ca_coords: dict[int, np.ndarray]
    raw_text: str
    fmt: Literal["pdb", "cif"]
    # CA coords/sequence per chain for every protein chain in the file (not
    # just the one selected above) -- used to render a whole-structure
    # "colored by chain" overview for multi-chain files (e.g. a PKD1-PKD2
    # complex). Chains with identical sequences (e.g. the 3 PKD2 copies)
    # are grouped under one color/legend entry rather than one each.
    all_chain_ca_coords: dict[str, dict[int, np.ndarray]] = field(default_factory=dict)
    all_chain_sequences: dict[str, str] = field(default_factory=dict)
    # CA atom B-factor per chain/resnum -- for real experimental structures
    # this is the crystallographic/EM B-factor, but AlphaFold (both the
    # official DB and the AlphaFold Server) writes per-atom pLDDT into this
    # same column, so it doubles as confidence data for AF-derived
    # structures. Callers decide which interpretation applies (see
    # pipeline.py's --confidence flag) -- this dataclass just carries
    # whatever value is actually in the file, unconditionally.
    all_chain_ca_bfactor: dict[str, dict[int, float]] = field(default_factory=dict)


@dataclass
class AlignmentResult:
    pos_to_resnum: dict[int, int] = field(default_factory=dict)
    identity: float = 0.0
    coverage: float = 0.0


def parse_source_spec(spec: str) -> StructureSource:
    if ":" not in spec:
        raise ValueError(
            f"structure spec {spec!r} must be 'alphafold:<ACC>', "
            f"'pdb:<ID>[:<CHAIN>]', or 'file:<path>'"
        )
    kind, _, rest = spec.partition(":")
    if kind == "alphafold":
        return StructureSource(kind="alphafold", accession_or_path=rest)
    if kind == "pdb":
        if ":" in rest:
            pdb_id, chain = rest.split(":", 1)
            return StructureSource(kind="pdb", accession_or_path=pdb_id, chain=chain)
        return StructureSource(kind="pdb", accession_or_path=rest)
    if kind == "file":
        return StructureSource(kind="file", accession_or_path=rest)
    raise ValueError(f"unknown structure source kind {kind!r} in spec {spec!r}")


def _manifest_path(cache_dir: Path) -> Path:
    return cache_dir / "manifest.json"


def _read_manifest(cache_dir: Path) -> dict:
    path = _manifest_path(cache_dir)
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def _write_manifest_entry(cache_dir: Path, key: str, entry: dict) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    manifest = _read_manifest(cache_dir)
    manifest[key] = entry
    _manifest_path(cache_dir).write_text(json.dumps(manifest, indent=2, sort_keys=True))


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _utcnow_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_write(dest: Path, data: bytes) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    tmp.write_bytes(data)
    tmp.replace(dest)


def fetch_structure(spec: str, cache_dir: Path, *, force: bool = False) -> Path:
    """NETWORK. Download and cache a structure file; record it in manifest.json."""
    cache_dir = Path(cache_dir)
    source = parse_source_spec(spec)
    manifest_key = f"{source.kind}:{source.accession_or_path}"

    if source.kind == "file":
        src_path = Path(source.accession_or_path).expanduser().resolve()
        if not src_path.exists():
            raise StructureFetchError(f"file source {src_path} does not exist")
        _write_manifest_entry(
            cache_dir,
            manifest_key,
            {
                "url": None,
                "local_path": str(src_path),
                "fetched_at": _utcnow_iso(),
                "sha256": _sha256(src_path.read_bytes()),
            },
        )
        return src_path

    manifest = _read_manifest(cache_dir)
    if not force and manifest_key in manifest:
        cached_path = cache_dir / manifest[manifest_key]["local_path"]
        if cached_path.exists():
            return cached_path

    if source.kind == "alphafold":
        acc = source.accession_or_path
        resp = requests.get(ALPHAFOLD_API.format(accession=acc), timeout=30)
        resp.raise_for_status()
        data = resp.json()
        if not data:
            raise StructureFetchError(
                f"no AlphaFold DB model exists for accession {acc!r} "
                f"(API returned no entries -- the protein may be too large "
                f"or excluded from AlphaFold DB, as is the case for BRCA2/"
                f"P51587). Use 'pdb:<ID>' or 'file:<path>' instead."
            )
        url = data[0]["pdbUrl"]
        dest = cache_dir / "structures" / f"AF-{acc}-F1.pdb"
    elif source.kind == "pdb":
        pdb_id = source.accession_or_path.upper()
        url = RCSB_DOWNLOAD.format(pdb_id=pdb_id)
        dest = cache_dir / "structures" / f"PDB-{pdb_id}.pdb"
    else:
        raise AssertionError(f"unhandled source kind {source.kind!r}")

    resp = requests.get(url, timeout=60)
    if resp.status_code != 200:
        raise StructureFetchError(
            f"failed to fetch structure from {url} (HTTP {resp.status_code})"
        )
    _atomic_write(dest, resp.content)
    _write_manifest_entry(
        cache_dir,
        manifest_key,
        {
            "url": url,
            "local_path": str(dest.relative_to(cache_dir)),
            "fetched_at": _utcnow_iso(),
            "sha256": _sha256(resp.content),
        },
    )
    return dest


def fetch_uniprot_reference(accession: str, cache_dir: Path) -> Path:
    """NETWORK. Fetch + cache UniProt JSON record and a derived FASTA."""
    cache_dir = Path(cache_dir)
    manifest_key = f"uniprot:{accession}"
    manifest = _read_manifest(cache_dir)
    json_dest = cache_dir / "uniprot" / f"{accession}.json"
    if manifest_key in manifest and json_dest.exists():
        return json_dest

    url = UNIPROT_REST.format(accession=accession)
    resp = requests.get(url, timeout=30)
    if resp.status_code != 200:
        raise StructureFetchError(
            f"failed to fetch UniProt record for {accession} from {url} "
            f"(HTTP {resp.status_code})"
        )
    _atomic_write(json_dest, resp.content)

    data = resp.json()
    sequence = data.get("sequence", {}).get("value", "")
    fasta_dest = cache_dir / "uniprot" / f"{accession}.fasta"
    _atomic_write(
        fasta_dest, f">{accession}\n{sequence}\n".encode()
    )

    _write_manifest_entry(
        cache_dir,
        manifest_key,
        {
            "url": url,
            "local_path": str(json_dest.relative_to(cache_dir)),
            "fasta_path": str(fasta_dest.relative_to(cache_dir)),
            "fetched_at": _utcnow_iso(),
            "sha256": _sha256(resp.content),
        },
    )
    return json_dest


def fetch_3dmol_js(cache_dir: Path) -> Path:
    """NETWORK. Fetch + cache 3Dmol.js so interactive HTML can be self-contained."""
    cache_dir = Path(cache_dir)
    dest = cache_dir / "js" / "3Dmol.min.js"
    if dest.exists():
        return dest
    resp = requests.get(THREEDMOL_JS_URL, timeout=60)
    if resp.status_code != 200:
        raise StructureFetchError(
            f"failed to fetch 3Dmol.js from {THREEDMOL_JS_URL} "
            f"(HTTP {resp.status_code})"
        )
    _atomic_write(dest, resp.content)
    return dest


def load_uniprot_sequence(cache_dir: Path, accession: str) -> str:
    """LOCAL ONLY. Read the canonical sequence cached by fetch_uniprot_reference."""
    fasta_path = Path(cache_dir) / "uniprot" / f"{accession}.fasta"
    if not fasta_path.exists():
        raise StructureLoadError(
            f"no cached UniProt reference for {accession} at {fasta_path} -- "
            f"run `protein-vis fetch --uniprot {accession}` on the login node first"
        )
    lines = fasta_path.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def load_structure(spec: str, cache_dir: Path) -> StructureData:
    """LOCAL ONLY. Load a structure previously cached by fetch_structure()."""
    cache_dir = Path(cache_dir)
    source = parse_source_spec(spec)
    manifest_key = f"{source.kind}:{source.accession_or_path}"
    manifest = _read_manifest(cache_dir)

    if manifest_key not in manifest:
        raise StructureLoadError(
            f"structure source {spec!r} has not been fetched -- run "
            f"`protein-vis fetch --structure {spec}` on the login node first"
        )
    local_path = manifest[manifest_key]["local_path"]
    path = Path(local_path) if Path(local_path).is_absolute() else cache_dir / local_path
    if not path.exists():
        raise StructureLoadError(f"cached structure file missing: {path}")

    fmt: Literal["pdb", "cif"] = "cif" if path.suffix.lower() in (".cif", ".mmcif") else "pdb"
    parser = MMCIFParser(QUIET=True) if fmt == "cif" else PDBParser(QUIET=True)
    parsed = parser.get_structure(source.accession_or_path, str(path))
    model = next(iter(parsed))

    chain = None
    if source.chain:
        try:
            chain = model[source.chain]
        except KeyError:
            raise StructureLoadError(
                f"chain {source.chain!r} not found in {path} "
                f"(available: {[c.id for c in model]})"
            )
    else:
        best = None
        best_len = -1
        for candidate in model:
            n_standard = sum(
                1 for res in candidate if res.id[0] == " " and res.resname in THREE_TO_ONE
            )
            if n_standard > best_len:
                best, best_len = candidate, n_standard
        chain = best

    if chain is None:
        raise StructureLoadError(f"no usable protein chain found in {path}")

    resnums: list[int] = []
    sequence_chars: list[str] = []
    ca_coords: dict[int, np.ndarray] = {}
    for res in chain:
        if res.id[0] != " " or res.resname not in THREE_TO_ONE:
            continue
        resnum = res.id[1]
        resnums.append(resnum)
        sequence_chars.append(THREE_TO_ONE[res.resname])
        if "CA" in res:
            ca_coords[resnum] = np.array(res["CA"].coord, dtype=float)

    all_chain_ca_coords: dict[str, dict[int, np.ndarray]] = {}
    all_chain_sequences: dict[str, str] = {}
    all_chain_ca_bfactor: dict[str, dict[int, float]] = {}
    for candidate in model:
        standard_res = [
            res for res in candidate if res.id[0] == " " and res.resname in THREE_TO_ONE
        ]
        coords = {res.id[1]: np.array(res["CA"].coord, dtype=float) for res in standard_res if "CA" in res}
        if coords:
            all_chain_ca_coords[candidate.id] = coords
            all_chain_sequences[candidate.id] = "".join(THREE_TO_ONE[res.resname] for res in standard_res)
            all_chain_ca_bfactor[candidate.id] = {
                res.id[1]: float(res["CA"].get_bfactor()) for res in standard_res if "CA" in res
            }

    return StructureData(
        chain_id=chain.id,
        resnums=resnums,
        sequence="".join(sequence_chars),
        ca_coords=ca_coords,
        raw_text=path.read_text(),
        fmt=fmt,
        all_chain_ca_coords=all_chain_ca_coords,
        all_chain_sequences=all_chain_sequences,
        all_chain_ca_bfactor=all_chain_ca_bfactor,
    )


def align_to_reference(
    struct: StructureData, reference_seq: str, *, min_identity: float = 0.4
) -> AlignmentResult:
    """Pairwise local alignment of the structure's sequence to a reference sequence.

    Builds a UniProt-position -> structure-resnum map. Required in general
    (numbering offsets, gaps, engineered constructs, ortholog structures are
    all common) -- not a special case for any one protein.
    """
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(struct.sequence, reference_seq)
    best = alignments[0]

    struct_blocks, ref_blocks = best.aligned
    pos_to_resnum: dict[int, int] = {}
    n_matches = 0
    n_aligned = 0
    for (s_start, s_end), (r_start, r_end) in zip(struct_blocks, ref_blocks):
        block_len = s_end - s_start
        for i in range(block_len):
            struct_idx = s_start + i
            ref_idx = r_start + i
            ref_pos = ref_idx + 1  # 1-based UniProt numbering
            pos_to_resnum[ref_pos] = struct.resnums[struct_idx]
            n_aligned += 1
            if struct.sequence[struct_idx] == reference_seq[ref_idx]:
                n_matches += 1

    identity = n_matches / n_aligned if n_aligned else 0.0
    coverage = n_aligned / len(reference_seq) if reference_seq else 0.0

    if identity < min_identity:
        raise AlignmentError(
            f"alignment identity ({identity:.1%}) is below the minimum "
            f"threshold ({min_identity:.1%}) -- this structure is likely the "
            f"wrong protein/chain, or min_identity needs tuning"
        )

    return AlignmentResult(pos_to_resnum=pos_to_resnum, identity=identity, coverage=coverage)
