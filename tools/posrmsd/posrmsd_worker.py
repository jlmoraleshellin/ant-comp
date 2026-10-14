#!/usr/bin/env python3
"""
CA RMSD between two structures over an EXPLICIT list of chain positions.

PREMISE -- read this before using any number below.

Structure A (the design / prediction, ``--col-a``) and structure B (the reference
/ backbone, ``--col-b``) are assumed to carry **the same residue ordering** on the
two chains being compared: position N of A's chain is the same residue as position
N of B's chain. That holds by construction when A is a structure prediction of a
sequence designed on B's backbone (Boltz of a ProteinMPNN sequence vs the rfd3
backbone chain it was designed on): same length, same order, no insertions.

Given that, the residue correspondence is just "the same index in both files", and
this tool never derives or searches for one. It:

1. selects ONE chain in each structure (``--chain-a`` / ``--chain-b``; when omitted
   the structure must hold exactly one polymer chain -- a chain is never guessed);
2. takes the Nth residue of that chain's polymer, **1-indexed, in FILE ORDER** --
   the Nth polymer residue as written, NOT the residue with PDB number N;
3. superposes A's CA atoms at those positions ONTO B's CA atoms at the same
   positions (SVD/Kabsch, reflections excluded) -- B is the fixed reference frame,
   A is moved;
4. reports the CA RMSD after that one fit, the largest single-CA deviation, and a
   per-segment RMSD for every comma-separated fragment of the position spec, ALL
   from the same single global fit (never a per-segment refit).

WHEN THE NUMBERS MEAN NOTHING

* A and B do not share the residue ordering (an insertion/deletion between them,
  a different construct, a domain-swapped chain, the wrong chain picked). The tool
  will still compare "position N with position N" and return a plausible RMSD for
  the wrong residue pairs. The trust columns below are what catch this --
  ``seq_match_frac`` above all; ``a_len``/``b_len`` differing sets a ``warn:``
  status.
* The position spec was written against a different numbering than file order
  (e.g. PDB residue numbers of a chain not numbered 1..N). See ``*_resnum_offset``
  and ``*_gapless``.

It does NOT answer: whether the rest of the chain folded (that is ``usalign``),
whether the positions are the ones you meant (the spec is taken as given), or
anything about side chains (CA only).

TRUST COLUMNS (read before the RMSD):
    a_chain, b_chain        chain actually indexed in each structure
    a_len, b_len            polymer residues on those chains; differing -> warn:
    n_positions             positions compared (every requested one, or an error)
    seq_match_frac          fraction of positions whose one-letter residue identity
                            is equal in A and B. 'X' (unknown) never counts as a
                            match, so a poly-UNK file cannot confirm anything
    a_motif_seq/b_motif_seq one-letter residues actually read, segments
                            comma-joined -- shows WHICH positions disagree
    a_resnum_offset         residue number minus file position, when constant over
    b_resnum_offset         the selected positions; NA when not constant
    a_gapless, b_gapless    1 = the whole chain is numbered consecutively with no
                            gaps and no insertion codes, 0 otherwise

Every refusal is loud: a missing chain, an ambiguous chain, a position beyond a
chain's end, a selected residue without a CA, fewer than 3 positions -- each is an
``error:`` status naming the structure, the chain and the index, with NA for every
numeric column. Never an RMSD over a quietly smaller set.

Writes a one-row TSV (name, status, sup_path, metrics...) and structure A moved
into B's frame as mmCIF. Normally invoked per task by ``posrmsd.sh``.

Pure gemmi + numpy; no prosapia import (the task image does not carry it). The
position spec arrives ALREADY PARSED, as JSON, from the manifest builder, which
parsed it with prosapia's own position parser -- there is no second grammar here.

Usage (standalone):
    python posrmsd_worker.py --name d0 --a pred.cif --b backbone.pdb \\
        --chain-a A --chain-b D \\
        --segments '{"labels":["11:13"],"segments":[[11,12,13]]}' \\
        --out-struct d0_sup.cif --result-tsv d0.tsv
"""

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np

# Kabsch needs three non-collinear points; below that the fit is meaningless.
MIN_POSITIONS = 3

# Metric columns of the per-design TSV, in order. Bare names; the collector adds
# the comparison prefix and the driver adds the tool leaf. The per-segment
# seg<i>_rmsd columns are NOT here -- their count is per design.
METRIC_COLUMNS = [
    # --- trust: read these before any RMSD below ---
    "a_chain",
    "b_chain",
    "a_len",
    "b_len",
    "n_positions",
    "seq_match_frac",
    "a_motif_seq",
    "b_motif_seq",
    "a_resnum_offset",
    "b_resnum_offset",
    "a_gapless",
    "b_gapless",
    # --- science ---
    "rmsd",
    "max_dev",
    "n_segments",
    "worst_seg_rmsd",
    "seg_ranges",
    "seg_rmsds",
]
RESULT_COLUMNS = ["name", "status", "sup_path", *METRIC_COLUMNS]

PEPTIDE_TYPES = {gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD}


# ---------------------------------------------------------------------------
# Structure bookkeeping
# ---------------------------------------------------------------------------


def one_letter(res: gemmi.Residue) -> str:
    """One-letter code, 'X' when gemmi does not know it. Non-standard residues
    gemmi tabulates in lower case (MSE -> 'm') are upper-cased (-> 'M')."""
    info = gemmi.find_tabulated_residue(res.name)
    code = (info.one_letter_code if info else "").upper()
    return code if len(code) == 1 and code.isalpha() else "X"


def read_structure(path: Path, what: str) -> gemmi.Structure:
    """Parse a PDB/mmCIF (optionally .gz) or raise with a message naming *what*."""
    if str(path) in ("", "."):
        raise ValueError(f"{what}: no path given")
    if not path.is_file():
        raise ValueError(f"{what} not found: {path}")
    try:
        st = gemmi.read_structure(str(path))
    except Exception as e:  # noqa: BLE001 - reported as data by the caller
        raise ValueError(f"{what} {path} could not be parsed: {e}") from e
    if len(st) == 0:
        raise ValueError(f"{what} {path} holds no model")
    # Deterministic CA pick when a crystal structure carries altlocs.
    st.remove_alternative_conformations()
    # A PDB chain letter that reappears after other chains (e.g. its waters
    # listed at the end) must be ONE chain. gemmi 0.7.5 already merges such parts
    # on read; this is an idempotent guard for other versions. Polymer residues
    # keep their file order.
    st.merge_chain_parts()
    st.setup_entities()
    return st


@dataclass
class SelectedChain:
    name: str
    residues: list[gemmi.Residue]  # the chain's polymer, in file order


def select_chain(
    st: gemmi.Structure, requested: str, what: str, flag: str
) -> SelectedChain:
    """The ONE chain to index into, from the first model.

    ``requested`` empty -> the structure must hold exactly one polymer chain.
    More than one is an error listing them: silently taking "the first chain" is
    how a comparison ends up target-vs-target (see the chainsel skill on
    ``usalign --ter 2``).
    """
    model = st[0]
    names = [c.name for c in model]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        raise ValueError(
            f"{what} has more than one chain named {', '.join(dupes)} after "
            f"merging chain parts; cannot tell which one is meant"
        )
    polymers = [(c.name, c.get_polymer()) for c in model]
    polymers = [(n, p) for n, p in polymers if len(p) > 0]

    if requested:
        found = [(n, p) for n, p in polymers if n == requested]
        if not found:
            if requested in names:
                raise ValueError(
                    f"{what}: chain {requested} ({flag}) has no polymer residues "
                    f"(polymer chains present: {', '.join(n for n, _ in polymers) or 'none'})"
                )
            raise ValueError(
                f"{what}: chain {requested} ({flag}) not found (chains present: "
                f"{', '.join(names) or 'none'})"
            )
        name, polymer = found[0]
    else:
        if not polymers:
            raise ValueError(f"{what} holds no polymer chain")
        if len(polymers) > 1:
            raise ValueError(
                f"{what} has {len(polymers)} polymer chains "
                f"({', '.join(n for n, _ in polymers)}) and {flag} was not given; "
                f"pass {flag} to choose one -- a chain is never guessed"
            )
        name, polymer = polymers[0]

    ptype = polymer.check_polymer_type()
    if ptype not in PEPTIDE_TYPES:
        raise ValueError(
            f"{what}: chain {name} is not a protein chain (gemmi polymer type "
            f"{ptype.name}); CA positions are undefined"
        )
    return SelectedChain(name, list(polymer))


def gapless(residues: list[gemmi.Residue]) -> bool:
    """True when residue numbers run first, first+1, ... with no insertion codes --
    the case in which file position, (resnum - first + 1) and ProteinMPNN's
    position all coincide."""
    first = residues[0].seqid.num
    return all(
        r.seqid.num == first + i and r.seqid.icode in (" ", "")
        for i, r in enumerate(residues)
    )


@dataclass
class Picked:
    xyz: np.ndarray  # (N, 3) CA coordinates, spec order
    seq: list[str]  # one-letter per position, spec order
    offset: int | None  # resnum - position, when constant over the picks


def pick_positions(chain: SelectedChain, positions: list[int], what: str) -> Picked:
    """CA atom and identity of every requested 1-indexed FILE-ORDER position.

    Any position beyond the chain's end, or a selected residue without a CA, is
    an error -- reported for all offending positions at once, never skipped.
    """
    n = len(chain.residues)
    beyond = [p for p in positions if p > n]
    if beyond:
        shown = ", ".join(map(str, beyond[:6])) + ("..." if len(beyond) > 6 else "")
        raise ValueError(
            f"{what}: chain {chain.name} has {n} polymer residues but the spec asks "
            f"for position(s) {shown} (max {max(positions)}); refusing a smaller fit"
        )
    xyz, seq, no_ca, offsets = [], [], [], set()
    for p in positions:
        res = chain.residues[p - 1]
        ca = res.get_ca()
        if ca is None:
            no_ca.append(f"{p} ({res.name} {res.seqid.num}{res.seqid.icode.strip()})")
            continue
        xyz.append([ca.pos.x, ca.pos.y, ca.pos.z])
        seq.append(one_letter(res))
        offsets.add(res.seqid.num - p)
    if no_ca:
        shown = ", ".join(no_ca[:6]) + ("..." if len(no_ca) > 6 else "")
        raise ValueError(
            f"{what}: chain {chain.name} position(s) {shown} have no CA atom"
        )
    return Picked(
        np.array(xyz, dtype=float), seq, offsets.pop() if len(offsets) == 1 else None
    )


# ---------------------------------------------------------------------------
# Geometry (same convention as motifrmsd)
# ---------------------------------------------------------------------------


def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Optimal rigid superposition of ``p`` onto ``q`` (both (N, 3), paired).

    Returns ``(rot, trans)`` with ``p @ rot.T + trans`` the moved ``p``. A
    reflection is excluded (the last singular vector is flipped when the naive
    rotation has determinant -1), so a mirror-image structure is never fitted.
    """
    if p.shape != q.shape or p.shape[0] < MIN_POSITIONS:
        raise ValueError(f"need >={MIN_POSITIONS} paired points, got {p.shape[0]}")
    p_mean, q_mean = p.mean(axis=0), q.mean(axis=0)
    pc, qc = p - p_mean, q - q_mean
    u, _, vt = np.linalg.svd(pc.T @ qc)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    trans = q_mean - rot @ p_mean
    return rot, trans


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD of two paired (N, 3) coordinate sets, as they stand (no fitting)."""
    return float(np.sqrt((((a - b) ** 2).sum(axis=1)).mean()))


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------


def parse_segments_json(text: str) -> tuple[list[str], list[list[int]]]:
    """The builder's already-parsed spec. Validated again here only for shape --
    a malformed payload is a builder bug and must not become a smaller fit."""
    try:
        data = json.loads(text)
        labels = [str(x) for x in data["labels"]]
        segments = [[int(p) for p in seg] for seg in data["segments"]]
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as e:
        raise ValueError(f"malformed --segments payload {text!r}: {e}") from e
    if len(labels) != len(segments) or not segments or any(not s for s in segments):
        raise ValueError(f"malformed --segments payload {text!r}")
    flat = [p for s in segments for p in s]
    if min(flat) < 1:
        raise ValueError(f"positions must be >= 1, got {min(flat)}")
    if len(set(flat)) != len(flat):
        raise ValueError("a position is listed more than once in the spec")
    if len(flat) < MIN_POSITIONS:
        raise ValueError(
            f"only {len(flat)} position(s); a superposition needs {MIN_POSITIONS}"
        )
    return labels, segments


@dataclass
class Result:
    metrics: dict[str, object]
    seg_rmsds: list[float]
    status: str
    moved: gemmi.Structure


def compute(args: "Args") -> Result:
    """Every posrmsd metric for one design, plus A moved into B's frame. Raises on
    anything unusable (the caller records the message as data)."""
    labels, segments = parse_segments_json(args.segments)
    positions = [p for s in segments for p in s]

    what_a = f"structure A (--col-a, {args.a})"
    what_b = f"structure B (--col-b, {args.b})"
    st_a = read_structure(args.a, "structure A (--col-a)")
    st_b = read_structure(args.b, "structure B (--col-b)")
    chain_a = select_chain(st_a, args.chain_a, what_a, "--chain-a")
    chain_b = select_chain(st_b, args.chain_b, what_b, "--chain-b")
    pick_a = pick_positions(chain_a, positions, what_a)
    pick_b = pick_positions(chain_b, positions, what_b)

    # One global fit of A's selected CA onto B's selected CA (B = reference frame).
    rot, trans = kabsch(pick_a.xyz, pick_b.xyz)
    moved = pick_a.xyz @ rot.T + trans
    deviations = np.sqrt(((moved - pick_b.xyz) ** 2).sum(axis=1))

    # Per-segment RMSDs, from that SAME fit (a loop that is locally right but
    # globally misplaced must show up; a per-segment refit would hide it).
    seg_rmsds, a_seqs, b_seqs = [], [], []
    start = 0
    for seg in segments:
        idx = slice(start, start + len(seg))
        seg_rmsds.append(rmsd(moved[idx], pick_b.xyz[idx]))
        a_seqs.append("".join(pick_a.seq[idx]))
        b_seqs.append("".join(pick_b.seq[idx]))
        start += len(seg)

    # 'X' never confirms a match: two unknowns agreeing says nothing about the
    # indexing, and a poly-UNK reference must not read as 1.0.
    matches = sum(a == b and a != "X" for a, b in zip(pick_a.seq, pick_b.seq))

    a_len, b_len = len(chain_a.residues), len(chain_b.residues)
    warnings = []
    if a_len != b_len:
        warnings.append(
            f"chain lengths differ (A chain {chain_a.name}: {a_len}, B chain "
            f"{chain_b.name}: {b_len}); the same-ordering premise is not guaranteed"
        )
    status = "OK" if not warnings else "warn: " + "; ".join(warnings)

    metrics: dict[str, object] = {
        "a_chain": chain_a.name,
        "b_chain": chain_b.name,
        "a_len": a_len,
        "b_len": b_len,
        "n_positions": len(positions),
        "seq_match_frac": matches / len(positions),
        "a_motif_seq": ",".join(a_seqs),
        "b_motif_seq": ",".join(b_seqs),
        "a_resnum_offset": pick_a.offset,
        "b_resnum_offset": pick_b.offset,
        "a_gapless": int(gapless(chain_a.residues)),
        "b_gapless": int(gapless(chain_b.residues)),
        "rmsd": rmsd(moved, pick_b.xyz),
        "max_dev": float(deviations.max()),
        "n_segments": len(segments),
        "worst_seg_rmsd": max(seg_rmsds),
        "seg_ranges": "/".join(labels),
        "seg_rmsds": "/".join(f"{v:.3f}" for v in seg_rmsds),
    }

    # Move the WHOLE of A into B's frame with the same rot/trans, for eyeballing.
    for model in st_a:
        for chain in model:
            for res in chain:
                for atom in res:
                    v = np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
                    x, y, z = rot @ v + trans
                    atom.pos = gemmi.Position(float(x), float(y), float(z))

    return Result(metrics, seg_rmsds, status, st_a)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def seg_columns(seg_rmsds: list[float]) -> dict[str, float]:
    """``seg0_rmsd``, ``seg1_rmsd``, ... -- one filterable column per segment."""
    return {f"seg{i}_rmsd": v for i, v in enumerate(seg_rmsds)}


def _fmt(value: object) -> str:
    """TSV cell: NA (None) becomes empty, floats get 3 decimals."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value).replace("\t", " ").replace("\n", " ")


def write_result(
    result_tsv: Path,
    name: str,
    status: str,
    sup_path: str,
    metrics: dict,
    extra_columns: dict[str, float] | None = None,
) -> None:
    """One-row TSV. Error rows carry every metric empty (NA), never 0."""
    extra_columns = extra_columns or {}
    columns = [*RESULT_COLUMNS, *extra_columns]
    result_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(result_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)
        writer.writerow(
            [name, _fmt(status), sup_path]
            + [_fmt(metrics.get(c)) for c in METRIC_COLUMNS]
            + [_fmt(v) for v in extra_columns.values()]
        )


class Args(argparse.Namespace):
    name: str
    a: Path
    b: Path
    chain_a: str
    chain_b: str
    segments: str
    out_struct: Path
    result_tsv: Path


def run(args: Args) -> None:
    """compute + write, with every failure recorded as an ``error:`` row."""
    try:
        result = compute(args)
        args.out_struct.parent.mkdir(parents=True, exist_ok=True)
        result.moved.setup_entities()
        result.moved.make_mmcif_document().write_file(str(args.out_struct))
        m = result.metrics
        print(
            f"{args.name}: rmsd {m['rmsd']:.3f} A over {m['n_positions']} CA "
            f"(A chain {m['a_chain']} vs B chain {m['b_chain']}), max_dev "
            f"{m['max_dev']:.3f}, segments {m['seg_ranges']} = {m['seg_rmsds']}, "
            f"seq_match_frac {m['seq_match_frac']:.3f}, status {result.status}"
        )
        write_result(
            args.result_tsv,
            args.name,
            result.status,
            str(args.out_struct),
            result.metrics,
            seg_columns(result.seg_rmsds),
        )
    except Exception as e:  # noqa: BLE001 - errors are recorded as data
        print(f"{args.name}: ERROR {e}", file=sys.stderr)
        write_result(args.result_tsv, args.name, f"error: {e}", "", {})


def main() -> None:
    ap = argparse.ArgumentParser(
        description="CA RMSD of structure A onto structure B over explicit "
        "1-indexed file-order chain positions."
    )
    ap.add_argument("--name", required=True)
    ap.add_argument("--a", type=Path, required=True, help="Structure A (moved).")
    ap.add_argument("--b", type=Path, required=True, help="Structure B (reference).")
    ap.add_argument("--chain-a", default="", help="Chain of A; empty = the only one.")
    ap.add_argument("--chain-b", default="", help="Chain of B; empty = the only one.")
    ap.add_argument(
        "--segments",
        required=True,
        help='Already-parsed spec as JSON: {"labels": [...], "segments": [[...]]}.',
    )
    ap.add_argument("--out-struct", type=Path, required=True)
    ap.add_argument("--result-tsv", type=Path, required=True)
    run(ap.parse_args(namespace=Args()))


if __name__ == "__main__":
    main()
