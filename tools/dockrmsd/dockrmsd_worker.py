#!/usr/bin/env python3
"""
DOCKING RMSD: did the model put the binder where the reference put it, allowing
for the target's symmetry?

THE MEASUREMENT, in four steps and no more:

1. Superpose the MODEL's TARGET chains onto the REFERENCE's TARGET chains
   (Kabsch, reflections excluded). **Target atoms only** -- the binder contributes
   nothing to the transform. That is what makes this a DOCKING metric and not a
   fold metric.
2. Apply that one transform T to the model's BINDER chain. ONCE.
3. RMSD between the transformed model binder and the reference binder, with NO
   further fitting.
4. The target is usually a homo-oligomer, so repeat 1-3 for EVERY permutation of
   the target chains and report the MINIMUM, naming the permutation that won.

PREMISE -- read before using any number.

* Both files hold the SAME target (same sequence) plus exactly one binder chain.
  The roles are taken from the flags and are NEVER inferred from chain length:
  a binder longer than a target protomer does not invert the measurement.
* THE TARGET'S PROTOMERS ARE PHYSICALLY INTERCHANGEABLE -- a true homo-oligomer.
  Minimising over chain permutations is only meaningful if every mapping names an
  EQUIVALENT binding site. On a PSEUDO-symmetric target (homologous but distinct
  protomers, similar enough to clear --min-target-identity) the minimum reports
  the best of several INEQUIVALENT sites, so a binder docked at site B of a
  reference designed for site A scores as though it hit the designed site. There
  is no way for this tool to detect that; pass --pin-mapping to evaluate only the
  chain order as given, which is the honest question on such a target.
* TARGET correspondence is PER ROLE and chosen with ``--target-match``. The
  binder's correspondence (file-order index, below) is NOT negotiable and is
  never sequence-aligned; applying ONE mode to both roles is precisely what
  disqualified DockQ for this comparison, and doing it in mirror image would be
  the same mistake.

  ``number``   -- intersect AUTHOR residue numbers (gemmi's ``seqid``, i.e.
                  auth_seq_id for mmCIF) per mapped chain pair. Right when the
                  two files share a numbering convention.
  ``sequence`` -- gapped global alignment of the two target chains. Right when
                  one file MODELS RESIDUES THE OTHER LACKS. Measured on this
                  campaign's pilot pair: the Boltz model numbers its target
                  chains 1-149 and the rfd3 reference numbers its 1-139, BOTH
                  FROM 1, and the model carries an extra 10-residue loop from
                  position 95. Under ``number`` the author intersection is 139
                  residues of which only 97 are identical (0.698) because the
                  last 45 are shifted by 10 -- so every permutation is skipped
                  and the row errors. No renumbering fixes that; only an
                  alignment does.

  There is no fallback between the modes: ``target_match`` records which one
  produced the pairing on every row.

  ``n_target_res_used`` is how many residue pairs the fit used (aligned pairs
  under ``sequence``, shared numbers under ``number``);
  ``n_target_res_used_per_chain`` shows how they split, and the floor is applied
  PER PAIR so one chain cannot carry the fit alone. ``target_seq_match_frac`` is
  the identity over those pairs -- the alignment's identity under ``sequence``.
  ``target_resnum_range_model`` / ``target_resnum_range_ref`` record the two
  numbering schemes outright.

  NOTE what the identity fraction can and cannot catch UNDER ``number``: it is
  ``matches / n_used`` over the intersection, so mis-registered pairs sit in the
  denominator and a numbering break at position p of N scores ``(p-1)/N``. A
  break near the middle is caught; a break near the C-terminus (say 120 of 139 ->
  0.86) clears a 0.8 floor and fits the last ~19 residues against the wrong
  partners. It degrades, it does not invert -- but read the two range columns
  rather than trusting the floor to notice, and prefer ``sequence`` whenever the
  two files can differ in modelled extent.
* BINDER correspondence is by FILE-ORDER INDEX: position N of the model binder is
  position N of the reference binder. Nothing is aligned. The two binders are
  expected to be the SAME BACKBONE carrying DIFFERENT SEQUENCES (an rfd3 backbone
  vs the sequence designed onto it), so sequence identity is NOT the relationship
  and is never required -- 0.45-0.65 is normal. The lengths MUST be equal; they
  are the only thing that makes index correspondence defined, so an inequality is
  an ``error:`` row, never a trimmed comparison.

WHEN THE NUMBERS MEAN NOTHING

* ``target_fit_rmsd`` is large: the model did not reproduce the target, so the
  frame the binder was moved into is not the reference frame, and the binder RMSD
  is meaningless at ANY magnitude. Above ``--max-target-fit-rmsd`` the status is
  ``warn:`` (numbers kept, excluded by ``== "OK"``).
* ``n_target_res_used`` is small relative to the chains: the two files disagree on
  author numbering, so the "shared" residues are few and possibly mispaired.
* The two files do not hold the same target, or the chain flags name the wrong
  molecules: ``target_seq_match_frac`` collapses and the row errors.
* The binder is not rigid between the two files (a different fold, a different
  construct): this tool measures PLACEMENT and will report the fold error folded
  into the placement error. Use ``usalign``/``posrmsd`` for the fold.

It does NOT answer: whether the binder folded (usalign), whether it touches the
intended epitope (epitope/ifacegeom), or anything about side chains, chemistry or
interface quality (cms, pyrosetta).

OUTPUT: one TSV per design (name, status, sup_path, metrics...) plus the whole
model moved into the reference frame under the winning transform, as mmCIF.

Pure gemmi + numpy; no prosapia import (the task image does not carry it).

Usage (standalone):
    python dockrmsd_worker.py --task-file task.tsv \\
        --binder-chain B --target-chains A,C,D \\
        --ref-binder-chain B --ref-target-chains A,C,D \\
        --results-dir out/
"""

import argparse
import csv
import itertools
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np

# Kabsch needs three non-collinear points.
MIN_FIT_POINTS = 3

PEPTIDE_TYPES = {gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD}

# Atom selections. 'ca' is the default and the only one validated on real data.
FIT_ATOM_SETS = {"ca": ("CA",), "backbone": ("N", "CA", "C", "O")}

# Metric columns of the per-design TSV, in order. Bare names; the collector adds
# the comparison prefix and the driver adds the tool leaf.
TRUST_COLUMNS = [
    "binder_chain",
    "ref_binder_chain",
    "target_chains",
    "ref_target_chains",
    "fit_atoms",
    "target_match",
    "pin_mapping",
    "model_chain_ids",
    "ref_chain_ids",
    "binder_n_res",
    "ref_binder_n_res",
    "binder_seq_match_frac",
    "target_fit_rmsd",
    "best_target_fit_rmsd",
    "n_target_res_used",
    "n_target_res_used_per_chain",
    "n_target_res_model",
    "n_target_res_ref",
    "target_res_used_frac",
    "target_seq_match_frac",
    "target_resnum_range_model",
    "target_resnum_range_ref",
    "n_fit_atoms",
    "n_binder_atoms",
    "n_permutations",
    "n_candidates",
    "permutation",
    "margin",
    "rmsd_all",
]
SCIENCE_COLUMNS = ["rmsd", "max_dev", "complex_rmsd"]
METRIC_COLUMNS = TRUST_COLUMNS + SCIENCE_COLUMNS
RESULT_COLUMNS = ["name", "status", "sup_path", *METRIC_COLUMNS]


class MeasureError(Exception):
    """A refusal that is recorded as an ``error:`` row. ``metrics`` carries the
    trust values already established (chain ids, binder lengths, ...) so the row
    still shows WHY it was refused; every science column stays NA."""

    def __init__(self, message: str, metrics: dict | None = None):
        super().__init__(message)
        self.metrics = metrics or {}


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
        raise MeasureError(f"{what}: no path given")
    if not path.is_file():
        raise MeasureError(f"{what} not found: {path}")
    try:
        st = gemmi.read_structure(str(path))
    except Exception as e:  # noqa: BLE001 - reported as data by the caller
        raise MeasureError(f"{what} {path} could not be parsed: {e}") from e
    if len(st) == 0:
        raise MeasureError(f"{what} {path} holds no model")
    if len(st) > 1:
        # Taking st[0] silently would measure one member of an ensemble and
        # report it as "the" structure. rpxdock's input guard is the precedent.
        raise MeasureError(
            f"{what} {path} holds {len(st)} models; dockrmsd measures ONE "
            f"structure and will not pick a member of an ensemble for you. Split "
            f"the file first"
        )
    st.remove_alternative_conformations()  # deterministic atom pick
    st.merge_chain_parts()
    st.setup_entities()
    return st


def polymer_chains(st: gemmi.Structure) -> dict[str, list[gemmi.Residue]]:
    """Every protein chain of the first model, name -> polymer residues in FILE
    ORDER. A chain name that appears twice after merging is a refusal: silently
    taking the first is how a comparison ends up on the wrong molecule."""
    model = st[0]
    names = [c.name for c in model]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        raise MeasureError(
            f"more than one chain named {', '.join(dupes)} after merging chain "
            f"parts; cannot tell which one is meant"
        )
    out: dict[str, list[gemmi.Residue]] = {}
    for chain in model:
        polymer = chain.get_polymer()
        if len(polymer) == 0:
            continue
        if polymer.check_polymer_type() not in PEPTIDE_TYPES:
            continue
        out[chain.name] = list(polymer)
    return out


def select(
    chains: dict[str, list[gemmi.Residue]], wanted: list[str], what: str, flag: str
) -> dict[str, list[gemmi.Residue]]:
    """The named protein chains, or a refusal listing what the file actually has."""
    missing = [c for c in wanted if c not in chains]
    if missing:
        raise MeasureError(
            f"{what}: chain(s) {', '.join(missing)} ({flag}) are not protein chains "
            f"of the file (protein chains present: "
            f"{', '.join(chains) or 'none'})"
        )
    return {c: chains[c] for c in wanted}


# Needleman-Wunsch with affine (Gotoh) gaps, identity scoring. Enough for "the
# same protein, one of the two modelling residues the other does not"; it is not
# a substitution-matrix aligner and is not meant to align homologues.
NW_MATCH, NW_MISMATCH, NW_GAP_OPEN, NW_GAP_EXTEND = 2, -1, -11, -1


def align_sequences(a: str, b: str) -> list[tuple[int, int]]:
    """Global alignment of two one-letter sequences; returns the 0-indexed
    (i_a, i_b) pairs of the aligned columns, gaps excluded, in order.

    Substituted columns are KEPT -- they are aligned positions, and dropping them
    (as DockQ's sequence mode does) would silently shrink the comparison to the
    coincidentally identical part.
    """
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return []
    neg = float("-inf")
    # M: a_i with b_j; X: a_i with a gap; Y: b_j with a gap.
    mat = [[neg] * (m + 1) for _ in range(n + 1)]
    xg = [[neg] * (m + 1) for _ in range(n + 1)]
    yg = [[neg] * (m + 1) for _ in range(n + 1)]
    ptr = [[0] * (m + 1) for _ in range(n + 1)]  # 1=M 2=X 3=Y, per layer below
    px = [[0] * (m + 1) for _ in range(n + 1)]
    py = [[0] * (m + 1) for _ in range(n + 1)]
    mat[0][0] = 0.0
    for i in range(1, n + 1):
        xg[i][0] = NW_GAP_OPEN + i * NW_GAP_EXTEND
        px[i][0] = 2
    for j in range(1, m + 1):
        yg[0][j] = NW_GAP_OPEN + j * NW_GAP_EXTEND
        py[0][j] = 3
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            score = NW_MATCH if ai == b[j - 1] else NW_MISMATCH
            best, which = mat[i - 1][j - 1], 1
            if xg[i - 1][j - 1] > best:
                best, which = xg[i - 1][j - 1], 2
            if yg[i - 1][j - 1] > best:
                best, which = yg[i - 1][j - 1], 3
            mat[i][j], ptr[i][j] = best + score, which

            open_x = mat[i - 1][j] + NW_GAP_OPEN + NW_GAP_EXTEND
            ext_x = xg[i - 1][j] + NW_GAP_EXTEND
            xg[i][j], px[i][j] = (open_x, 1) if open_x >= ext_x else (ext_x, 2)

            open_y = mat[i][j - 1] + NW_GAP_OPEN + NW_GAP_EXTEND
            ext_y = yg[i][j - 1] + NW_GAP_EXTEND
            yg[i][j], py[i][j] = (open_y, 1) if open_y >= ext_y else (ext_y, 3)

    layer = max((mat[n][m], 1), (xg[n][m], 2), (yg[n][m], 3))[1]
    pairs: list[tuple[int, int]] = []
    i, j = n, m
    while i > 0 or j > 0:
        if layer == 1:
            pairs.append((i - 1, j - 1))
            layer = ptr[i][j]
            i, j = i - 1, j - 1
        elif layer == 2:
            layer = px[i][j]
            i -= 1
        else:
            layer = py[i][j]
            j -= 1
        if i == 0 and j > 0:
            layer = 3
        elif j == 0 and i > 0:
            layer = 2
    pairs.reverse()
    return pairs


def correspondence(
    m_residues: list[gemmi.Residue],
    r_residues: list[gemmi.Residue],
    mode: str,
    m_chain: str,
    r_chain: str,
) -> list[tuple[gemmi.Residue, gemmi.Residue]]:
    """Which model target residue answers to which reference target residue.

    ``number``  -- intersect AUTHOR residue numbers (+ insertion code). Correct
                   when the two files share a numbering convention; it reports
                   nothing at all when they do not.
    ``sequence`` -- gapped global alignment of the two chains' sequences. Correct
                   when one file MODELS RESIDUES THE OTHER LACKS (an unresolved
                   loop), which no renumbering can fix because both files may
                   already be numbered 1..N.

    There is NO fallback between the two: a caller who asked for ``number`` and
    silently got an alignment would have no way to know which question was
    answered. ``target_match`` records the mode on every row.
    """
    if mode == "number":
        by_key: dict[tuple[int, str], gemmi.Residue] = {}
        for res in m_residues:
            key = res_key(res)
            if key in by_key:
                # Keeping the first silently would pair the reference residue
                # with whichever copy happened to be written first.
                raise MeasureError(
                    f"model chain {m_chain} has more than one residue numbered "
                    f"{key[0]}{key[1]}; author numbering must identify a residue "
                    f"uniquely for --target-match number to mean anything"
                )
            by_key[key] = res
        seen: set[tuple[int, str]] = set()
        out = []
        for r_res in r_residues:
            key = res_key(r_res)
            if key in seen:
                raise MeasureError(
                    f"reference chain {r_chain} has more than one residue numbered "
                    f"{key[0]}{key[1]}; author numbering must identify a residue "
                    f"uniquely for --target-match number to mean anything"
                )
            seen.add(key)
            if (m_res := by_key.get(key)) is not None:
                out.append((m_res, r_res))
        return out

    if mode == "sequence":
        m_seq = "".join(one_letter(r) for r in m_residues)
        r_seq = "".join(one_letter(r) for r in r_residues)
        return [
            (m_residues[i], r_residues[j]) for i, j in align_sequences(m_seq, r_seq)
        ]

    raise MeasureError(f"unknown --target-match mode {mode!r}")


def resnum_ranges(chains: dict[str, list[gemmi.Residue]], order: list[str]) -> str:
    """``'9-157;9-157;9-157'`` -- the author numbering each named chain spans,
    in the order given. An OBSERVATION of the numbering scheme, not an inference
    from how much of it matched."""
    out = []
    for name in order:
        residues = chains[name]
        out.append(f"{residues[0].seqid.num}-{residues[-1].seqid.num}")
    return ";".join(out)


def res_key(res: gemmi.Residue) -> tuple[int, str]:
    """Author residue number + insertion code -- the identity a target residue is
    intersected on."""
    return (res.seqid.num, res.seqid.icode.strip())


def atoms_of(res: gemmi.Residue, names: tuple[str, ...]) -> list[list[float]] | None:
    """The requested atoms of one residue, or None when any of them is absent."""
    picked = []
    for atom_name in names:
        atom = res.find_atom(atom_name, "*")
        if atom is None:
            return None
        picked.append([atom.pos.x, atom.pos.y, atom.pos.z])
    return picked


# ---------------------------------------------------------------------------
# Geometry (same convention as posrmsd / motifrmsd)
# ---------------------------------------------------------------------------


def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Optimal rigid superposition of ``p`` onto ``q`` (both (N, 3), paired).

    Returns ``(rot, trans)`` with ``p @ rot.T + trans`` the moved ``p``. A
    reflection is excluded, so a mirror image is never fitted.
    """
    if p.shape != q.shape or p.shape[0] < MIN_FIT_POINTS:
        raise MeasureError(f"need >={MIN_FIT_POINTS} paired points, got {p.shape[0]}")
    p_mean, q_mean = p.mean(axis=0), q.mean(axis=0)
    pc, qc = p - p_mean, q - q_mean
    u, _, vt = np.linalg.svd(pc.T @ qc)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    trans = q_mean - rot @ p_mean
    return rot, trans


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD of two paired (N, 3) coordinate sets AS THEY STAND (no fitting)."""
    return float(np.sqrt((((a - b) ** 2).sum(axis=1)).mean()))


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------


@dataclass
class Binder:
    """One binder chain, atom-selected by file-order index."""

    chain: str
    n_res: int
    xyz: np.ndarray
    seq: str


@dataclass
class Permutation:
    """One assignment of model target chains to reference target chains."""

    mapping: str
    pairs: list[tuple[str, str]]
    skipped: str = ""  # non-empty = excluded from the minimum, with the reason
    candidate: bool = False  # survived the fit tolerance; marked '*' in the audit
    per_pair: list[int] = field(default_factory=list)
    identity: float = float("nan")
    fit_rmsd: float = float("nan")
    binder_rmsd: float = float("nan")
    rot: np.ndarray = field(default_factory=lambda: np.eye(3))
    trans: np.ndarray = field(default_factory=lambda: np.zeros(3))
    model_xyz: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    ref_xyz: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    @property
    def n_res(self) -> int:
        """Residues entering the fit, summed over the chain pairs. Only ever read
        AFTER the per-pair floor has been applied."""
        return sum(self.per_pair)

    def audit(self) -> str:
        """One permutation's entry in ``rmsd_all``.

        A leading ``*`` marks a CANDIDATE -- a permutation that competed for the
        minimum -- so the candidate set is read off the audit rather than
        re-derived from the fit values by hand.
        """
        per_chain = "/".join(str(c) for c in self.per_pair)
        if self.skipped:
            return f"{self.mapping}:{self.skipped}:per_chain={per_chain}"
        star = "*" if self.candidate else ""
        return (
            f"{star}{self.mapping}:fit={self.fit_rmsd:.3f}"
            f":binder={self.binder_rmsd:.3f}:n={self.n_res}"
            f":per_chain={per_chain}:id={self.identity:.3f}"
        )


def take_binder(
    chains: dict[str, list[gemmi.Residue]],
    chain: str,
    atom_names: tuple[str, ...],
    what: str,
    flag: str,
) -> Binder:
    """Every atom of the binder chain, in file order. A residue missing any
    requested atom is a refusal: index correspondence needs all of them."""
    residues = select(chains, [chain], what, flag)[chain]
    xyz, seq, incomplete = [], [], []
    for i, res in enumerate(residues, start=1):
        picked = atoms_of(res, atom_names)
        if picked is None:
            incomplete.append(f"{i} ({res.name} {res.seqid.num})")
            continue
        xyz.extend(picked)
        seq.append(one_letter(res))
    if incomplete:
        shown = ", ".join(incomplete[:6]) + ("..." if len(incomplete) > 6 else "")
        raise MeasureError(
            f"{what}: binder chain {chain} position(s) {shown} lack one of the "
            f"{'/'.join(atom_names)} atoms; refusing a comparison over a smaller set"
        )
    return Binder(chain, len(residues), np.array(xyz, dtype=float), "".join(seq))


def pair_target(
    corr: dict[tuple[str, str], list[tuple[gemmi.Residue, gemmi.Residue]]],
    pairs: list[tuple[str, str]],
    atom_names: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray, list[int], float]:
    """Paired target atoms for one chain mapping, from the already-computed
    per-chain-pair correspondence (cached, so an alignment is never redone for a
    chain pair that reappears in another permutation).

    A residue missing a requested atom on either side is left out.

    Returns (model_xyz, ref_xyz, per_pair_counts, identity_fraction). The counts
    are PER CHAIN PAIR and never summed here: a pooled total hides the case where
    one pair carries the whole correspondence and the others contribute nothing,
    which turns the "whole-target" fit into a single-protomer fit (see compute()).
    """
    m_list, r_list, per_pair, matches, n_used = [], [], [], 0, 0
    for m_chain, r_chain in pairs:
        pair_count = 0
        for m_res, r_res in corr[(m_chain, r_chain)]:
            m_atoms = atoms_of(m_res, atom_names)
            r_atoms = atoms_of(r_res, atom_names)
            if m_atoms is None or r_atoms is None:
                continue
            m_list.extend(m_atoms)
            r_list.extend(r_atoms)
            pair_count += 1
            n_used += 1
            m_code, r_code = one_letter(m_res), one_letter(r_res)
            # 'X' never confirms a match: two unknowns agreeing says nothing.
            matches += int(m_code == r_code and m_code != "X")
        per_pair.append(pair_count)
    identity = matches / n_used if n_used else float("nan")
    return (
        np.array(m_list, dtype=float).reshape(-1, 3),
        np.array(r_list, dtype=float).reshape(-1, 3),
        per_pair,
        identity,
    )


@dataclass
class Result:
    metrics: dict[str, object]
    status: str
    moved: gemmi.Structure


def compute(
    model_path: Path,
    ref_path: Path,
    opt: "Options",
) -> Result:
    """Every dockrmsd metric for one design, plus the model moved into the
    reference frame. Raises MeasureError (with partial trust metrics) on anything
    unusable; the caller records the message as data."""
    metrics: dict[str, object] = opt.echo()
    atom_names = FIT_ATOM_SETS[opt.fit_atoms]

    def fail(msg: str):
        return MeasureError(msg, dict(metrics))

    try:
        st_model = read_structure(model_path, "model (--col-a)")
        st_ref = read_structure(ref_path, "reference (--col-b)")
        model_chains = polymer_chains(st_model)
        ref_chains = polymer_chains(st_ref)
    except MeasureError as e:
        raise MeasureError(str(e), dict(metrics)) from None
    metrics["model_chain_ids"] = ",".join(model_chains) or "none"
    metrics["ref_chain_ids"] = ",".join(ref_chains) or "none"

    try:
        model_binder = take_binder(
            model_chains, opt.binder_chain, atom_names, "model (--col-a)", "--binder-chain"
        )
        ref_binder = take_binder(
            ref_chains,
            opt.ref_binder_chain,
            atom_names,
            "reference (--col-b)",
            "--ref-binder-chain",
        )
    except MeasureError as e:
        raise MeasureError(str(e), dict(metrics)) from None
    metrics["binder_n_res"] = model_binder.n_res
    metrics["ref_binder_n_res"] = ref_binder.n_res

    # THE equal-length check. Index correspondence is undefined otherwise, and a
    # trimmed comparison would be a plausible number for a question nobody asked.
    if model_binder.n_res != ref_binder.n_res:
        raise fail(
            f"binder length mismatch: model chain {model_binder.chain} has "
            f"{model_binder.n_res} residues, reference chain {ref_binder.chain} has "
            f"{ref_binder.n_res}. dockrmsd pairs binder residues by FILE-ORDER INDEX, "
            f"which is undefined unless the lengths are equal -- it never trims. Check "
            f"the chain flags, or use chainsel to cut matching chains"
        )
    if model_binder.xyz.shape[0] == 0:
        raise fail(f"binder chain {model_binder.chain} has no residues")
    metrics["n_binder_atoms"] = int(model_binder.xyz.shape[0])
    binder_matches = sum(
        a == b and a != "X" for a, b in zip(model_binder.seq, ref_binder.seq)
    )
    metrics["binder_seq_match_frac"] = binder_matches / len(model_binder.seq)

    try:
        model_target = select(
            model_chains, opt.target_chains, "model (--col-a)", "--target-chains"
        )
        ref_target = select(
            ref_chains, opt.ref_target_chains, "reference (--col-b)", "--ref-target-chains"
        )
    except MeasureError as e:
        raise MeasureError(str(e), dict(metrics)) from None
    metrics["n_target_res_model"] = sum(len(v) for v in model_target.values())
    metrics["n_target_res_ref"] = sum(len(v) for v in ref_target.values())
    # Numbering PROVENANCE, not an inference from it: the author range each
    # target chain actually spans, so "the two files use different numbering
    # schemes" is an observation in the table rather than something a reader has
    # to deduce from a shrunken intersection. '9-157;9-157;9-157' against
    # '1-139;1-139;1-139' is the rfd3-renumbering signature.
    metrics["target_resnum_range_model"] = resnum_ranges(model_target, opt.target_chains)
    metrics["target_resnum_range_ref"] = resnum_ranges(ref_target, opt.ref_target_chains)

    # ---- the chain mappings to try ----
    # --pin-mapping evaluates ONLY the order as given: the honest question on a
    # PSEUDO-symmetric target, whose protomers are similar enough to pass the
    # identity floor but are NOT physically interchangeable, so a minimum over
    # mappings would report the best of several INEQUIVALENT sites.
    orders = (
        [tuple(opt.target_chains)]
        if opt.pin_mapping
        else list(itertools.permutations(opt.target_chains))
    )
    # The correspondence depends only on the CHAIN PAIR, never on the permutation
    # it appears in, so it is computed once per pair and reused -- 9 alignments
    # for a trimer rather than 18.
    corr: dict[tuple[str, str], list[tuple[gemmi.Residue, gemmi.Residue]]] = {}
    try:
        for m_chain in opt.target_chains:
            for r_chain in opt.ref_target_chains:
                corr[(m_chain, r_chain)] = correspondence(
                    model_target[m_chain], ref_target[r_chain],
                    opt.target_match, m_chain, r_chain,
                )
    except MeasureError as e:
        raise MeasureError(str(e), dict(metrics)) from None

    perms: list[Permutation] = []
    for order in orders:
        pairs = list(zip(order, opt.ref_target_chains))
        perm = Permutation(
            mapping=",".join(f"{m}>{r}" for m, r in pairs), pairs=pairs
        )
        m_xyz, r_xyz, per_pair, identity = pair_target(corr, pairs, atom_names)
        perm.per_pair, perm.identity = per_pair, identity
        n_used = sum(per_pair)
        # EVERY chain pair must clear the floor, not their pooled sum. Pooling
        # admits a mapping in which one pair carries the whole intersection and
        # the others contribute nothing -- a SINGLE-PROTOMER superposition
        # reported under a column name that promises the whole target frame, and
        # one that also defeats --fit-tolerance, because with one protomer in the
        # fit every permutation fits near-equally well.
        if min(per_pair) < opt.min_target_res:
            perm.skipped = (
                f"skipped_few_shared_residues(min_pair={min(per_pair)}"
                f"<{opt.min_target_res})"
            )
        elif m_xyz.shape[0] < MIN_FIT_POINTS:
            perm.skipped = f"skipped_few_shared_residues(n={n_used})"
        elif identity < opt.min_target_identity:
            perm.skipped = f"skipped_identity({identity:.3f})"
        else:
            rot, trans = kabsch(m_xyz, r_xyz)
            moved_target = m_xyz @ rot.T + trans
            perm.rot, perm.trans = rot, trans
            perm.fit_rmsd = rmsd(moved_target, r_xyz)
            moved_binder = model_binder.xyz @ rot.T + trans
            perm.binder_rmsd = rmsd(moved_binder, ref_binder.xyz)
            perm.model_xyz, perm.ref_xyz = m_xyz, r_xyz
        perms.append(perm)

    metrics["n_permutations"] = len(perms)
    metrics["rmsd_all"] = ";".join(p.audit() for p in perms)

    usable = [p for p in perms if not p.skipped]
    if not usable:
        raise fail(
            f"no chain permutation gave a usable target superposition: "
            f"{metrics['rmsd_all']}. Either the two files do not hold the same "
            f"target (--min-target-identity {opt.min_target_identity}), or some "
            f"chain pair shares too few author residue numbers (--min-target-res "
            f"{opt.min_target_res}, applied PER CHAIN PAIR -- see per_chain= in "
            f"the audit); check the chain flags and the numbering, and compare "
            f"target_resnum_range_model with target_resnum_range_ref"
        )

    best_fit = min(p.fit_rmsd for p in usable)
    metrics["best_target_fit_rmsd"] = best_fit
    # A permutation is a CANDIDATE only if its target fit is as good as the best
    # one found, within --fit-tolerance. On a Cn homo-oligomer that admits the n
    # genuine symmetry operations and rejects the n!-n impostor mappings, whose
    # fit is tens of angstroms worse -- so the minimum is taken over real
    # symmetry copies and not over an arbitrary rotation that happened to land
    # the binder near the reference.
    candidates = [p for p in usable if p.fit_rmsd <= best_fit + opt.fit_tolerance]
    for p in candidates:
        p.candidate = True
    metrics["n_candidates"] = len(candidates)
    # Rebuilt now that the candidates are known, so the audit carries the '*'
    # marks and a reader never re-derives the candidate set by hand.
    metrics["rmsd_all"] = ";".join(p.audit() for p in perms)

    ordered = sorted(candidates, key=lambda p: p.binder_rmsd)
    winner = ordered[0]
    metrics["permutation"] = winner.mapping
    metrics["target_fit_rmsd"] = winner.fit_rmsd
    metrics["n_target_res_used"] = winner.n_res
    metrics["n_target_res_used_per_chain"] = "/".join(str(c) for c in winner.per_pair)
    metrics["target_seq_match_frac"] = winner.identity
    metrics["n_fit_atoms"] = int(winner.model_xyz.shape[0])
    denom = min(int(metrics["n_target_res_model"]), int(metrics["n_target_res_ref"]))
    metrics["target_res_used_frac"] = winner.n_res / denom if denom else None
    # NA, never 0: with one candidate there is no second-best to be close to, and
    # a 0 here would read as "two placements tied" -- the opposite of the truth.
    metrics["margin"] = (
        ordered[1].binder_rmsd - winner.binder_rmsd if len(ordered) > 1 else None
    )

    moved_binder = model_binder.xyz @ winner.rot.T + winner.trans
    deviations = np.sqrt(((moved_binder - ref_binder.xyz) ** 2).sum(axis=1))
    metrics["rmsd"] = winner.binder_rmsd
    metrics["max_dev"] = float(deviations.max())
    # EXACTLY the atom-count-weighted quadrature mean of target_fit_rmsd and
    # rmsd -- it carries NO independent information, and it is NOT "dominated by
    # the target": the terms are quadratic, so once rmsd exceeds ~2 A the binder
    # term leads. It always sits BELOW rmsd (a compressed, monotone copy), so a
    # filter written against it keeps designs far worse than it appears to admit.
    # Kept only because that exactness makes it a free arithmetic self-check.
    moved_target = winner.model_xyz @ winner.rot.T + winner.trans
    metrics["complex_rmsd"] = rmsd(
        np.vstack([moved_target, moved_binder]),
        np.vstack([winner.ref_xyz, ref_binder.xyz]),
    )

    warnings = []
    if winner.fit_rmsd > opt.max_target_fit_rmsd:
        warnings.append(
            f"target_fit_rmsd {winner.fit_rmsd:.3f} A exceeds "
            f"--max-target-fit-rmsd {opt.max_target_fit_rmsd}; the model did not "
            f"reproduce the target, so the binder RMSD is not a docking number"
        )
    used_frac = metrics["target_res_used_frac"]
    if used_frac is not None and float(used_frac) < opt.min_target_res_frac:
        warnings.append(
            f"only {float(used_frac):.3f} of the shorter target entered the fit "
            f"(--min-target-res-frac {opt.min_target_res_frac}); the two files "
            f"agree on author numbering over part of the target only -- compare "
            f"target_resnum_range_model with target_resnum_range_ref"
        )
    status = "OK" if not warnings else "warn: " + "; ".join(warnings)

    # The WHOLE model moved into the reference frame, for eyeballing.
    for model in st_model:
        for chain in model:
            for res in chain:
                for atom in res:
                    v = np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
                    x, y, z = winner.rot @ v + winner.trans
                    atom.pos = gemmi.Position(float(x), float(y), float(z))
    return Result(metrics, status, st_model)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def _fmt(value: object) -> str:
    """TSV cell: NA (None/NaN) becomes empty, floats get 4 decimals."""
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return ""
        return f"{value:.4f}"
    return str(value).replace("\t", " ").replace("\n", " ")


def write_result(
    result_tsv: Path, name: str, status: str, sup_path: str, metrics: dict
) -> None:
    """One-row TSV. A non-OK row carries NA (never 0) in every science column."""
    result_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(result_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(RESULT_COLUMNS)
        writer.writerow(
            [name, _fmt(status), sup_path]
            + [_fmt(metrics.get(c)) for c in METRIC_COLUMNS]
        )


@dataclass
class Options:
    binder_chain: str
    target_chains: list[str]
    ref_binder_chain: str
    ref_target_chains: list[str]
    fit_atoms: str
    target_match: str
    min_target_res: int
    min_target_res_frac: float
    min_target_identity: float
    fit_tolerance: float
    max_target_fit_rmsd: float
    pin_mapping: bool = False

    def echo(self) -> dict[str, object]:
        """The flag echo, available BEFORE anything is parsed -- so an unexpected
        failure still records which chains were asked for."""
        return {
            "binder_chain": self.binder_chain,
            "ref_binder_chain": self.ref_binder_chain,
            "target_chains": ",".join(self.target_chains),
            "ref_target_chains": ",".join(self.ref_target_chains),
            "fit_atoms": self.fit_atoms,
            "target_match": self.target_match,
            "pin_mapping": int(self.pin_mapping),
        }


def one_design(
    name: str, model: str, ref: str, opt: Options, results_dir: Path
) -> None:
    """Measure one design and write its TSV, recording any refusal as data."""
    result_tsv = results_dir / f"{name}.tsv"
    sup = results_dir / f"{name}_sup.cif"
    try:
        result = compute(Path(model), Path(ref), opt)
    except MeasureError as e:
        print(f"{name}: ERROR {e}", file=sys.stderr)
        write_result(result_tsv, name, f"error: {e}", "", e.metrics)
        return
    except Exception as e:  # noqa: BLE001 - errors are recorded as data
        # An UNEXPECTED failure is exactly where a diagnosis is scarcest, so the
        # flag echo goes on the row too -- never an empty metrics dict.
        print(f"{name}: ERROR {e}", file=sys.stderr)
        write_result(result_tsv, name, f"error: {type(e).__name__}: {e}", "", opt.echo())
        return
    sup.parent.mkdir(parents=True, exist_ok=True)
    result.moved.setup_entities()
    result.moved.make_mmcif_document().write_file(str(sup))
    m = result.metrics
    print(
        f"{name}: docking rmsd {m['rmsd']:.3f} A ({m['permutation']}), target fit "
        f"{m['target_fit_rmsd']:.3f} A over {m['n_target_res_used']} residues, "
        f"margin {_fmt(m['margin']) or 'NA'}, {m['n_candidates']}/"
        f"{m['n_permutations']} candidate permutations, status {result.status}"
    )
    write_result(result_tsv, name, result.status, str(sup), result.metrics)


def _chains(spec: str) -> list[str]:
    return [c.strip() for c in spec.split(",") if c.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Docking RMSD: fit on the target, move the binder once, "
        "minimise over target-chain permutations."
    )
    ap.add_argument("--task-file", type=Path, required=True,
                    help="TSV: name<TAB>model<TAB>reference, one design per line.")
    ap.add_argument("--binder-chain", required=True, help="Binder chain in the MODEL.")
    ap.add_argument("--target-chains", required=True, help="Target chains in the MODEL (comma).")
    ap.add_argument("--ref-binder-chain", required=True, help="Binder chain in the REFERENCE.")
    ap.add_argument("--ref-target-chains", required=True,
                    help="Target chains in the REFERENCE (comma), same count.")
    ap.add_argument("--fit-atoms", choices=sorted(FIT_ATOM_SETS), default="ca")
    ap.add_argument("--target-match", choices=("number", "sequence"), default="number")
    ap.add_argument("--min-target-res", type=int, default=10,
                    help="Per CHAIN PAIR, not pooled.")
    ap.add_argument("--min-target-res-frac", type=float, default=0.75)
    ap.add_argument("--min-target-identity", type=float, default=0.8)
    ap.add_argument("--fit-tolerance", type=float, default=1.0)
    ap.add_argument("--max-target-fit-rmsd", type=float, default=5.0)
    ap.add_argument("--pin-mapping", default="0",
                    help="'1' = evaluate only the chain order as given.")
    ap.add_argument("--results-dir", type=Path, required=True)
    args = ap.parse_args()

    opt = Options(
        binder_chain=args.binder_chain.strip(),
        target_chains=_chains(args.target_chains),
        ref_binder_chain=args.ref_binder_chain.strip(),
        ref_target_chains=_chains(args.ref_target_chains),
        fit_atoms=args.fit_atoms,
        target_match=args.target_match,
        min_target_res=args.min_target_res,
        min_target_res_frac=args.min_target_res_frac,
        min_target_identity=args.min_target_identity,
        fit_tolerance=args.fit_tolerance,
        max_target_fit_rmsd=args.max_target_fit_rmsd,
        pin_mapping=str(args.pin_mapping).strip() not in ("", "0", "false", "False"),
    )
    if len(opt.target_chains) != len(opt.ref_target_chains):
        raise SystemExit(
            f"--target-chains ({len(opt.target_chains)}) and --ref-target-chains "
            f"({len(opt.ref_target_chains)}) must name the same number of chains."
        )
    args.results_dir.mkdir(parents=True, exist_ok=True)

    with open(args.task_file) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            name = fields[0]
            model = fields[1] if len(fields) > 1 else ""
            ref = fields[2] if len(fields) > 2 else ""
            one_design(name, model, ref, opt, args.results_dir)


if __name__ == "__main__":
    main()
