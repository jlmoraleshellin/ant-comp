#!/usr/bin/env python3
"""
Motif CA RMSD of one rfdiffusion3 motif-scaffolding design against its reference.

PREMISE -- read this before using any number below.

This measures **one thing**: did rfd3 rebuild the grafted motif in the reference's
3D geometry? It answers it by taking rfd3's **own** per-design record of which
reference residue became which output residue (``diffused_index_map`` in the
sidecar JSON written beside each ``.cif.gz``), lifting the CA atom of every mapped
pair out of both structures, superposing the design's motif CA onto the
reference's motif CA with a least-squares (SVD/Kabsch) fit, and reporting the RMSD
of those same atoms after that fit.

It is therefore **scoped to rfd3 motif-scaffolding output and nothing else**:

* No sidecar JSON, or a sidecar with no ``diffused_index_map`` -> an error status
  and NA everywhere. The tool does **not** fall back to a whole-structure
  alignment. Comparing two structures that have no externally supplied residue
  correspondence is ``usalign``'s job, which is a different premise.
* The correspondence is read, never re-derived. The contig string is not parsed
  and ``metrics.max_ca_deviation`` is not read: that field is a maximum rather
  than an RMSD and its reference frame is not documented anywhere this repo can
  check, so it is neither used nor trusted here.
* The numbers say nothing about whether the motif is *bindable*, whether the
  scaffold is designable, or whether the non-motif residues are sensible. They
  only say whether the motif's own geometry survived.

Every refusal is loud. A missing sidecar, an unparseable residue key, a map that
is not one-to-one (two reference residues sent to the same design residue), a
reference file that lacks a mapped chain/residue, fewer than 3 usable pairs, or a
count mismatch between the map and the atoms actually found all produce a status
starting with ``error:`` and NA for every numeric column -- never a plausible
number computed over a smaller or wrong atom set, and never a 0.0.

A segment is a run of motif residues contiguous on BOTH sides of the map --
consecutive in the reference AND consecutive in the design -- so two motif pieces
that abut in reference numbering but were placed apart by the contig count as two.

Per-segment RMSDs are reported after **one global superposition** shared by every
segment (never a per-segment refit): a loop that is locally perfect but globally
misplaced relative to the other loops is exactly the failure this tool has to be
able to show, and a per-segment refit would hide it.

The same correspondence also yields ``design_fixed_positions``: those segments
printed from the DESIGN side instead of the reference side, in ProteinMPNN's
``--fixed-positions`` syntax for one chain group, so a downstream sequence-design
run can hold exactly the grafted motif fixed. It is NA -- never a partial string,
and never with any other column or the status changed -- when the motif spans
more than one design chain, when that chain is not numbered 1..N in file order
(ProteinMPNN positions are 1-indexed chain POSITIONS, not residue numbers; the
two coincide only for a gapless chain numbered from 1), or when fewer motif
residues matched than the map names (an --allow-partial row: a spec over the
matched subset would leave the rest of the motif free to be redesigned).

``--motif-chain <letter>`` exists for the one case where that multi-chain refusal
fires on every design of a well-formed campaign: **binder scaffolding against a
target held fixed**. rfd3's ``diffused_index_map`` records correspondence for
EVERY held-fixed residue, so such a contig maps the whole target (chains A, B, C)
*and* the grafted motif on the binder chain (D) -- the motif spans four design
chains, and ``design_fixed_positions`` is correctly but uselessly NA. (Measured on
``outputs/20261009_130629_5M2J_vhh_redmotif_lowtemp_scaleup``: ``design_chains``
was ``A,B,C,D`` on all 5 designs.) Given a REFERENCE-frame chain letter, the pairs
feeding ``design_fixed_positions`` are filtered to that reference chain's residues
BEFORE the same ``segments()`` grouping runs, and the same checks are then applied
to the restricted result.

The flag changes **one column** (plus the design-side letter reported alongside
it in ``fixed_positions_chain``). ``rmsd``, ``max_dev``, ``seg_ranges``,
``seg_rmsds``, ``n_segments``, ``n_expected``, ``n_matched``, ``seq_match_frac``,
``design_chains`` and every other column keep measuring the WHOLE mapping, flag or
no flag -- so ``design_chains`` still honestly lists all four chains the RMSD
covers, and ``fixed_positions_chain`` is what names the one chain the
``--fixed-positions`` string applies to. A well-formed letter that names no
reference residue in this design's motif is a caller mistake in the flag, not a
measurement failure: that design's column is NA with a reason, and nothing else
about the row changes.

Writes a one-row TSV (name, status, aligned_path, then the metrics) for
collect_motifrmsd.py, and the design moved into the reference frame as
``<name>_onmotif.pdb`` so the graft can be eyeballed against the reference.

Normally invoked per array task by the motifrmsd tool -- run ``sapia run motifrmsd``
rather than calling this directly.

Usage (standalone, single design):
    python tools/motifrmsd/motifrmsd_worker.py \\
        --name design_0 \\
        --design design_0.pdb \\
        --sidecar results_shard_0/shard_0_x_diff_0_model_0.json \\
        --ref 5M2J_D_vhh.pdb \\
        --result-tsv design_0.tsv --out-pdb design_0_onmotif.pdb
"""

import argparse
import csv
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np

# Residue key format of ``diffused_index_map``: a chain id made of letters
# followed by a residue number, e.g. "D24" / "A23". Anything else (an insertion
# code, a numeric chain id, a blank) is refused rather than guessed at -- a
# mis-parsed key would silently score the wrong residue.
RES_KEY = re.compile(r"^([A-Za-z]+)(-?\d+)$")

# Kabsch needs three non-collinear points; below that the fit is meaningless.
MIN_PAIRS = 3

# Metric columns of the per-design TSV, in order. Bare names: collect_motifrmsd.py
# hands them to the driver, which leaf-prefixes them to motifrmsd_<name>.
# The per-segment seg<i>_rmsd columns are NOT here -- their number depends on the
# design, so they are appended per row (see ``seg_columns``).
METRIC_COLUMNS = [
    # --- trust: read these before any RMSD below ---
    "n_expected",
    "n_matched",
    "seq_match_frac",
    "ref_source",
    "ref_matches_spec",
    "ref_used",
    "spec_input",
    "ref_chains",
    "design_chains",
    "ref_motif_seq",
    "design_motif_seq",
    "sampled_contig",
    # --- science ---
    "rmsd",
    "max_dev",
    "n_segments",
    "worst_seg_rmsd",
    "seg_ranges",
    "seg_rmsds",
    "design_fixed_positions",
    "fixed_positions_chain",
]
RESULT_COLUMNS = ["name", "status", "aligned_path", *METRIC_COLUMNS]


# ---------------------------------------------------------------------------
# Structure bookkeeping
# ---------------------------------------------------------------------------


def _is_amino_acid(res: gemmi.Residue) -> bool:
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and info.is_amino_acid())


def one_letter(res: gemmi.Residue) -> str:
    """One-letter code of an amino acid, 'X' when gemmi does not know it."""
    info = gemmi.find_tabulated_residue(res.name)
    code = (info.one_letter_code if info else "").upper()
    return code if code.isalpha() else "X"


@dataclass
class CaIndex:
    """CA atoms of the first model, keyed by ``(chain name, residue number)``.

    Built once per structure; ``code[key]`` is the residue's one-letter code and
    ``coord[key]`` its CA position. Duplicate keys (insertion codes, a chain
    listed twice) are refused by the builder rather than resolved silently --
    "last one wins" would pick a residue nobody chose.
    """

    coord: dict[tuple[str, int], np.ndarray]
    code: dict[tuple[str, int], str]
    chains: list[str]


def read_structure(path: Path, what: str) -> gemmi.Structure:
    """Parse a PDB/mmCIF (optionally .gz) or raise with a message naming *what*."""
    if not path.is_file():
        raise ValueError(f"{what} not found: {path}")
    try:
        st = gemmi.read_structure(str(path))
    except Exception as e:  # noqa: BLE001 - reported as data by the caller
        raise ValueError(f"{what} {path} could not be parsed: {e}") from e
    if len(st) == 0:
        raise ValueError(f"{what} {path} holds no model")
    # Deterministic CA pick when a crystal structure carries altlocs.
    if hasattr(st, "remove_alternative_conformations"):
        st.remove_alternative_conformations()
    return st


def ca_index(st: gemmi.Structure, what: str) -> CaIndex:
    """Index the first model's amino-acid CA atoms by (chain, resnum)."""
    coord: dict[tuple[str, int], np.ndarray] = {}
    code: dict[tuple[str, int], str] = {}
    chains: list[str] = []
    duplicates: list[str] = []
    for chain in st[0]:
        chains.append(chain.name)
        for res in chain:
            if not _is_amino_acid(res):
                continue
            ca = next((a for a in res if a.name == "CA"), None)
            if ca is None:
                continue
            key = (chain.name, res.seqid.num)
            if key in coord:
                duplicates.append(f"{chain.name}{res.seqid.num}")
                continue
            coord[key] = np.array([ca.pos.x, ca.pos.y, ca.pos.z], dtype=float)
            code[key] = one_letter(res)
    if duplicates:
        preview = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(
            f"{what} has duplicate residue keys (insertion codes or a repeated "
            f"chain): {preview}. The motif correspondence cannot be resolved "
            f"unambiguously."
        )
    return CaIndex(coord, code, chains)


# ---------------------------------------------------------------------------
# The correspondence
# ---------------------------------------------------------------------------


def parse_res_key(token: str, where: str) -> tuple[str, int]:
    """``"D24"`` -> ``("D", 24)``. Raises on anything that is not that shape."""
    m = RES_KEY.match(token.strip())
    if not m:
        raise ValueError(
            f"{where}: residue key {token!r} is not <chain><resnum> (e.g. 'D24'); "
            f"insertion codes and numeric chain ids are not supported"
        )
    return m.group(1), int(m.group(2))


def read_sidecar(path: Path) -> dict:
    """The rfd3 per-design sidecar JSON, or a message saying why it is unusable."""
    # The manifest carries an empty field when the submit step could not resolve
    # a sidecar for this design; argparse turns that into Path('.').
    if str(path) in ("", "."):
        raise ValueError(
            "rfd3 sidecar JSON could not be located for this design at submit "
            "time (the run log names the route that failed). Without it there is "
            "no motif correspondence and nothing to measure."
        )
    if not path.is_file():
        raise ValueError(f"rfd3 sidecar JSON not found: {path}")
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        raise ValueError(f"rfd3 sidecar JSON {path} is unreadable: {e}") from e
    if not isinstance(data, dict):
        raise ValueError(f"rfd3 sidecar JSON {path} is not an object")
    return data


def index_map(sidecar: dict, path: Path) -> dict[str, str]:
    """``diffused_index_map`` -- the authoritative motif correspondence.

    Absent or empty is an error, not an empty measurement: a design with no
    recorded motif is a design this tool has nothing to say about.
    """
    raw = sidecar.get("diffused_index_map")
    if raw is None:
        raise ValueError(
            f"rfd3 sidecar {path} has no 'diffused_index_map' (top-level keys: "
            f"{sorted(sidecar)}). This is not motif-scaffolding output."
        )
    if not isinstance(raw, dict) or not raw:
        raise ValueError(f"rfd3 sidecar {path} has an empty/non-dict "
                         f"'diffused_index_map'")
    return {str(k): str(v) for k, v in raw.items()}


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Optimal rigid superposition of ``p`` onto ``q`` (both (N, 3), paired).

    Returns ``(rot, trans)`` with ``p @ rot.T + trans`` the moved ``p``.
    Reflections are excluded the usual way (flip the last singular vector when
    the naive rotation has determinant -1), so this is a rotation, never a
    mirror image that would flatter a wrong-handed design.
    """
    if p.shape != q.shape or p.shape[0] < MIN_PAIRS:
        raise ValueError(f"need >={MIN_PAIRS} paired points, got {p.shape[0]}")
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


def segments(
    ref_keys: list[tuple[str, int]], design_keys: list[tuple[str, int]]
) -> list[list[int]]:
    """Group matched pairs into runs that are contiguous on BOTH sides.

    ``ref_keys[i]`` and ``design_keys[i]`` are the two ends of pair ``i``. The
    grouping sorts by the reference ``(chain, resnum)`` first, so segment order is
    reference order and does not depend on the JSON's key order or on the contig.

    A pair extends the current run only when BOTH of these hold relative to the
    previous pair in that order:

    * reference side: same chain, and resnum exactly one higher;
    * design side: same chain, and resnum exactly one higher.

    Either check failing starts a new segment. The design-side check is what
    separates two motif pieces that happen to be numbered consecutively in the
    reference (e.g. ``D26-32`` and ``D33-40``) but that the contig placed apart in
    the design, with scaffold between them: they are two loops, and merging them
    would turn ``worst_seg_rmsd`` into a blended average that cannot isolate the
    bad one. Returns lists of indices into ``ref_keys``.
    """
    if len(ref_keys) != len(design_keys):
        raise ValueError(
            f"segments: {len(ref_keys)} reference keys but {len(design_keys)} "
            f"design keys -- the pairs are not aligned"
        )
    order = sorted(range(len(ref_keys)), key=lambda i: ref_keys[i])
    runs: list[list[int]] = []
    for i in order:
        chain, num = ref_keys[i]
        d_chain, d_num = design_keys[i]
        if runs:
            prev = runs[-1][-1]
            p_chain, p_num = ref_keys[prev]
            pd_chain, pd_num = design_keys[prev]
            ref_contiguous = p_chain == chain and num == p_num + 1
            design_contiguous = pd_chain == d_chain and d_num == pd_num + 1
            if ref_contiguous and design_contiguous:
                runs[-1].append(i)
                continue
        runs.append([i])
    return runs


def _mpnn_indexable(chain: str, design_ca: CaIndex) -> str:
    """``""`` when this design chain's residue NUMBERS are also ProteinMPNN's
    1-indexed chain POSITIONS; otherwise the reason they are not.

    ProteinMPNN's ``--position_list`` is positional: ``make_fixed_positions_dict.py``
    takes the i-th chain of ``--chain_list`` and treats each number as a 1-based
    index into that chain's parsed sequence, NOT as a PDB residue number
    (prosapia's own position mini-language says the same -- "positions are
    1-indexed within each chain"). The two coincide exactly when the chain is
    numbered ``1..N`` in file order with no gaps, which is what rfd3 writes
    because it renumbers every output chain from 1.

    So this is the check that makes ``design_fixed_positions`` safe to paste into
    a ProteinMPNN flag. If it fails the column is NA rather than a plausible
    string, because an offset or gapped chain would fix the WRONG residues --
    silently, at exit 0, with a redesigned motif.
    """
    nums = [n for (c, n) in design_ca.coord if c == chain]
    if not nums:
        return f"design chain {chain!r} holds no CA atoms"
    if nums != list(range(1, len(nums) + 1)):
        return (
            f"design chain {chain!r} is numbered {nums[0]}..{nums[-1]} over "
            f"{len(nums)} CA residues rather than 1..{len(nums)} in file order, "
            f"so its residue numbers are not ProteinMPNN's 1-indexed chain "
            f"positions"
        )
    return ""


def restrict_to_ref_chain(
    ref_keys: list[tuple[str, int]],
    design_keys: list[tuple[str, int]],
    chain: str,
) -> tuple[list[list[int]], list[tuple[str, int]]]:
    """The matched pairs whose REFERENCE residue lies on *chain*, segmented.

    Returns ``(runs, design_keys)`` for that subset alone, re-indexed from 0, in
    exactly the shape ``design_fixed_positions`` consumes.

    The filtering happens BEFORE grouping and the grouping is the shared
    ``segments()`` -- not a second implementation of contiguity. Dropping the
    other chains first is what makes the result right rather than merely smaller:
    the pairs that remain are regrouped from scratch, so a run is never left
    straddling a hole where a filtered-out chain used to be.

    Only ``design_fixed_positions`` is computed from this. ``seg_ranges``,
    ``seg_rmsds``, ``n_segments`` and the RMSDs are all computed from the
    UNFILTERED ``segments(ref_keys, design_keys)``, whatever this returns.
    """
    keep = [i for i, (c, _) in enumerate(ref_keys) if c == chain]
    ref_sub = [ref_keys[i] for i in keep]
    design_sub = [design_keys[i] for i in keep]
    return segments(ref_sub, design_sub), design_sub


def design_fixed_positions(
    runs: list[list[int]],
    design_keys: list[tuple[str, int]],
    design_ca: CaIndex,
) -> tuple[str | None, str | None, str]:
    """The motif's DESIGN-side residue ranges as a ProteinMPNN ``--fixed-positions``
    value for ONE chain group, or ``(None, None, reason)``.

    Returns ``(spec, design_chain, reason)``: the second element is the design-side
    chain letter the spec addresses -- the value to pass as
    ``--chains-to-design`` -- so a caller never has to guess which of the
    ``design_chains`` letters it belongs to. It is non-None exactly when ``spec``
    is.

    Same segments as ``seg_ranges`` -- this consumes the very ``runs`` that
    ``segments()`` produced, so the contiguity rule (consecutive on BOTH the
    reference and the design side) is shared, not re-implemented. The only
    difference is which end of each pair is printed: ``seg_ranges`` prints the
    reference residues, this prints the design's own output residues.

    Format: ``"28:35,50:59,98:103"`` -- inclusive ``start:end`` fragments,
    comma-joined, **no trailing slash**. In prosapia's position mini-language
    ``/`` breaks chains and the groups map one-to-one onto ``--chains-to-design``,
    so for a single designed chain a bare group is the whole value; a trailing
    slash would declare a second, empty chain group. Fragment *i* is segment *i*,
    in the same reference order as ``seg_ranges``/``seg_rmsds``, so the three
    columns can be read side by side.

    ``(None, reason)`` -- the column is NA, nothing else about the row changes --
    when the string could not be made honestly:

    * the motif -- after any ``--motif-chain`` restriction the caller asked for --
      maps into more than one design chain, so no single-chain group can
      describe it;
    * the design chain is not numbered ``1..N`` in file order, so residue numbers
      are not ProteinMPNN positions (see ``_mpnn_indexable``);
    * a segment's design residue span does not cover exactly its mapped residues
      (an invariant of ``segments()``; checked rather than assumed).

    The caller, ``compute()``, adds one more NA condition before calling this:
    an incomplete match (``n_matched < n_expected``, only reachable under
    ``--allow-partial``), since ``runs`` then describes only part of the motif.
    """
    chains = {design_keys[i][0] for run in runs for i in run}
    if len(chains) != 1:
        return None, None, (
            f"the motif maps into {len(chains)} design chains "
            f"({', '.join(sorted(chains))}); a single --fixed-positions chain "
            f"group cannot describe that"
            + ("" if chains else " (no mapped residues remain)")
            + ". Pass --motif-chain <reference chain letter> to restrict this "
            "column to the residues lifted from one reference chain"
        )
    chain = chains.pop()
    # The 1..N check is about the WHOLE design chain as it stands in the file --
    # an absolute property of that chain, not of the motif subset -- so it is
    # applied unchanged to whichever single chain the (possibly filtered)
    # residues landed on.
    if why := _mpnn_indexable(chain, design_ca):
        return None, None, why

    fragments: list[str] = []
    for run in runs:
        start = design_keys[run[0]][1]
        end = design_keys[run[-1]][1]
        if end - start + 1 != len(run):
            return None, None, (
                f"design segment {chain}{start}-{end} spans "
                f"{end - start + 1} positions but holds {len(run)} mapped "
                f"residues; the segment grouping is not contiguous as assumed"
            )
        fragments.append(f"{start}:{end}")
    return ",".join(fragments), chain, ""


# ---------------------------------------------------------------------------
# Reference selection
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def choose_reference(
    ref_flag: str, sidecar: dict
) -> tuple[Path, str, int | None, str]:
    """Decide which reference file to measure against.

    Returns ``(path, source, matches_spec, spec_input)``.

    **``--ref-structure`` wins when given** -- an explicit flag is the caller's
    stated intent, and the path rfd3 recorded may no longer exist or may be a
    stale copy. The sidecar's own ``specification.input`` is always reported as
    ``spec_input``, and the two are compared by **content** (sha256), not by path
    string, because the same file reaches a task under several spellings:

        1  the flag and the sidecar's input are the same bytes -- the normal case
        0  they are different files -- the RMSD is measured against something
           rfd3 did not diffuse, so it answers a different question than it looks
           like it does (the status is downgraded to ``warn:``)
        NA the sidecar's input is not readable from here, so the question cannot
           be settled -- NOT an implied agreement

    With no ``--ref-structure``, the sidecar's input is used and ``matches_spec``
    is trivially 1.
    """
    spec = sidecar.get("specification")
    spec_input = ""
    if isinstance(spec, dict) and spec.get("input") is not None:
        spec_input = str(spec["input"])

    if not ref_flag:
        if not spec_input:
            raise ValueError(
                "no --ref-structure given and the sidecar records no "
                "specification.input, so there is no reference to measure against"
            )
        path = Path(spec_input)
        if not path.is_file():
            raise ValueError(
                f"no --ref-structure given and the sidecar's specification.input "
                f"{path} is not readable from this task"
            )
        return path, "sidecar", 1, spec_input

    path = Path(ref_flag)
    if not path.is_file():
        raise ValueError(f"--ref-structure not found: {path}")

    matches: int | None = None
    if spec_input:
        spec_path = Path(spec_input)
        if spec_path.is_file():
            try:
                matches = int(_sha256(spec_path) == _sha256(path))
            except OSError:
                matches = None
    return path, "flag", matches, spec_input


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------


@dataclass
class Result:
    metrics: dict[str, object]
    seg_rmsds: list[float]
    status: str
    structure: gemmi.Structure | None


def compute(args: "Args") -> Result:
    """Every motifrmsd metric for one design, plus the design moved into the
    reference frame. Raises on anything unusable (the caller records the message
    as data)."""
    sidecar = read_sidecar(args.sidecar)
    mapping = index_map(sidecar, args.sidecar)
    ref_path, ref_source, matches_spec, spec_input = choose_reference(
        args.ref, sidecar
    )

    extra = {}
    spec = sidecar.get("specification")
    if isinstance(spec, dict) and isinstance(spec.get("extra"), dict):
        extra = spec["extra"]
    sampled_contig = str(extra.get("sampled_contig", "") or "")

    design_st = read_structure(args.design, "design structure")
    ref_st = read_structure(ref_path, "reference structure")
    design_ca = ca_index(design_st, f"design structure {args.design}")
    ref_ca = ca_index(ref_st, f"reference structure {ref_path}")

    n_expected = len(mapping)
    ref_keys: list[tuple[str, int]] = []
    design_keys: list[tuple[str, int]] = []
    missing_ref: list[str] = []
    missing_design: list[str] = []

    # Every parsed key -> the map tokens that produced it, over the WHOLE map
    # (matched or not). Used for the one-to-one check below.
    ref_seen: dict[tuple[str, int], list[str]] = {}
    design_seen: dict[tuple[str, int], list[str]] = {}

    for ref_token, design_token in mapping.items():
        r_key = parse_res_key(ref_token, f"diffused_index_map key in {args.sidecar}")
        d_key = parse_res_key(
            design_token, f"diffused_index_map value in {args.sidecar}"
        )
        ref_seen.setdefault(r_key, []).append(ref_token)
        design_seen.setdefault(d_key, []).append(ref_token)
        ok = True
        if r_key not in ref_ca.coord:
            missing_ref.append(ref_token)
            ok = False
        if d_key not in design_ca.coord:
            missing_design.append(design_token)
            ok = False
        if ok:
            ref_keys.append(r_key)
            design_keys.append(d_key)

    # The map must be ONE-TO-ONE. A design residue named by two reference
    # residues (or two key spellings that parse to the same reference residue,
    # e.g. "D24"/"D024") would put one CA into the fit twice: every count check
    # below still passes, and the RMSD is silently wrong. Not overridable by
    # --allow-partial -- this is a broken correspondence, not a missing residue.
    dup_design = {k: v for k, v in design_seen.items() if len(v) > 1}
    dup_ref = {k: v for k, v in ref_seen.items() if len(v) > 1}
    if dup_design or dup_ref:
        detail = []
        if dup_design:
            detail.append(
                "design residue(s) mapped more than once: "
                + "; ".join(
                    f"{c}{n} <- {', '.join(toks)}"
                    for (c, n), toks in sorted(dup_design.items())[:4]
                )
                + ("..." if len(dup_design) > 4 else "")
            )
        if dup_ref:
            detail.append(
                "reference residue(s) keyed more than once: "
                + "; ".join(
                    f"{c}{n} <- {', '.join(toks)}"
                    for (c, n), toks in sorted(dup_ref.items())[:4]
                )
                + ("..." if len(dup_ref) > 4 else "")
            )
        raise ValueError(
            f"diffused_index_map in {args.sidecar} is not one-to-one: "
            + " | ".join(detail)
        )

    n_matched = len(ref_keys)
    warnings: list[str] = []

    # A count mismatch is a mapping failure, not a smaller measurement. Refuse by
    # default; --allow-partial is an explicit, recorded opt-in that still keeps
    # the row out of any `== "OK"` selection.
    if n_matched < n_expected:
        detail = []
        if missing_ref:
            detail.append(
                f"{len(missing_ref)} absent from reference {ref_path.name} "
                f"({', '.join(missing_ref[:4])}{'...' if len(missing_ref) > 4 else ''})"
            )
        if missing_design:
            detail.append(
                f"{len(missing_design)} absent from design "
                f"({', '.join(missing_design[:4])}{'...' if len(missing_design) > 4 else ''})"
            )
        message = (
            f"matched {n_matched} of {n_expected} motif residues: "
            + "; ".join(detail)
        )
        if not args.allow_partial:
            raise ValueError(message + " (pass --allow-partial to measure the subset)")
        warnings.append(message)

    if n_matched < MIN_PAIRS:
        raise ValueError(
            f"only {n_matched} motif CA pair(s) resolved in both structures; "
            f"a superposition needs at least {MIN_PAIRS}"
        )

    ref_xyz = np.array([ref_ca.coord[k] for k in ref_keys], dtype=float)
    design_xyz = np.array([design_ca.coord[k] for k in design_keys], dtype=float)

    # One global fit of the design's motif CA onto the reference's motif CA.
    rot, trans = kabsch(design_xyz, ref_xyz)
    moved = design_xyz @ rot.T + trans

    global_rmsd = rmsd(moved, ref_xyz)
    deviations = np.sqrt(((moved - ref_xyz) ** 2).sum(axis=1))
    max_dev = float(deviations.max())

    # Per-segment RMSDs, all from that SAME global fit (never a per-segment refit).
    runs = segments(ref_keys, design_keys)
    seg_rmsds = [rmsd(moved[idx], ref_xyz[idx]) for idx in map(np.asarray, runs)]
    seg_ranges = [
        f"{ref_keys[idx[0]][0]}{ref_keys[idx[0]][1]}-{ref_keys[idx[-1]][1]}"
        for idx in runs
    ]

    # The SAME segments, printed from the design side and in ProteinMPNN's own
    # --fixed-positions syntax, so a downstream run can hold exactly the grafted
    # motif fixed. NA (never a partial or guessed string) when it cannot be made
    # honestly -- the reason goes to stderr; no other column and no status is
    # affected by it.
    #
    # It requires a COMPLETE match. Under --allow-partial the other columns are
    # measured on the matched subset, but a fixed-positions spec covering only
    # that subset would leave the unmatched motif residues free to be redesigned
    # -- the opposite of this column's purpose -- so it is NA instead.
    # --motif-chain restricts THIS COLUMN ONLY (see the module docstring): the
    # pairs are filtered to one reference chain before the shared segments()
    # regroups them. Every column above is already computed, over the full map,
    # and nothing below touches them.
    motif_chain = (args.motif_chain or "").strip()
    map_ref_chains = sorted({c for c, _ in ref_seen})
    fixed_chain: str | None = None
    if motif_chain and motif_chain not in map_ref_chains:
        # A caller mistake in the flag, not a measurement failure: NA with a
        # reason, the status and every other column untouched.
        fixed_positions: str | None = None
        fixed_reason = (
            f"--motif-chain {motif_chain} given but no reference residue in "
            f"chain {motif_chain} is part of this design's motif "
            f"(diffused_index_map names reference chain(s) "
            f"{', '.join(map_ref_chains)})"
        )
    elif n_matched < n_expected:
        fixed_positions = None
        fixed_reason = (
            f"only {n_matched} of {n_expected} motif residues matched "
            f"(--allow-partial); a spec over the matched subset would leave "
            f"{n_expected - n_matched} motif residue(s) free to be redesigned"
        )
    else:
        fp_runs, fp_design_keys = (
            restrict_to_ref_chain(ref_keys, design_keys, motif_chain)
            if motif_chain
            else (runs, design_keys)
        )
        fixed_positions, fixed_chain, fixed_reason = design_fixed_positions(
            fp_runs, fp_design_keys, design_ca
        )
    if fixed_positions is None:
        print(
            f"design_fixed_positions not emitted (NA): {fixed_reason}",
            file=sys.stderr,
        )

    # Residue identity over the matched pairs, in segment order: the trust column
    # that catches a right-looking map pointing at the wrong residues.
    seq_order = [i for run in runs for i in run]
    ref_seq = "".join(ref_ca.code[ref_keys[i]] for i in seq_order)
    design_seq = "".join(design_ca.code[design_keys[i]] for i in seq_order)
    seq_match_frac = sum(a == b for a, b in zip(ref_seq, design_seq)) / len(ref_seq)

    # Measuring against a file rfd3 did not diffuse answers a different question
    # than the column name suggests, so it never reads as OK.
    if matches_spec == 0:
        warnings.append(
            f"--ref-structure {ref_path.name} differs in content from the "
            f"sidecar's specification.input {Path(spec_input).name}"
        )

    status = "OK" if not warnings else "warn: " + "; ".join(warnings)

    metrics: dict[str, object] = {
        "n_expected": n_expected,
        "n_matched": n_matched,
        "seq_match_frac": seq_match_frac,
        "ref_source": ref_source,
        "ref_matches_spec": matches_spec,
        "ref_used": str(ref_path),
        "spec_input": spec_input,
        "ref_chains": ",".join(sorted({k[0] for k in ref_keys})),
        "design_chains": ",".join(sorted({k[0] for k in design_keys})),
        "ref_motif_seq": ref_seq,
        "design_motif_seq": design_seq,
        "sampled_contig": sampled_contig,
        "rmsd": global_rmsd,
        "max_dev": max_dev,
        "n_segments": len(runs),
        "worst_seg_rmsd": max(seg_rmsds),
        "seg_ranges": "/".join(seg_ranges),
        "seg_rmsds": "/".join(f"{v:.3f}" for v in seg_rmsds),
        "design_fixed_positions": fixed_positions,
        # The design chain design_fixed_positions addresses -- i.e. the value to
        # pass as --chains-to-design. NA exactly when that column is NA. NOT the
        # same as design_chains, which keeps listing every chain the RMSD covers.
        "fixed_positions_chain": fixed_chain,
    }

    # Move the WHOLE design into the reference frame, so the graft can be loaded
    # next to the reference and looked at. Applied with the same rot/trans the
    # RMSD was computed from -- plain numpy, so there is one transform in play
    # and no second convention to get wrong.
    for model in design_st:
        for chain in model:
            for res in chain:
                for atom in res:
                    v = np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
                    x, y, z = rot @ v + trans
                    atom.pos = gemmi.Position(float(x), float(y), float(z))

    return Result(metrics, seg_rmsds, status, design_st)


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
    return str(value)


def write_result(
    result_tsv: Path,
    name: str,
    status: str,
    aligned_path: str,
    metrics: dict,
    extra_columns: dict[str, float] | None = None,
) -> None:
    """One-row TSV. Error rows carry the fixed columns empty (NA), never 0."""
    extra_columns = extra_columns or {}
    columns = [*RESULT_COLUMNS, *extra_columns]
    result_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(result_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)
        writer.writerow(
            [name, status, aligned_path]
            + [_fmt(metrics.get(c)) for c in METRIC_COLUMNS]
            + [_fmt(v) for v in extra_columns.values()]
        )


class Args(argparse.Namespace):
    name: str
    design: Path
    sidecar: Path
    ref: str
    allow_partial: bool
    motif_chain: str
    out_pdb: Path
    result_tsv: Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Motif CA RMSD of one rfd3 motif-scaffolding design against "
        "its reference, over rfd3's own diffused_index_map correspondence."
    )
    ap.add_argument("--name", required=True)
    ap.add_argument(
        "--design", type=Path, required=True, help="Design structure (PDB or mmCIF)."
    )
    ap.add_argument(
        "--sidecar",
        type=Path,
        required=True,
        help="The design's rfd3 sidecar JSON (carries diffused_index_map).",
    )
    ap.add_argument(
        "--ref",
        default="",
        help="Reference structure the motif was taken from. Empty (default) "
        "falls back to the sidecar's specification.input.",
    )
    ap.add_argument(
        "--allow-partial",
        action="store_true",
        help="Measure the matched subset when some mapped residues are absent "
        "from either structure, with a 'warn:' status. Off by default: a count "
        "mismatch is a mapping failure, not a smaller measurement.",
    )
    ap.add_argument(
        "--motif-chain",
        default="",
        help="Restrict design_fixed_positions to the motif residues lifted from "
        "this ONE REFERENCE chain (the key side of diffused_index_map), e.g. 'D'. "
        "Affects that column and fixed_positions_chain only -- every other column "
        "keeps measuring the whole mapping. Empty (default) = unrestricted.",
    )
    ap.add_argument(
        "--out-pdb", type=Path, required=True, help="Superposed design to write."
    )
    ap.add_argument(
        "--result-tsv", type=Path, required=True, help="Per-design result TSV to write."
    )
    args = ap.parse_args(namespace=Args())

    try:
        result = compute(args)
        assert result.structure is not None
        args.out_pdb.parent.mkdir(parents=True, exist_ok=True)
        result.structure.setup_entities()
        result.structure.write_pdb(str(args.out_pdb))
        m = result.metrics
        print(
            f"{args.name}: rmsd {m['rmsd']:.3f} A over {m['n_matched']}/"
            f"{m['n_expected']} motif CA, max_dev {m['max_dev']:.3f} A, "
            f"segments {m['seg_ranges']} = {m['seg_rmsds']}, "
            f"seq_match_frac {m['seq_match_frac']:.3f}, ref_source {m['ref_source']}"
        )
        write_result(
            args.result_tsv,
            args.name,
            result.status,
            str(args.out_pdb),
            result.metrics,
            seg_columns(result.seg_rmsds),
        )
    except Exception as e:  # noqa: BLE001 - errors are recorded as data
        print(f"{args.name}: ERROR {e}", file=sys.stderr)
        write_result(args.result_tsv, args.name, f"error: {e}", "", {})


if __name__ == "__main__":
    main()
