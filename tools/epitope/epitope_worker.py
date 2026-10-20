#!/usr/bin/env python3
"""
Which target residues does this binder actually touch, and did it hit the epitope?

This is the per-array-task step of the epitope tool (and works standalone). One task
scores several designs. For each design it counts heavy-atom contacts between
``--design-chains`` (the binder) and every residue of ``--target-chains`` (the
target), then reports the epitope that was ACTUALLY used and how much of the
REQUESTED one it covers.

It measures CONTACT and nothing else -- no buried surface area (that is cms), no
energy (that is pyrosetta), no assembly fit (that is ringfit), and no judgement at
all about whether the coordinates it was handed are trustworthy.

DEFINITIONS, all at ``--contact-cutoff`` A between heavy atoms (hydrogens and
waters excluded; a ligand sitting in a named chain counts as part of that chain):

    contact pair        one design heavy atom within the cutoff of one target
                        heavy atom. n_contacts is the total over the interface.
    interface residue   a target residue with >= 1 contact pair. The list of them
                        IS the epitope the binder used.
    hotspot hit         a hotspot that is an interface residue.
    hotspot_recall      hits / listed hotspots.
    hotspot_n_contacts  contact pairs on the LISTED hotspots only -- the same
                        pairs, the same cutoff, just summed over fewer residues.
    hotspot_contact_frac
                        hotspot_n_contacts / n_contacts: the SHARE of the whole
                        interface the listed residues carry, 0.0-1.0. Undefined
                        (NA) when n_contacts is 0 -- never 0.0.
    hotspot_contacts    the per-hotspot pair counts in the order listed,
                        'A24=0,A36=3,A99=21'.

STERIC OVERLAP (added 2026-10-08). The same heavy atoms, the same two chain sets,
the same residue universe (hydrogens and waters out), read for overlap instead of
contact:

    min_dist_any        the smallest heavy-atom distance between ANY design atom
                        and ANY target atom, over all target residues -- the
                        minimum of the per-residue file's min_dist column, so it
                        is always <= min_dist_hotspot.
    n_clash             design/target heavy-atom PAIRS strictly closer than
                        --clash-cutoff (default 2.5 A -- the same cutoff VALUE as
                        ringfit's default, but NOT the same count: the unit here is
                        atom PAIRS, while ringfit_n_clash counts binder ATOMS, so
                        the numbers are never comparable).
                        Non-bonded atoms of different chains below 2.5 A overlap;
                        hydrogen-bond and salt-bridge distances (2.5-3.0 A) are not
                        counted. 0 is a real, measured zero on an OK row, and, when
                        both are numbers, n_clash == 0 exactly when min_dist_any >=
                        the cutoff (both are NA together when the target chains hold
                        no heavy atoms).
    n_clash_design_res  DISTINCT design residues with >= 1 clash pair (residues,
                        not atoms) -- the same quantity as ringfit_n_clash_res.
    clash_cutoff        the echo of the cutoff used.
    min_dist_bb         as min_dist_any, but with the DESIGN side restricted to its
                        BACKBONE atoms (N, CA, C, O, OXT of amino-acid residues);
                        the target side is still every heavy atom. Always >=
                        min_dist_any.
    n_clash_bb          as n_clash over those pairs; always <= n_clash, and 0
                        exactly when min_dist_bb >= the cutoff (when both are
                        numbers). Both NA -- not 0 -- when the design chains have no
                        backbone-named atoms. On a generator backbone the design's
                        side chains are generator artifacts that sequence design
                        replaces, while its backbone persists and the target's side
                        chains are real and fixed: backbone-vs-any-target-atom is
                        the overlap no sequence design can remove. On a
                        backbone-only file these equal min_dist_any and n_clash.
                        CB is NOT in this set -- this is the FLOOR (lenient).
    min_dist_bbcb       as min_dist_bb with CB added to the design set (N, CA, C,
    n_clash_bbcb        O, OXT, CB). CB's position is fixed by N, CA and C, so a
                        CB overlap persists under every sequence except GLY -- a
                        real design cost at an interface. The honest persistent
                        set: min_dist_any <= min_dist_bbcb <= min_dist_bb and
                        n_clash_bb <= n_clash_bbcb <= n_clash; n_clash_bbcb is 0
                        exactly when min_dist_bbcb >= the cutoff (when both are
                        numbers). NA when the design chains have no atom of the
                        set. With no CB anywhere (poly-GLY) they equal the bb pair.

This is geometry: not an energy, and not a statement that the target (or the
binder) is valid. A clash count is only as good as the atoms present -- a
backbone-only binder has no side chains to clash with.

WEIGHT IS NOT PRESENCE. hotspot_recall answers 'did it touch the epitope'; the
three columns above answer 'how much of the binding is on it'. On a target where
nearly every binder grazes one residue, recall discriminates nothing and the share
does. Read the share as a SHARE: it is a fraction of this design's own interface,
so a 300-pair interface and an 80-pair one can carry the same 0.15, and a rising
share can mean the hotspot got more buried OR that the rest of the interface got
smaller. It is also a count of ATOM PAIRS, so a large, deeply inserted side chain
(ARG, TRP, a histidine wedged into a pocket) contributes more pairs than a small
one at the same burial -- the number is weighted by atom count, not by energy.

BOTH AXES OF THE SAME MATRIX. The contact matrix has a target side and a design
side, and both are reported. The target side answers "which epitope did it use"
(n_iface_target_res, iface_target_resnums, the hotspot_* columns); the design side
answers "which of MY residues are doing the binding" (n_iface_design_res,
iface_design_resnums, iface_design_resnames). Same cutoff, same matrix, one pass.
``iface_design_resnames`` is to the binder what ``hotspot_resnames`` is to the
target: 'B12=HIS,B15=TYR,B19=GLU', eyeballable, and a frame shift is obvious in it.

WHERE THE DESIGN-CHAIN IDENTITIES COME FROM. By default, from the coordinates.
With ``--seq-source <column name>`` the task file carries a one-chain sequence per
design and the identities are taken from THAT, mapped positionally onto the design
chains' amino-acid residues in structure order. The contact GEOMETRY always comes
from the structure. The reason: a sequence designer writes its sequence into a
CHILD table while the backbone that defines the geometry lives in the PARENT, and
that backbone's own residue names are generator artifacts, not the designed
sequence. Two hard contracts, because a positional mapping that is off by one
misreports the interface composition of every row with no other symptom:

    len(sequence) != <design-chain amino-acid residues>
                        -> 'error: ...' naming BOTH numbers. Never truncate,
                           never pad, never score a partial mapping.
    a '/' in the sequence
                        -> 'error: ...'. This tool wants the design chain alone;
                           a multi-chain string means the wrong column was named.

``seq_len``, ``n_design_res_struct`` and ``seq_source`` record what was actually
used, so the table always says where the residue names came from.

THE NUMBERING CONTRACT. A hotspot label that is not a residue of the target chains
is an ERROR for that design (status 'error: ...', every column NA), never a recall
of 0.0. The error message reports the residue-number span actually present in each
target chain, because the cause is nearly always a renumbered chain. On success,
``hotspot_resnames`` records the residue NAME found at each listed position, in the
order listed -- the cheap way to see that A96 really is the lysine you meant.

NA IS NEVER 0. A design that failed has NA columns. A binder that touched the
target but missed every hotspot has ``hotspot_recall 0.0`` and ``hotspot_hits
'none'``; a binder that touched nothing at all has ``n_iface_target_res 0`` and
``iface_target_resnums 'none'``. The string 'none' is deliberate: an empty cell
would be read back as NA, and "measured, and the answer is nothing" must stay
distinguishable from "not measured".

Writes, per design, ``<name>.tsv`` (one row: name, status, path, metrics...) and,
on success, ``<name>_per_residue.tsv`` -- EVERY target residue, contacted or not,
with chain, resnum, resname, n_contacts, min_dist and is_hotspot. The full target is
written on purpose: it is the per-design record in which a numbering mismatch or a
second, unintended binding site is visible at a glance.

Errors are recorded as data (status 'error: ...') so a partial array still collects.

Normally invoked per array task by the epitope tool -- run ``sapia run epitope``
rather than calling this directly.

Usage (standalone):
    printf 'design_0\\tcomplex.pdb\\tA96,A99,A101,A155\\t\\n' > task.tsv
    python tools/epitope/epitope_worker.py --task-file task.tsv \\
        --design-chains B --target-chains A --contact-cutoff 5.0 --out-dir out/

    # identities from a designed sequence instead of the coordinates:
    printf 'design_0\\tcomplex.pdb\\tA96,A99\\tMHHKLEDG...\\n' > task.tsv
    python tools/epitope/epitope_worker.py --task-file task.tsv \\
        --design-chains B --target-chains A --seq-source atomium_sequence \\
        --contact-cutoff 5.0 --out-dir out/
"""

import argparse
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}
# A canonical hotspot label as the manifest writes it: chain + number + optional
# insertion code. Ranges were already expanded at submit time (run_epitope.py).
_LABEL = re.compile(r"^([A-Za-z])(\d+)([A-Za-z]?)$")

# The --seq-source value meaning "identities came from the coordinates" -- the
# default, and exactly what the tool did before --sequence-column existed.
SEQ_SOURCE_STRUCTURE = "structure"

# One-letter -> three-letter, so a sequence-derived identity is formatted exactly
# like a structure-derived one and iface_design_resnames reads the same either way
# ('B12=HIS'). A letter outside the standard 20 is kept AS the letter rather than
# guessed into a three-letter code -- visible, and never silently a real residue.
ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

# Metric columns of the per-design TSV, in order (bare names: collect_epitope.py
# hands them to the driver, which leaf-prefixes them to epitope_<name>).
#
# APPEND-ONLY. Live tables in this campaign carry epitope_bb_*, epitope_bbnterm_*
# and epitope_pred_* columns from earlier runs; nothing above the marker may be
# renamed, reordered out of existence or change meaning, or those numbers stop
# matching the ones already relied on.
METRIC_COLUMNS = [
    "hotspot_recall",
    "hotspot_hits",
    "n_hotspots",
    "hotspot_resnames",
    "n_iface_target_res",
    "iface_target_resnums",
    "min_dist_hotspot",
    "n_contacts",
    "design_chains",
    "target_chains",
    "contact_cutoff",
    # --- added 2026-09-30: the design side of the same contact matrix ---
    "n_iface_design_res",
    "iface_design_resnums",
    "iface_design_resnames",
    # --- added 2026-09-30: where the design identities came from ---
    "seq_len",
    "n_design_res_struct",
    "seq_source",
    # --- added 2026-10-02: HOW MUCH of the interface sits on the listed
    # hotspots, not merely whether they were touched ---
    "hotspot_n_contacts",
    "hotspot_contact_frac",
    "hotspot_contacts",
    # --- added 2026-10-08: STERIC OVERLAP between the two sides, from the same
    # heavy atoms. Geometry only -- not an energy, not a validity statement ---
    "min_dist_any",
    "n_clash",
    "n_clash_design_res",
    "clash_cutoff",
    # --- added 2026-10-08 (second pass): the DESIGN BACKBONE against every target
    # heavy atom -- the overlap sequence design cannot remove ---
    "min_dist_bb",
    "n_clash_bb",
    # --- added 2026-10-08 (third pass): backbone PLUS CB -- CB is fixed by N, CA,
    # C, so its overlap persists under every sequence except GLY ---
    "min_dist_bbcb",
    "n_clash_bbcb",
]

# Design-side atoms that survive sequence design.
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
# ... plus CB, whose position the backbone fixes (removable only by choosing GLY).
BACKBONE_CB_ATOMS = frozenset(BACKBONE_ATOMS | {"CB"})
RESULT_COLUMNS = ["name", "status", "path", *METRIC_COLUMNS]

# Default --clash-cutoff (A): the same cutoff VALUE as ringfit's default (but a
# different unit: epitope counts atom PAIRS, ringfit_n_clash counts binder atoms). A
# pair of heavy atoms from DIFFERENT chains closer than this overlaps; hydrogen-bond
# and salt-bridge distances (2.5-3.0 A) are deliberately NOT counted.
DEFAULT_CLASH_CUTOFF = 2.5
PER_RESIDUE_COLUMNS = [
    "chain",
    "resnum",
    "resname",
    "n_contacts",
    "min_dist",
    "is_hotspot",
]
# "measured, and the answer is an empty list" -- never a blank cell, which collect
# reads back as NA.
NONE = "none"


def split_list(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


# ---------------------------------------------------------------------------
# Structure bookkeeping
# ---------------------------------------------------------------------------


@dataclass
class Res:
    """One residue of the target, flattened for numpy."""

    chain: str
    resnum: str  # seqid number + insertion code, e.g. "96" or "96A"
    resname: str
    heavy: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))

    @property
    def label(self) -> str:
        """Residue key in the ``A96`` style used by --hotspots."""
        return f"{self.chain}{self.resnum}"

    @property
    def is_amino_acid(self) -> bool:
        return _is_amino_acid(self.resname)

    @property
    def sort_key(self) -> tuple[str, int, str]:
        """(chain, number, insertion code) -- ascending residue order, which is what
        the iface_design_* columns are reported in so that resnums and resnames stay
        index-aligned and a human can read one against the other."""
        digits = re.sub(r"[^0-9-]", "", self.resnum)
        icode = re.sub(r"[0-9-]", "", self.resnum)
        return (self.chain, int(digits) if digits else 0, icode)


def _is_heavy(atom: gemmi.Atom) -> bool:
    """Non-hydrogen. Tested on the element NAME because ``Element.is_hydrogen`` is a
    method in some gemmi versions and a property in others (0.7.x)."""
    return atom.element.name.upper() not in ("H", "D")


def _is_water(res: gemmi.Residue) -> bool:
    return res.name.strip().upper() in WATER_NAMES


def _is_amino_acid(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_amino_acid())


def read_model(path: Path) -> gemmi.Model:
    """First model of a PDB or mmCIF file, one altloc per atom, no hydrogens."""
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        raise ValueError(f"no model in {path.name}")
    return st[0]


def chain_ids(model: gemmi.Model) -> list[str]:
    return [chain.name for chain in model]


def residues_of(model: gemmi.Model, chain_names: list[str]) -> list[Res]:
    """Non-water residues of the named chains, in the order the chains are named and
    the residues appear in the file. Ligands and modified residues are kept: a
    hotspot may legitimately point at one, and a ligand at the interface is a real
    contact partner."""
    out: list[Res] = []
    for name in chain_names:
        for chain in model:
            if chain.name != name:
                continue
            for res in chain:
                if _is_water(res):
                    continue
                coords = [
                    [a.pos.x, a.pos.y, a.pos.z] for a in res if _is_heavy(a)
                ]
                icode = res.seqid.icode.strip()
                out.append(
                    Res(
                        chain=chain.name,
                        resnum=f"{res.seqid.num}{icode}",
                        resname=res.name.strip(),
                        heavy=np.asarray(coords, dtype=float).reshape(-1, 3),
                    )
                )
    return out


def heavy_atoms(residues: list[Res]) -> np.ndarray:
    """(M, 3) heavy atoms of the given residues."""
    blocks = [r.heavy for r in residues if len(r.heavy)]
    if not blocks:
        return np.empty((0, 3), dtype=float)
    return np.concatenate(blocks, axis=0)


def _span(residues: list[Res], chain: str) -> str:
    """'96-288 (193 residues)' for a chain -- the diagnostic a numbering mismatch
    needs, printed in the error rather than left for someone to go and look up."""
    nums = [
        int(re.sub(r"[A-Za-z]", "", r.resnum))
        for r in residues
        if r.chain == chain and re.sub(r"[A-Za-z]", "", r.resnum)
    ]
    if not nums:
        return "no residues"
    return f"{min(nums)}-{max(nums)} ({len(nums)} residues)"


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------


def contacts_per_residue(
    target: list[Res], design_xyz: np.ndarray, cutoff: float
) -> tuple[list[int], list[float | None]]:
    """Per target residue: how many (design atom, target atom) pairs are within the
    cutoff, and the closest approach. A residue with no heavy atoms gets 0 pairs and
    a min_dist of None (unmeasurable, not 'far away')."""
    counts: list[int] = []
    min_dists: list[float | None] = []
    for res in target:
        if len(res.heavy) == 0 or len(design_xyz) == 0:
            counts.append(0)
            min_dists.append(None)
            continue
        d = np.linalg.norm(res.heavy[:, None, :] - design_xyz[None, :, :], axis=-1)
        counts.append(int((d <= cutoff).sum()))
        min_dists.append(round(float(d.min()), 3))
    return counts, min_dists


def overlap(
    design: list[Res], target_xyz: np.ndarray, cutoff: float
) -> tuple[int | None, int | None]:
    """(n_clash, n_clash_design_res): heavy-atom PAIRS (design atom, target atom)
    closer than ``cutoff``, and the number of DISTINCT design residues taking part
    in at least one such pair. (None, None) when the target has no heavy atoms
    (nothing was measured -- never a zero).

    "Closer than" is decided on the distance ROUNDED TO 0.001 A, exactly as every
    distance this tool writes (the per-residue min_dist, min_dist_hotspot,
    min_dist_any) is rounded -- so n_clash == 0 if and only if min_dist_any >=
    cutoff, with no disagreement at the third decimal. The rounding is Python's
    round(), applied only to the few pairs within 0.001 A of the cutoff; every other
    pair is unambiguous on the raw distance.
    """
    if len(target_xyz) == 0:
        return None, None
    n_pairs = 0
    n_res = 0
    for res in design:
        if len(res.heavy) == 0:
            continue
        d = np.linalg.norm(res.heavy[:, None, :] - target_xyz[None, :, :], axis=-1)
        n = count_below(d, cutoff)
        n_pairs += n
        n_res += int(n > 0)
    return n_pairs, n_res


def count_below(d: np.ndarray, cutoff: float) -> int:
    """Distances whose value ROUNDED TO 0.001 A is strictly below ``cutoff``.
    Python's round() is applied only to the pairs within 0.001 A of the cutoff;
    every other pair is unambiguous on the raw distance."""
    lo, hi = cutoff - 0.001, cutoff + 0.001
    n = int((d < lo).sum())
    border = d[(d >= lo) & (d < hi)]
    return n + sum(1 for x in border if round(float(x), 3) < cutoff)


def backbone_atoms(
    model: gemmi.Model,
    chain_names: list[str],
    # A frozen COPY: an in-place change to BACKBONE_ATOMS cannot redefine n_clash_bb.
    atom_names: frozenset[str] = frozenset(BACKBONE_ATOMS),
) -> np.ndarray:
    """(M, 3) BACKBONE atoms (N, CA, C, O, OXT) of the amino-acid residues of the
    named chains, from the same model the rest of the measurement uses (altlocs and
    hydrogens already removed, waters skipped). A ligand or ion in a design chain
    has no backbone and contributes nothing here."""
    coords = []
    for name in chain_names:
        for chain in model:
            if chain.name != name:
                continue
            for res in chain:
                if _is_water(res) or not _is_amino_acid(res.name.strip()):
                    continue
                coords.extend(
                    [a.pos.x, a.pos.y, a.pos.z]
                    for a in res
                    if a.name.strip() in atom_names and _is_heavy(a)
                )
    return np.asarray(coords, dtype=float).reshape(-1, 3)


def backbone_overlap(
    bb_xyz: np.ndarray, target_xyz: np.ndarray, cutoff: float
) -> tuple[float | None, int | None]:
    """(min_dist_bb, n_clash_bb): design-backbone vs every target heavy atom, with
    the same rounding rules as min_dist_any / n_clash. (None, None) -- never 0 --
    when either side has no atoms."""
    if len(bb_xyz) == 0 or len(target_xyz) == 0:
        return None, None
    d = np.linalg.norm(bb_xyz[:, None, :] - target_xyz[None, :, :], axis=-1)
    return round(float(d.min()), 3), count_below(d, cutoff)


def design_identities(
    design: list[Res], sequence: str, seq_source: str, src_name: str
) -> tuple[dict[str, str], int, int]:
    """(label -> three-letter identity, seq_len, n_design_res_struct).

    With ``seq_source == 'structure'`` the identities are simply the structural
    residue names -- exactly what the tool did before ``--sequence-column`` existed.
    Otherwise ``sequence`` is mapped POSITIONALLY onto the design chains'
    amino-acid residues in structure order, and the mapping is checked rather than
    assumed. Non-amino-acid residues in a design chain (an ion, a ligand) are not
    part of the mapping and keep their structural name.
    """
    design_aa = [r for r in design if r.is_amino_acid]
    n_struct = len(design_aa)

    if seq_source == SEQ_SOURCE_STRUCTURE:
        return {r.label: r.resname for r in design}, n_struct, n_struct

    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError(
            f"--sequence-column {seq_source!r} is empty for this design, so the "
            f"design-chain identities cannot be read. Either that column does not "
            f"exist on this row or any ancestor, or the sequence step failed there"
        )
    if "/" in sequence:
        raise ValueError(
            f"--sequence-column {seq_source!r} holds a MULTI-CHAIN sequence for this "
            f"design ({sequence[:40]}...): it contains '/'. epitope maps the design "
            f"chain alone -- name the single designed chain's column (e.g. "
            f"proteinmpnn_sequence / atomium_sequence), not a complex sequence such "
            f"as mkcomplex_sequence"
        )
    if len(sequence) != n_struct:
        raise ValueError(
            f"sequence/structure length mismatch: seq_len={len(sequence)} from "
            f"column {seq_source!r} but n_design_res_struct={n_struct} amino-acid "
            f"residues in design chain(s) of {src_name}. Refusing to map -- a "
            f"positional mapping off by even one residue misreports the interface "
            f"composition of every row with no other symptom. Check that "
            f"{seq_source!r} is the single designed chain and that --design-chains "
            f"names the chain it belongs to"
        )

    names = {r.label: r.resname for r in design}  # non-AA keep the structural name
    for res, letter in zip(design_aa, sequence):
        names[res.label] = ONE_TO_THREE.get(letter, letter)
    return names, len(sequence), n_struct


def score(
    src: Path,
    design_chains: list[str],
    target_chains: list[str],
    hotspots: list[str],
    contact_cutoff: float,
    sequence: str = "",
    seq_source: str = SEQ_SOURCE_STRUCTURE,
    clash_cutoff: float = DEFAULT_CLASH_CUTOFF,
) -> tuple[str, dict, list[list]]:
    """One structure -> (status, metrics, per-target-residue rows). Raises on
    anything that makes the measurement impossible or dishonest; main() records the
    exception as the design's status."""
    model = read_model(src)
    present = chain_ids(model)

    missing_design = [c for c in design_chains if c not in present]
    if missing_design:
        raise ValueError(
            f"design chain(s) {','.join(missing_design)} absent from {src.name} "
            f"(have: {','.join(present)})"
        )
    missing_target = [c for c in target_chains if c not in present]
    if missing_target:
        raise ValueError(
            f"target chain(s) {','.join(missing_target)} absent from {src.name} "
            f"(have: {','.join(present)})"
        )

    design = residues_of(model, design_chains)
    design_xyz = heavy_atoms(design)
    if len(design_xyz) == 0:
        raise ValueError(
            f"design chain(s) {','.join(design_chains)} have no heavy atoms in "
            f"{src.name}"
        )

    target = residues_of(model, target_chains)
    if not target:
        raise ValueError(
            f"target chain(s) {','.join(target_chains)} have no non-water residues "
            f"in {src.name}"
        )

    # Hotspot resolution, BEFORE any number is computed: a hotspot that is not in
    # the target is an error, never a miss. (The chain check is also done at submit
    # time; repeated here because the worker is standalone-usable.)
    by_label: dict[str, Res] = {}
    ambiguous: list[str] = []
    for res in target:
        if res.label in by_label:
            ambiguous.append(res.label)
        else:
            by_label[res.label] = res

    stray = sorted({lab[0] for lab in hotspots} - set(target_chains))
    if stray:
        raise ValueError(
            f"hotspot chain(s) {','.join(stray)} are not target chain(s) "
            f"{','.join(target_chains)}"
        )
    unknown = [lab for lab in hotspots if lab not in by_label]
    if unknown:
        spans = "; ".join(f"{c}: {_span(target, c)}" for c in target_chains)
        raise ValueError(
            f"hotspot(s) {','.join(unknown)} do not exist in the target chain(s) of "
            f"{src.name} -- the numbering does not match. Target spans {spans}. "
            f"NOT scored as a miss: a renumbered chain would give a meaningless "
            f"hotspot_recall"
        )
    clashing = sorted(set(ambiguous) & set(hotspots))
    if clashing:
        raise ValueError(
            f"residue label(s) {','.join(clashing)} occur more than once in the "
            f"target chain(s) of {src.name}; a hotspot must name one residue"
        )

    # The design-chain identities, resolved (and length-checked) BEFORE any number
    # is computed, so a bad sequence column errors instead of producing plausible
    # resnames at the wrong positions.
    design_names, seq_len, n_design_res_struct = design_identities(
        design, sequence, seq_source, src.name
    )

    counts, min_dists = contacts_per_residue(target, design_xyz, contact_cutoff)
    by_index = {res.label: i for i, res in enumerate(target)}

    # The SAME matrix, read along its other axis: which design residues touch the
    # target. Reported in ascending residue order so resnums and resnames line up.
    target_xyz = heavy_atoms(target)
    design_counts, _ = contacts_per_residue(design, target_xyz, contact_cutoff)
    iface_design = sorted(
        (res for res, n in zip(design, design_counts) if n > 0),
        key=lambda r: r.sort_key,
    )

    iface = [res.label for res, n in zip(target, counts) if n > 0]
    hits = [lab for lab in hotspots if counts[by_index[lab]] > 0]
    hotspot_dists = [
        d
        for lab in hotspots
        if (d := min_dists[by_index[lab]]) is not None
    ]

    # WEIGHT, not just presence. The same `counts` vector, read on the listed
    # residues only: how many of the interface's contact pairs the hotspots carry.
    # No second distance computation and no second cutoff -- these are the numbers
    # already in `counts`, which is also what the per-residue TSV is written from.
    hotspot_counts = [counts[by_index[lab]] for lab in hotspots]
    hotspot_n_contacts = int(sum(hotspot_counts))
    n_contacts_total = int(sum(counts))
    # NA, never 0: with no interface at all the SHARE is undefined (0/0), while
    # hotspot_n_contacts is a real, measured 0.
    hotspot_contact_frac = (
        round(hotspot_n_contacts / n_contacts_total, 4)
        if n_contacts_total > 0
        else None
    )

    metrics = {
        "hotspot_recall": round(len(hits) / len(hotspots), 4),
        "hotspot_hits": ",".join(hits) if hits else NONE,
        "n_hotspots": len(hotspots),
        "hotspot_resnames": ",".join(
            f"{lab}={by_label[lab].resname}" for lab in hotspots
        ),
        "n_iface_target_res": len(iface),
        "iface_target_resnums": ",".join(iface) if iface else NONE,
        # NA only when every hotspot residue is atom-less -- not a distance of 0.
        "min_dist_hotspot": min(hotspot_dists) if hotspot_dists else None,
        "n_contacts": n_contacts_total,
        "design_chains": ",".join(design_chains),
        "target_chains": ",".join(target_chains),
        "contact_cutoff": contact_cutoff,
        # The design side of the matrix. 'none' (not a blank cell) for "measured,
        # and nothing touched", exactly as iface_target_resnums and hotspot_hits do
        # -- a blank would be read back as NA, i.e. "not measured".
        "n_iface_design_res": len(iface_design),
        "iface_design_resnums": (
            ",".join(r.label for r in iface_design) if iface_design else NONE
        ),
        "iface_design_resnames": (
            ",".join(f"{r.label}={design_names[r.label]}" for r in iface_design)
            if iface_design
            else NONE
        ),
        # Where the design identities came from. Read these before believing
        # iface_design_resnames.
        "seq_len": seq_len,
        "n_design_res_struct": n_design_res_struct,
        "seq_source": seq_source,
        # How much of the interface the LISTED hotspots carry. hotspot_recall says
        # 'touched or not'; these say 'how central'. hotspot_contacts is formatted
        # like hotspot_resnames and in the SAME order --hotspots was given in, so
        # the two can be read against each other residue by residue.
        "hotspot_n_contacts": hotspot_n_contacts,
        "hotspot_contact_frac": hotspot_contact_frac,
        "hotspot_contacts": ",".join(
            f"{lab}={n}" for lab, n in zip(hotspots, hotspot_counts)
        ),
    }

    # STERIC OVERLAP (added 2026-10-08). min_dist_any is the minimum of the very
    # per-residue min_dist values written to the per-residue file (each already
    # rounded to 0.001 A), so it equals that file's column minimum exactly and can
    # never exceed min_dist_hotspot (a minimum over a subset of the same values).
    measured = [d for d in min_dists if d is not None]
    n_clash, n_clash_design_res = overlap(design, target_xyz, clash_cutoff)
    metrics["min_dist_any"] = min(measured) if measured else None
    metrics["n_clash"] = n_clash
    metrics["n_clash_design_res"] = n_clash_design_res
    metrics["clash_cutoff"] = clash_cutoff
    # The design BACKBONE against every target heavy atom (second 2026-10-08 pass).
    # A subset of the pairs above, so min_dist_bb >= min_dist_any and n_clash_bb <=
    # n_clash; on a backbone-only design chain the pair sets are identical.
    bb_xyz = backbone_atoms(model, design_chains)
    if len(bb_xyz) == 0:
        print(
            f"{src.name}: design chain(s) {','.join(design_chains)} have no "
            f"backbone-named atoms (N/CA/C/O/OXT of amino-acid residues) -- "
            f"min_dist_bb and n_clash_bb are NA for this design"
        )
    metrics["min_dist_bb"], metrics["n_clash_bb"] = backbone_overlap(
        bb_xyz, target_xyz, clash_cutoff
    )
    # Backbone PLUS CB (third 2026-10-08 pass): a superset of the bb pairs and a
    # subset of all pairs, so the three readings bracket each other exactly.
    bbcb_xyz = backbone_atoms(model, design_chains, BACKBONE_CB_ATOMS)
    if len(bbcb_xyz) == 0:
        print(
            f"{src.name}: design chain(s) {','.join(design_chains)} have no "
            f"N/CA/C/O/OXT/CB atoms on amino-acid residues -- min_dist_bbcb and "
            f"n_clash_bbcb are NA for this design"
        )
    metrics["min_dist_bbcb"], metrics["n_clash_bbcb"] = backbone_overlap(
        bbcb_xyz, target_xyz, clash_cutoff
    )

    hotspot_set = set(hotspots)
    rows = [
        [
            res.chain,
            res.resnum,
            res.resname,
            n,
            "" if d is None else d,
            int(res.label in hotspot_set),
        ]
        for res, n, d in zip(target, counts, min_dists)
    ]
    return "OK", metrics, rows


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def write_result(out_dir: Path, name: str, status: str, path: str, metrics: dict):
    """One-row TSV. A metric that did not apply is written as NA (an empty cell),
    never as 0 -- a hotspot_recall of 0.0 is a real, informative value (the binder
    bound elsewhere) and must stay distinguishable from 'nothing was measured'."""
    write_tsv(
        out_dir / f"{name}.tsv",
        RESULT_COLUMNS,
        [
            [name, status, path]
            + [
                "" if metrics.get(c) is None else str(metrics[c])
                for c in METRIC_COLUMNS
            ]
        ],
    )


def positive_float(text: str) -> float:
    """argparse type: a finite float > 0."""
    value = float(text)
    if not np.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError(f"must be a positive number, got {text!r}")
    return value


class Args(argparse.Namespace):
    task_file: Path
    design_chains: str
    target_chains: str
    contact_cutoff: float
    clash_cutoff: float
    seq_source: str
    input_column: str
    out_dir: Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Per-design epitope contact map and hotspot recall."
    )
    ap.add_argument(
        "--task-file",
        type=Path,
        required=True,
        help="Tab-separated 'name<TAB>structure<TAB>hotspots<TAB>sequence' rows, one "
        "per design. The hotspot list is per design (already expanded and "
        "{expr}-resolved by the manifest builder), e.g. 'A96,A99,A101,A155'. The "
        "trailing sequence field is empty unless --seq-source names a column; a "
        "3-field row (the pre-2026-09-30 layout) is read as an empty sequence.",
    )
    ap.add_argument("--design-chains", required=True, help="Comma-joined chain IDs.")
    ap.add_argument("--target-chains", required=True, help="Comma-joined chain IDs.")
    ap.add_argument("--contact-cutoff", type=float, default=5.0)
    ap.add_argument(
        "--clash-cutoff",
        type=positive_float,
        default=DEFAULT_CLASH_CUTOFF,
        help="Heavy-atom distance (A) below which a design atom and a target atom "
        "overlap (strictly closer than, on distances rounded to 0.001 A). Default "
        f"{DEFAULT_CLASH_CUTOFF}. Feeds n_clash / n_clash_design_res / n_clash_bb / "
        "n_clash_bbcb only; every pre-existing column is independent of it.",
    )
    ap.add_argument(
        "--seq-source",
        default=SEQ_SOURCE_STRUCTURE,
        help="'structure' (default -- identities from the coordinates, the original "
        "behaviour), or the name of the table column the sequence field came from. "
        "Echoed into the seq_source column.",
    )
    ap.add_argument(
        "--input-column",
        default="the input structure column",
        help="Name of the table column the structure paths came from. Used only to "
        "make the 'could not resolve' error message name the right column.",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(namespace=Args())

    design_chains = split_list(args.design_chains)
    target_chains = split_list(args.target_chains)

    with open(args.task_file) as f:
        members = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    for member in members:
        if len(member) < 3:
            # Fail the whole task loudly: a task row without its hotspot list means
            # the manifest is wrong, and guessing a hotspot list would be worse than
            # stopping. The .sh wrapper then writes an error TSV per design.
            raise ValueError(
                f"malformed task row {member!r}: expected "
                f"'name<TAB>structure<TAB>hotspots'"
            )
        # The 4th field (the design sequence) may be absent or empty: absent is the
        # pre-2026-09-30 manifest layout, empty is "identities from the structure".
        name, src, hotspot_spec = member[0], member[1], member[2]
        sequence = member[3] if len(member) > 3 else ""
        per_res_path = args.out_dir / f"{name}_per_residue.tsv"
        try:
            if not src.strip():
                # The manifest builder could resolve no structure path for this
                # design, here or up its lineage. It submitted the row anyway, on
                # purpose: a design silently dropped from the manifest looks like a
                # smaller table, while an error status is a fact in the table.
                raise ValueError(
                    f"could not resolve a structure path from column "
                    f"{args.input_column!r} for this design -- neither on its own "
                    f"row nor on any ancestor (parent_name/parent_table). Check "
                    f"-i/--input-column, and that the parent generation actually "
                    f"produced that column"
                )
            hotspots = split_list(hotspot_spec)
            if not hotspots:
                raise ValueError("no hotspots given for this design")
            bad = [lab for lab in hotspots if not _LABEL.match(lab)]
            if bad:
                raise ValueError(
                    f"malformed hotspot label(s) {','.join(bad)} (expected 'A96')"
                )
            status, metrics, rows = score(
                Path(src),
                design_chains,
                target_chains,
                hotspots,
                args.contact_cutoff,
                sequence,
                args.seq_source,
                args.clash_cutoff,
            )
            write_tsv(per_res_path, PER_RESIDUE_COLUMNS, rows)
            write_result(args.out_dir, name, status, str(per_res_path), metrics)
            print(
                f"{name}: recall={metrics['hotspot_recall']} "
                f"hits={metrics['hotspot_hits']} "
                f"n_iface={metrics['n_iface_target_res']} "
                f"n_iface_design={metrics['n_iface_design_res']} "
                f"min_dist_hotspot={metrics['min_dist_hotspot']} "
                f"({metrics['seq_source']}) {status}"
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            print(f"{name}: ERROR {e}", file=sys.stderr)
            per_res_path.unlink(missing_ok=True)
            write_result(args.out_dir, name, f"error: {e}", "", {})


if __name__ == "__main__":
    main()
