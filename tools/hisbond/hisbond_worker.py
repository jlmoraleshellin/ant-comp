#!/usr/bin/env python3
"""
Does a designed histidine hydrogen-bond the TARGET? The per-array-task step of the
hisbond tool. Interpretation, traps and how to read the columns:
``.claude/skills/hisbond/SKILL.md``.

PREMISE. One predicted complex per design, holding a binder chain (``--binder-chain``),
target chains (``--target-chains``) and a list of centre POSITIONS. A position is
1-BASED ALONG THE BINDER CHAIN IN FILE ORDER (position 1 is the chain's first
amino-acid residue), NOT an auth residue number and not an index. Every index in this
tool's output -- the ``per_centre*`` cell prefixes included -- is that same 1-based
position, so ``centre<k>_idx1`` IS the join key into every cell.

The positions come from the caller, as plain integers; this tool knows no other tool's
column names. Per-design values reach it through prosapia's ``{expr}`` islands, e.g.
``--centers '{atomiumph_centre1_idx1},{atomiumph_centre2_idx1}'``, resolved at
manifest-build time. THE RESIDUE AT EACH POSITION MUST BE A HISTIDINE or the design
errors: that identity check against the structure is the whole reason positions are
given rather than histidines auto-discovered, and it is what makes a table of Asp/Glu
centres (an atomiumph acid_p / all_p run) error on every row.

THE MEASUREMENT IS A HEAVY-ATOM PROXY. Predictions carry no hydrogens, so a contact
between a ring N (ND1 or NE2) and a partner N/O atom qualifies when::

    --min-dist <= dist <= --max-dist                      (defaults 2.4 and 3.5 A)
    angle(external ring bisector at N, N->partner) <= --max-angle  (default 60 deg)

The lower bound is load-bearing: a partner under --min-dist is a clash, excluded from
both counts and reported as n_too_close_*, while ``min_dist_target`` stays the
UNCONDITIONAL minimum. Hydrogens are dropped from the contact geometry even when the
file has them, so counts stay comparable between a prediction and a relaxed pose; they
are read only for the tautomer. ND1 and NE2 are also reported separately, from the same
per-nitrogen masks the merged (unique-atom) counts are the OR of.

``centre<k>_idx1`` is the k-th supplied position, read back off the residue the
measurement actually ran on (the numbering atomium/ProteinMPNN use, every chain
renumbered 1..L), and ``centre<k>_resnum`` that residue's auth_seq_id. ``_ordinal_in``
returns the supplied position by construction, so what makes idx1 trustworthy is the
His-name check and ``--expect-binder-len-col``; pass it whenever the positions are
used. ``centers_spec`` carries every supplied position, past the column cap too.

ROLE OF EACH CONTACT. The tautomer is read from hydrogens when the file has them (an
H/D of the centre's OWN residue within TAUTOMER_H_DIST of a ring N, by distance not by
name): the N carrying the H donates, the other accepts. Otherwise the role comes from
the partner (PARTNER_CLASSES). A contact is ACCEPTOR-ROLE, DONOR-ROLE or INCOHERENT (a
trust signal, tautomer mode only); ``role_basis`` says which evidence was used. The
acceptor role is the pH-sensitive one (acid protonates the lone-pair N), which assumes
bind at neutral and release in acid. ``n_centres_switch_candidate_clean`` is what a
shortlist gates on.

The same roles are read for INTRA-BINDER contacts (``*_binder``), over the contacts
``n_centres_hb_binder`` counts; the two sides are disjoint by construction. A designed
His as partner reads like any His, which can inflate the candidate count, never hide
one ('*' marks it). All of it is a heavy-atom proxy: a role is not a pKa, and in
partner_class mode ``n_contacts_incoherent_binder == 0`` is a tautology.

``centre<k>_n_acidic_env`` / ``_n_basic_env`` count charged ATOMS of the TARGET chains
within ``--env-radius`` of the ring centroid (pKa context, not a pKa).

THE MAPPING IS CHECKED. A centre that is not on a histidine fails the design, and
``--expect-binder-len-col`` checks the binder's residue count against a table column.
NA IS NEVER 0: an unmeasurable design gets ``status = 'error: ...'`` and empty columns.

Normally run per array task by hisbond.sh; run ``sapia run hisbond`` rather than this.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import NamedTuple

import gemmi
import numpy as np

# Residue names that ARE a histidine; the protonation variants share the ring geometry.
HIS_NAMES = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}

# Ring N -> neighbours; insertion order (ND1 first) is the per-nitrogen reporting order.
RING_NEIGHBOURS = {"ND1": ("CG", "CE1"), "NE2": ("CD2", "CE1")}

# Lower-case column suffix per nitrogen (..._nd1, not ..._ND1).
NITROGEN_SUFFIX = {"ND1": "nd1", "NE2": "ne2"}

# Fixed centre<k>_* column set; extra centres are measured but warned about.
# 8, raised from 4: a his_s/his_p atomiumph run pins --centers-per-type centres (default
# 2) and a caller may raise that, while acid_p/all_p tables error here on residue
# identity anyway. Each centre costs FIVE columns (idx1, resnum, tautomer, n_acidic_env,
# n_basic_env), so 8 is already 40 columns on every table hisbond touches. Past the cap
# nothing is lost: centers_spec lists every position, every n_centres_*/per_centre*
# column covers every centre, and centres_reported < centers_n makes the cap visible.
# KEEP IN SYNC BY HAND with collect_hisbond.py (no shared imports).
MAX_CENTRE_COLUMNS = 8

# Polar heavy partner elements; S is deliberately absent.
POLAR_ELEMENTS = {"N", "O"}

# The five ring atoms whose mean is the centroid the pH-environment counts use.
RING_ATOMS = ("CG", "ND1", "CD2", "CE1", "NE2")

# N-H is ~1.0 A, the next H >= 2 A: 1.2 admits only the bonded H. Read BY DISTANCE.
TAUTOMER_H_DIST = 1.2

# Conventional ring-N H names, used only to cross-check the distance reading (warn).
EXPECTED_H_NAMES = {"ND1": {"HD1", "DD1"}, "NE2": {"HE2", "DE2"}}

# Standard-rules donor/acceptor table (no geometry); absent atoms are UNCLASSIFIED.
DONOR_ONLY = "donor"
ACCEPTOR_ONLY = "acceptor"
BOTH_CAPABLE = "both"
UNCLASSIFIED = "unclass"

PARTNER_CLASSES = {
    # Backbone (any residue). PRO N is special-cased in _partner_class.
    ("*", "N"): DONOR_ONLY,
    ("*", "O"): ACCEPTOR_ONLY,
    ("*", "OXT"): ACCEPTOR_ONLY,
    # Side-chain donors only.
    ("LYS", "NZ"): DONOR_ONLY,
    ("ARG", "NE"): DONOR_ONLY,
    ("ARG", "NH1"): DONOR_ONLY,
    ("ARG", "NH2"): DONOR_ONLY,
    ("ASN", "ND2"): DONOR_ONLY,
    ("GLN", "NE2"): DONOR_ONLY,
    ("TRP", "NE1"): DONOR_ONLY,
    # Side-chain acceptors only.
    ("ASP", "OD1"): ACCEPTOR_ONLY,
    ("ASP", "OD2"): ACCEPTOR_ONLY,
    ("GLU", "OE1"): ACCEPTOR_ONLY,
    ("GLU", "OE2"): ACCEPTOR_ONLY,
    ("ASN", "OD1"): ACCEPTOR_ONLY,
    ("GLN", "OE1"): ACCEPTOR_ONLY,
    # Hydroxyls and a target histidine: can do either.
    ("SER", "OG"): BOTH_CAPABLE,
    ("THR", "OG1"): BOTH_CAPABLE,
    ("TYR", "OH"): BOTH_CAPABLE,
    ("HIS", "ND1"): BOTH_CAPABLE,
    ("HIS", "NE2"): BOTH_CAPABLE,
    # A DOUBLY PROTONATED histidine partner is DONOR-ONLY at BOTH ring nitrogens: it
    # carries an H on each, so it has no lone pair to offer. Named explicitly because
    # HIP/HSP is the one variant name that is evidence rather than a spelling, and
    # folding it into HIS would score `don` where `incoh` is owed -- hiding exactly the
    # incoherence n_centres_switch_candidate_clean gates on (see _canonical_partner).
    ("HIP", "ND1"): DONOR_ONLY,
    ("HIP", "NE2"): DONOR_ONLY,
    ("HSP", "ND1"): DONOR_ONLY,
    ("HSP", "NE2"): DONOR_ONLY,
}

# His variant names that are only a SPELLING of neutral His: the tautomer they encode
# says which N carries the H, which this table does not model per nitrogen, so both
# nitrogens stay both-capable. HIP/HSP are deliberately ABSENT -- they are evidence.
HIS_SPELLING_ALIASES = frozenset({"HID", "HIE", "HSD", "HSE"})

# Charged side-chain atoms of the pH environment: atoms, not residues (Asp 2, Arg 3).
ENV_ACIDIC_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
ENV_BASIC_ATOMS = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")}

# Ring-centroid radius (A) for those counts: context only, echoed as env_radius.
DEFAULT_ENV_RADIUS = 6.0

# Backbone atom names, excluded from the side-chain pLDDT (B-factor) minimum.
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

# Set by run_hisbond.py in the sub-manifest; kept in sync by hand (no shared imports).
UNRESOLVED_EXPECT = "__unresolved__"

# A max B-factor at or below this is on the 0-1 scale, not 0-100 (as in boltzconf).
PLDDT_UNIT_MAX = 1.0

# Bare metric names; collect_hisbond.py hands them to the driver (trust columns first).
TRUST_COLUMNS = [
    "binder_chain",
    "target_chains",
    # The --centers template verbatim (provenance) and this design's resolved
    # positions; centers_spec covers every centre, past MAX_CENTRE_COLUMNS too.
    "centers_arg",
    "centers_spec",
    # 0 = every island came from the design's own table; otherwise the DEEPEST one.
    "centers_lineage_depth",
    "centers_n",
    "centers_his_ok",
    "binder_chain_n_res",
    "expect_binder_len",
    "target_chain_n_res",
    "target_n_het_atoms",
    "chain_ids_found",
    "structure_has_hydrogens",
    "plddt_scale",
    # 'tautomer', 'partner_class' or 'mixed'; read before any role column (NA on error).
    "role_basis",
    "min_dist",
    "max_dist",
    "max_angle",
    "env_radius",
]
SCIENCE_COLUMNS = [
    "per_centre",
    "n_centres_hb_target",
    "any_target_hb_proxy",
    "min_dist_target",
    # Per-nitrogen split of the same masks the merged columns are the OR of.
    "n_centres_hb_target_nd1",
    "n_centres_hb_target_ne2",
    "any_target_hb_proxy_nd1",
    "any_target_hb_proxy_ne2",
    "min_dist_target_nd1",
    "min_dist_target_ne2",
    "per_centre_nd1_target",
    "per_centre_ne2_target",
    # What each qualifying contact IS (same entries as the cells above, no new contact).
    "per_centre_role_target",
    "n_centres_acceptor_role_target",
    "n_centres_donor_role_target",
    "n_centres_switch_candidate",
    "n_centres_switch_candidate_clean",
    "any_switch_candidate",
    "n_contacts_incoherent_target",
    "n_contacts_ambiguous_target",
    "n_contacts_unclassified_target",
    # Row totals of the three reading defects (all centres, expected 0).
    "n_centres_no_ring_h",
    "n_centres_hip",
    "n_centres_foreign_h",
    # Centre positions for pinning, read off the residue actually found.
    "centres_reported",
    "centres_idx1_eq_resnum",
    *[
        f"centre{k}_{field}"
        for k in range(1, MAX_CENTRE_COLUMNS + 1)
        for field in ("idx1", "resnum")
    ],
    # Tautomer and pH environment of each reported centre (same cap and order).
    *[f"centre{k}_tautomer" for k in range(1, MAX_CENTRE_COLUMNS + 1)],
    *[
        f"centre{k}_{field}"
        for k in range(1, MAX_CENTRE_COLUMNS + 1)
        for field in ("n_acidic_env", "n_basic_env")
    ],
    "n_centres_hb_binder",
    "min_dist_binder",
    # The intra-binder reading of those contacts, mirroring the target-side role block.
    "per_centre_role_binder",
    "n_centres_acceptor_role_binder",
    "n_centres_donor_role_binder",
    "n_centres_switch_candidate_binder",
    "n_centres_switch_candidate_clean_binder",
    "any_switch_candidate_binder",
    "n_contacts_incoherent_binder",
    "n_contacts_ambiguous_binder",
    "n_contacts_unclassified_binder",
    "n_contacts_binder_entries",
    "n_centres_switch_candidate_clean_any",
    "any_switch_candidate_clean_any",
    "per_centre_nd1_binder",
    "per_centre_ne2_binder",
    "n_centres_polar_dist_only_target",
    "n_centres_too_close_target",
    "n_centres_too_close_binder",
    "his_plddt_min",
]
METRIC_COLUMNS = TRUST_COLUMNS + SCIENCE_COLUMNS
RESULT_COLUMNS = ["name", "status", "path", *METRIC_COLUMNS]

# Echo of the question rather than a measurement: kept on an error row.
ECHO_COLUMNS = {
    "binder_chain",
    "target_chains",
    "centers_arg",
    "min_dist",
    "max_dist",
    "max_angle",
    "env_radius",
}

_CIF_SUFFIXES = (".cif", ".mmcif", ".cif.gz", ".mmcif.gz")


class Args(argparse.Namespace):
    task_file: Path
    binder_chain: str
    target_chains: str
    min_dist: float
    max_dist: float
    max_angle: float
    env_radius: float
    input_column: str
    centers_arg: str
    out_dir: Path


class StructureError(ValueError):
    """The file was read but holds nothing measurable, or could not be parsed."""


def split_list(value: str) -> list[str]:
    return [tok.strip() for tok in value.split(",") if tok.strip()]


# --- atom records (not gemmi's reader: it merges repeated residue numbers) ---


class Residue:
    """One residue, as the file lists it: name, key, and its atoms."""

    __slots__ = ("chain", "resname", "seq", "icode", "names", "elements", "xyz", "bfac")

    def __init__(self, chain: str, resname: str, seq: str, icode: str) -> None:
        self.chain = chain
        self.resname = resname
        self.seq = seq
        self.icode = icode
        self.names: list[str] = []
        self.elements: list[str] = []
        self.xyz: list[tuple[float, float, float]] = []
        self.bfac: list[float | None] = []

    def add(
        self,
        name: str,
        element: str,
        xyz: tuple[float, float, float],
        bfac: float | None,
    ) -> None:
        """Add an atom, keeping the FIRST conformer of each atom name.

        A repeat is an altloc copy; first-seen also works for conformers labelled
        'B'/'C', where an altloc-letter rule would drop the residue's atoms entirely.
        """
        if name in self.names:
            return
        self.names.append(name)
        self.elements.append(element)
        self.xyz.append(xyz)
        self.bfac.append(bfac)

    def atom(self, name: str) -> tuple[float, float, float] | None:
        """First atom with this exact name, or None."""
        for i, n in enumerate(self.names):
            if n == name:
                return self.xyz[i]
        return None

    def label(self) -> str:
        return f"{self.chain}{self.seq}{self.icode}"


_TAB_CACHE: dict[str, bool] = {}


def _excluded_residue(resname: str) -> bool:
    """True for a residue gemmi tabulates as something other than an amino acid.

    Untabulated names are not excluded: an unusual designed residue keeps its index.
    """
    key = resname.strip().upper()
    if key not in _TAB_CACHE:
        info = gemmi.find_tabulated_residue(key)
        _TAB_CACHE[key] = bool(info) and not info.is_amino_acid()
    return _TAB_CACHE[key]


def _element_of(name: str, element: str) -> str:
    """The upper-cased element symbol; guessed from the atom name if blank."""
    if element:
        return element.upper()
    stripped = name.strip()
    # Blank element column: guess the first letter (' CA ' is C, 'CA  ' is calcium).
    for ch in stripped:
        if ch.isalpha():
            return ch.upper()
    return ""


class _Collector:
    """Accumulates residues per chain, in file order."""

    def __init__(self) -> None:
        self.chains_present: list[str] = []
        self.residues: dict[str, list[Residue]] = {}
        self._last_key: dict[str, tuple[str, str, str]] = {}
        self.has_hydrogens = False
        # H/D atoms stay out of the residues; used only for the tautomer reading.
        self.h_xyz: list[tuple[float, float, float]] = []
        self.h_names: list[str] = []
        # Owner of each H/D: only the centre's own hydrogens may set its tautomer.
        self.h_keys: list[tuple[str, str, str]] = []
        self.max_bfac: float | None = None

    def saw_chain(self, chain: str) -> None:
        if chain not in self.residues:
            self.chains_present.append(chain)
            self.residues[chain] = []

    def add_atom(
        self,
        chain: str,
        resname: str,
        seq: str,
        icode: str,
        altloc: str,
        name: str,
        element: str,
        xyz: tuple[float, float, float],
        bfac: float | None,
    ) -> None:
        self.saw_chain(chain)
        el = _element_of(name, element)
        if el in ("H", "D"):
            # Kept aside: contacts are heavy-atom only, H/D serve the tautomer reading.
            self.has_hydrogens = True
            self.h_xyz.append(xyz)
            self.h_names.append(name.strip().upper())
            # Keyed by (chain, seq, icode); a water's H is therefore always FOREIGN.
            self.h_keys.append((chain, seq, icode))
            return
        # Residue.add keeps the first conformer per atom name, so altloc is unused.
        key = (resname, seq, icode)
        if self._last_key.get(chain) != key:
            self._last_key[chain] = key
            self.residues[chain].append(Residue(chain, resname, seq, icode))
        if bfac is not None:
            self.max_bfac = bfac if self.max_bfac is None else max(self.max_bfac, bfac)
        self.residues[chain][-1].add(name, el, xyz, bfac)


def _is_cif(path: Path) -> bool:
    return path.name.lower().endswith(_CIF_SUFFIXES)


def _read_text(path: Path) -> str:
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return f.read()
    return path.read_text()


def _read_pdb(path: Path) -> _Collector:
    out = _Collector()
    n_atoms = 0
    for line in _read_text(path).splitlines():
        record = line[:6]
        if record == "ENDMDL" and n_atoms:
            break  # first model only
        if record not in ("ATOM  ", "HETATM"):
            continue
        if len(line) < 54:
            raise StructureError(f"truncated {record.strip()} record: {line!r}")
        n_atoms += 1
        bfac: float | None = None
        if len(line) >= 66:
            try:
                bfac = float(line[60:66])
            except ValueError:
                bfac = None
        out.add_atom(
            chain=line[20:22].strip(),
            resname=line[17:21].strip(),
            seq=line[22:26].strip(),
            icode=line[26].strip(),
            altloc=line[16].strip(),
            name=line[12:16].strip(),
            element=line[76:78].strip() if len(line) >= 78 else "",
            xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
            bfac=bfac,
        )
    if n_atoms == 0:
        raise StructureError(f"no ATOM/HETATM records in {path.name}")
    return out


def _read_cif(path: Path) -> _Collector:
    block = gemmi.cif.read(str(path)).sole_block()

    def col(*tags: str) -> list[str] | None:
        for tag in tags:
            values = block.find_values(f"_atom_site.{tag}")
            if len(values):
                return [gemmi.cif.as_string(v) for v in values]
        return None

    x, y, z = col("Cartn_x"), col("Cartn_y"), col("Cartn_z")
    chain_col = col("auth_asym_id", "label_asym_id")
    atom_col = col("auth_atom_id", "label_atom_id")
    comp_col = col("auth_comp_id", "label_comp_id")
    if (
        x is None
        or y is None
        or z is None
        or chain_col is None
        or atom_col is None
        or comp_col is None
    ):
        raise StructureError(f"no usable _atom_site loop in {path.name}")
    seq_col = col("auth_seq_id", "label_seq_id")
    icode_col = col("pdbx_PDB_ins_code")
    alt_col = col("label_alt_id")
    elem_col = col("type_symbol")
    bfac_col = col("B_iso_or_equiv")
    model_col = col("pdbx_PDB_model_num")

    def blank(v: str) -> str:
        return "" if v in ("", ".", "?") else v

    out = _Collector()
    first_model = model_col[0] if model_col else None
    for i in range(len(x)):
        if model_col is not None and model_col[i] != first_model:
            continue  # first model only
        bfac: float | None = None
        if bfac_col is not None and blank(bfac_col[i]):
            try:
                bfac = float(bfac_col[i])
            except ValueError:
                bfac = None
        out.add_atom(
            chain=chain_col[i],
            resname=comp_col[i],
            seq=blank(seq_col[i]) if seq_col else "",
            icode=blank(icode_col[i]) if icode_col else "",
            altloc=blank(alt_col[i]) if alt_col else "",
            name=atom_col[i],
            element=blank(elem_col[i]) if elem_col else "",
            xyz=(float(x[i]), float(y[i]), float(z[i])),
            bfac=bfac,
        )
    if not out.chains_present:
        raise StructureError(f"no atoms in the _atom_site loop of {path.name}")
    return out


def read_structure(path: Path) -> _Collector:
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    try:
        return _read_cif(path) if _is_cif(path) else _read_pdb(path)
    except StructureError:
        raise
    except Exception as e:  # noqa: BLE001 - any parse failure is "unreadable"
        raise StructureError(f"unreadable structure {path.name}: {e}") from e


# --- centres ---


def parse_centres(text: str) -> list[int]:
    """``'12,45'`` -> ``[12, 45]``: 1-BASED positions along the binder chain.

    Order is preserved and IS the order of ``centre<k>_*`` and of every ``;``-joined
    ``per_centre*`` cell. Raises on anything malformed: the position list defines the
    question. Kept in sync BY HAND with ``parse_centres`` in ``run_hisbond.py`` (no
    shared imports); the builder validates before submitting, this validates what
    actually reached the task file.
    """
    positions: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            raise ValueError(
                f"empty position in centre list {text!r}; expected comma-separated "
                f"1-based positions, e.g. '12,45'"
            )
        try:
            value = int(part)
        except ValueError:
            raise ValueError(
                f"centre position {part!r} in {text!r} is not an integer"
            ) from None
        if value < 1:
            raise ValueError(
                f"centre position {value} in {text!r} is not >= 1; positions are "
                f"1-BASED along the binder chain, so 0 is never one of them"
            )
        positions.append(value)
    if not positions:
        raise ValueError("the centre list is empty: there is nothing to measure")
    if len(set(positions)) != len(positions):
        dup = sorted({p for p in positions if positions.count(p) > 1})
        raise ValueError(
            f"centre position(s) {dup} appear twice in {text!r}; one residue measured "
            f"twice would be counted twice in every n_centres_* column"
        )
    return positions


# --- the measurement ---


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n <= 0:
        raise ValueError("zero-length vector")
    return v / n


def bisector(n_xyz: np.ndarray, nb1: np.ndarray, nb2: np.ndarray) -> np.ndarray:
    """The EXTERNAL bisector of the ring angle at an imidazole N.

    ``-(u1 + u2)`` normalised, ``u1``/``u2`` the unit vectors to its ring neighbours:
    where that N's H or lone pair points. A degenerate sum (collinear ring) raises."""
    s = _unit(nb1 - n_xyz) + _unit(nb2 - n_xyz)
    if float(np.linalg.norm(s)) < 1e-6:
        raise ValueError("imidazole ring neighbours are collinear with the N")
    return _unit(-s)


def _partner_class(resname: str, atom: str) -> str:
    """DONOR_ONLY, ACCEPTOR_ONLY, BOTH_CAPABLE or UNCLASSIFIED for one partner atom.

    A lookup in PARTNER_CLASSES, not a calculation: an exact (resname, atom) entry beats
    the ('*', atom) backbone rule, PRO's backbone N is UNCLASSIFIED, and anything the
    table does not name is counted as UNCLASSIFIED rather than given a role.

    NEUTRAL HISTIDINE SPELLINGS ARE CANONICALISED TO 'HIS'; HIP/HSP ARE NOT.
    HIS_NAMES accepts HID/HIE/HIP/HSD/HSE/HSP for a CENTRE, but PARTNER_CLASSES keys
    only ('HIS', ND1/NE2); without canonicalising, a His PARTNER written under a variant
    name (what a CHARMM/AMBER-style relaxer emits) falls through to the ('*', name)
    backbone rule, which has no ND1/NE2 entry, and scores UNCLASSIFIED instead of
    BOTH_CAPABLE. That is exactly the designed-His-to-designed-His intra-binder contact
    the '*' mark exists to flag, so it would subtract from `acc` and make
    n_centres_switch_candidate_clean_binder and _clean_any UNDER-report.

    ONLY HIS_SPELLING_ALIASES (HID/HIE/HSD/HSE) are folded in. HIP/HSP name a DOUBLY
    PROTONATED ring and are the only protonation evidence the file offers about a
    PARTNER, so they keep their own DONOR_ONLY entries. Folding them into HIS would be
    one-sided and in the dangerous direction: at a centre's protonated nitrogen a
    both-capable partner scores `don` where donor-only scores `incoh` (_contact_roles),
    so a HIP partner would HIDE an incoherence -- and `incoh == 0` is exactly what
    n_centres_switch_candidate_clean and _clean_any gate on. At the lone-pair nitrogen
    both classes give `acc`, so the error would be FALSE-CLEAN, never false-dirty."""
    res = resname.strip().upper()
    name = atom.strip().upper()
    if res in HIS_SPELLING_ALIASES:
        res = "HIS"
    if name == "N" and res == "PRO":
        return UNCLASSIFIED
    exact = PARTNER_CLASSES.get((res, name))
    if exact is not None:
        return exact
    return PARTNER_CLASSES.get(("*", name), UNCLASSIFIED)


def _polar_atoms(residues: list[Residue]) -> tuple[np.ndarray, list[str], list[str]]:
    """Coordinates, labels and partner classes of every N/O heavy atom in ``residues``.

    A label is ``<chain>:<resnum><icode>:<resname>:<atom>`` in the file's own auth
    namespace (what a viewer shows), not a 0-based index; a missing residue number is
    '?'. The class is computed from the same (resname, atom) the label prints, so a role
    and its label cannot drift apart."""
    xyz: list[tuple[float, float, float]] = []
    labels: list[str] = []
    classes: list[str] = []
    for res in residues:
        for i, el in enumerate(res.elements):
            if el in POLAR_ELEMENTS:
                xyz.append(res.xyz[i])
                seq = f"{res.seq}{res.icode}" or "?"
                labels.append(f"{res.chain}:{seq}:{res.resname}:{res.names[i]}")
                classes.append(_partner_class(res.resname, res.names[i]))
    return np.asarray(xyz, dtype=float).reshape(-1, 3), labels, classes


# Marks a binder partner that is itself a designed centre (by list position).
DESIGNED_PARTNER_MARK = "*"


def _binder_partner_atoms(
    binder_res: list[Residue], idx: int, designed: frozenset[int] = frozenset()
) -> tuple[np.ndarray, list[str], list[str]]:
    """N/O atoms of the binder chain available to the centre at offset ``idx``.

    ``idx`` and ``designed`` are 0-BASED OFFSETS into ``binder_res``, not the 1-based
    positions the columns report: they only index this list.

    Returns xyz, labels and partner classes. All atoms of the centre residue are
    excluded, and only the BACKBONE of i-1 / i+1 (at H-bond distance by local geometry);
    adjacency is by list position. A label whose residue is itself a designed centre
    (``designed``) ends in DESIGNED_PARTNER_MARK; the class ignores the mark."""
    xyz: list[tuple[float, float, float]] = []
    labels: list[str] = []
    classes: list[str] = []
    for j, res in enumerate(binder_res):
        if j == idx:
            continue
        drop_backbone = j in (idx - 1, idx + 1)
        for i, el in enumerate(res.elements):
            if el not in POLAR_ELEMENTS:
                continue
            if drop_backbone and res.names[i] in BACKBONE_ATOMS:
                continue
            xyz.append(res.xyz[i])
            seq = f"{res.seq}{res.icode}" or "?"
            mark = DESIGNED_PARTNER_MARK if j in designed else ""
            labels.append(f"{res.chain}:{seq}:{res.resname}:{res.names[i]}{mark}")
            classes.append(_partner_class(res.resname, res.names[i]))
    return np.asarray(xyz, dtype=float).reshape(-1, 3), labels, classes


class NContacts(NamedTuple):
    """One imidazole nitrogen against one side's partner atoms, atom by atom.

    ``near``/``hb``/``too_close`` are masks over the partner array, ``dist``/``angle``
    the per-atom distance (A) and angle to the ring bisector (deg), ``min_d`` the
    UNCONDITIONAL minimum distance (None with no partner atoms: NA, not 0)."""

    near: np.ndarray
    hb: np.ndarray
    too_close: np.ndarray
    dist: np.ndarray
    angle: np.ndarray
    min_d: float | None


def _contacts(
    n_xyz: np.ndarray,
    direction: np.ndarray,
    partners: np.ndarray,
    min_dist: float,
    max_dist: float,
    max_angle: float,
) -> NContacts:
    """Per partner atom: inside the distance window, inside it AND the cone, too close.

    The floor is part of the contact test (a sub-minimum atom never reaches a count) and
    ``too_close`` makes the exclusion visible. ``min_d`` is unconditional, None if
    empty."""
    if len(partners) == 0:
        z = np.zeros(0, dtype=bool)
        f = np.zeros(0, dtype=float)
        return NContacts(z, z.copy(), z.copy(), f, f.copy(), None)
    v = partners - n_xyz
    d = np.linalg.norm(v, axis=1)
    too_close = d < min_dist
    near = (d >= min_dist) & (d <= max_dist)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos = (v @ direction) / d
    cos = np.clip(cos, -1.0, 1.0)
    ang = np.degrees(np.arccos(cos))
    return NContacts(near, near & (ang <= max_angle), too_close, d, ang, float(d.min()))


def _centre_nitrogens(
    his: Residue,
    partners: np.ndarray,
    min_dist: float,
    max_dist: float,
    max_angle: float,
) -> dict[str, NContacts]:
    """One centre against one side, PER NITROGEN: the single source of both views.

    The merged counts are the OR of these masks (``_merge_nitrogens``) and the
    per-nitrogen columns count them one at a time. Raises on an incomplete imidazole.
    """
    out: dict[str, NContacts] = {}
    for n_name, (nb1_name, nb2_name) in RING_NEIGHBOURS.items():
        n_pos = his.atom(n_name)
        nb1 = his.atom(nb1_name)
        nb2 = his.atom(nb2_name)
        if n_pos is None:
            raise ValueError(
                f"histidine {his.label()} has no {n_name} atom; the imidazole is "
                f"incomplete and no H-bond can be computed"
            )
        if nb1 is None or nb2 is None:
            missing = ",".join(
                nm for nm, p in ((nb1_name, nb1), (nb2_name, nb2)) if p is None
            )
            raise ValueError(
                f"histidine {his.label()} has no {missing} atom; the orientation of "
                f"{n_name} cannot be computed"
            )
        n_arr = np.asarray(n_pos, dtype=float)
        direction = bisector(
            n_arr, np.asarray(nb1, dtype=float), np.asarray(nb2, dtype=float)
        )
        out[n_name] = _contacts(
            n_arr, direction, partners, min_dist, max_dist, max_angle
        )
    return out


def _merge_nitrogens(
    per_n: dict[str, NContacts], n_partners: int
) -> tuple[int, int, int, float | None]:
    """The ND1/NE2 views of one centre+side merged.

    Returns (n_partner_atoms, n_near, n_too_close, min_distance): counts of UNIQUE
    partner atoms qualifying with either nitrogen (an atom seen by both is one contact
    here and appears in both per-nitrogen cells); atoms, not residues."""
    near_any = np.zeros(n_partners, dtype=bool)
    hb_any = np.zeros(n_partners, dtype=bool)
    close_any = np.zeros(n_partners, dtype=bool)
    min_d: float | None = None
    for c in per_n.values():
        near_any |= c.near
        hb_any |= c.hb
        close_any |= c.too_close
        if c.min_d is not None:
            min_d = c.min_d if min_d is None else min(min_d, c.min_d)
    return int(hb_any.sum()), int(near_any.sum()), int(close_any.sum()), min_d


def _partner_entries(idx1: int, contacts: NContacts, labels: list[str]) -> str:
    """The qualifying partner atoms of ONE centre through ONE nitrogen, as text.

    ``<idx1>:<chain>:<resnum>:<resname>:<atom>=<dist>/<angle>`` per atom, joined by
    ';'; ``<idx1>:none`` when nothing qualified, so that stays distinct from an absent
    centre. ``<idx1>`` is the 1-BASED binder position, i.e. ``centre<k>_idx1`` itself:
    that is the join key between these cells and the centre columns."""
    hits = [
        f"{idx1}:{labels[i]}={contacts.dist[i]:.2f}/{contacts.angle[i]:.1f}"
        for i in np.flatnonzero(contacts.hb)
    ]
    return ";".join(hits) if hits else f"{idx1}:none"


# Order of the five role counts in per_centre_role_*; they partition a centre's entries.
ROLE_FIELDS = ("acc", "don", "amb", "incoh", "unclass")

# Looser probe only to warn when an H's name and its geometry disagree.
NAME_CHECK_DIST = 2.0


def _h_names_within(
    point: tuple[float, float, float],
    h_xyz: np.ndarray,
    h_names: list[str],
    h_keys: list[tuple[str, str, str]],
    own_key: tuple[str, str, str],
    cutoff: float,
) -> tuple[list[str], list[str]]:
    """H/D atoms within ``cutoff`` of ``point``, split by ownership.

    Returns (names of ``own_key``'s residue, names of any other). Only an own hydrogen
    can protonate the ring; a foreign one within 1.2 A is a defect
    (n_centres_foreign_h)."""
    if len(h_xyz) == 0:
        return [], []
    d = np.linalg.norm(h_xyz - np.asarray(point, dtype=float), axis=1)
    own: list[str] = []
    foreign: list[str] = []
    for i in np.flatnonzero(d <= cutoff):
        (own if h_keys[i] == own_key else foreign).append(h_names[i])
    return own, foreign


def _tautomer_of(
    his: Residue,
    h_xyz: np.ndarray,
    h_names: list[str],
    h_keys: list[tuple[str, str, str]],
) -> tuple[str, dict[str, bool] | None, bool]:
    """Which ring nitrogen carries a hydrogen OF ITS OWN RESIDUE, BY DISTANCE.

    Returns (label, carries, saw_foreign_h). ``label`` is 'HD1', 'HE2', 'HD1+HE2' or
    'none'. ``carries`` maps each ring N to "carries an H" and decides the roles, but is
    None when NEITHER does: a ring with no H in a file that has hydrogens is a broken
    reading (counted in n_centres_no_ring_h) and falls back to the partner class,
    whereas HIP is a real state (both donors, n_centres_hip). ``saw_foreign_h`` flags
    another residue's H within TAUTOMER_H_DIST, never used. Call only when the file HAS
    hydrogens. Atom names are cross-checked, never trusted: a name/geometry disagreement
    warns on stderr only."""
    own_key = (his.chain, his.seq, his.icode)
    carries: dict[str, bool] = {}
    saw_foreign = False
    for n_name in RING_NEIGHBOURS:
        pos = his.atom(n_name)
        if pos is None:
            # Unreachable: _centre_nitrogens has already raised for a missing ring N.
            raise ValueError(
                f"histidine {his.label()} has no {n_name} atom; the tautomer cannot "
                f"be read"
            )
        bonded, foreign = _h_names_within(
            pos, h_xyz, h_names, h_keys, own_key, TAUTOMER_H_DIST
        )
        carries[n_name] = bool(bonded)
        if foreign:
            saw_foreign = True
            print(
                f"WARNING: histidine {his.label()} has {','.join(sorted(set(foreign)))}"
                f" of ANOTHER residue within {TAUTOMER_H_DIST} A of {n_name}; it is "
                f"NOT used for the tautomer (see n_centres_foreign_h) -- check that "
                f"centre for a clash",
                file=sys.stderr,
            )
        expected = EXPECTED_H_NAMES[n_name]
        if bonded and not (set(bonded) & expected):
            print(
                f"WARNING: histidine {his.label()} carries {','.join(sorted(bonded))} "
                f"within {TAUTOMER_H_DIST} A of {n_name}, none of them named "
                f"{'/'.join(sorted(expected))}; the tautomer was read BY DISTANCE and "
                f"says {n_name} is protonated",
                file=sys.stderr,
            )
        if not bonded:
            wider, _ = _h_names_within(
                pos, h_xyz, h_names, h_keys, own_key, NAME_CHECK_DIST
            )
            named = set(wider) & expected
            if named:
                print(
                    f"WARNING: histidine {his.label()} has {','.join(sorted(named))} "
                    f"between {TAUTOMER_H_DIST} and {NAME_CHECK_DIST} A of {n_name}; "
                    f"the name says that nitrogen is protonated and the geometry says "
                    f"it is not -- the distance reading wins and {n_name} is read as "
                    f"UNprotonated",
                    file=sys.stderr,
                )
    if carries["ND1"] and carries["NE2"]:
        return "HD1+HE2", carries, saw_foreign
    if carries["ND1"]:
        return "HD1", carries, saw_foreign
    if carries["NE2"]:
        return "HE2", carries, saw_foreign
    # No ring hydrogen at all: report the defect, but give the role reader NOTHING.
    print(
        f"WARNING: histidine {his.label()} carries no hydrogen within "
        f"{TAUTOMER_H_DIST} A of ND1 or NE2 although this file has hydrogens; its "
        f"roles fall back to the partner class (see n_centres_no_ring_h)",
        file=sys.stderr,
    )
    return "none", None, saw_foreign


def _contact_roles(
    per_n: dict[str, NContacts],
    classes: list[str],
    carries: dict[str, bool] | None,
) -> dict[str, int]:
    """What each QUALIFYING contact of one centre is, counted by role.

    The entries are those in per_centre_nd1_target / _ne2_target, so the five counts sum
    to that centre's entries. With ``carries`` the His N is the donor if it carries the
    H, else the acceptor, and a same-kind pairing is INCOHERENT (a trust signal, no
    role). Without it the partner class decides: donor-only means the His accepted,
    acceptor-only that it donated, both-capable is AMBIGUOUS. Unclassified partners are
    counted, given no role."""
    roles = {field: 0 for field in ROLE_FIELDS}
    for n_name, contacts in per_n.items():
        n_carries_h = None if carries is None else carries[n_name]
        for i in np.flatnonzero(contacts.hb):
            partner = classes[i]
            if partner == UNCLASSIFIED:
                roles["unclass"] += 1
            elif n_carries_h is None:
                # Partner-class mode: the His N's own state is unknown.
                if partner == DONOR_ONLY:
                    roles["acc"] += 1
                elif partner == ACCEPTOR_ONLY:
                    roles["don"] += 1
                else:
                    roles["amb"] += 1
            elif n_carries_h:
                # The His N carries the H: it can only donate.
                roles[
                    "don" if partner in (ACCEPTOR_ONLY, BOTH_CAPABLE) else "incoh"
                ] += 1
            else:
                # Lone-pair N: acid protonates it, so acc is the pH-sensitive role.
                roles["acc" if partner in (DONOR_ONLY, BOTH_CAPABLE) else "incoh"] += 1
    return roles


def _ring_centroid(his: Residue) -> np.ndarray:
    """Mean of the five imidazole ring atoms (all required by _centre_nitrogens)."""
    points = []
    for name in RING_ATOMS:
        pos = his.atom(name)
        if pos is None:
            raise ValueError(
                f"histidine {his.label()} has no {name} atom; the ring centroid and "
                f"its pH environment cannot be computed"
            )
        points.append(pos)
    return np.asarray(points, dtype=float).mean(axis=0)


def _charged_atoms(residues: list[Residue]) -> tuple[np.ndarray, np.ndarray]:
    """Acidic and basic side-chain ATOM coordinates: an Asp counts two, an Arg three."""
    acidic: list[tuple[float, float, float]] = []
    basic: list[tuple[float, float, float]] = []
    for res in residues:
        key = res.resname.strip().upper()
        for i, name in enumerate(res.names):
            pair = (key, name.strip().upper())
            if pair in ENV_ACIDIC_ATOMS:
                acidic.append(res.xyz[i])
            elif pair in ENV_BASIC_ATOMS:
                basic.append(res.xyz[i])
    return (
        np.asarray(acidic, dtype=float).reshape(-1, 3),
        np.asarray(basic, dtype=float).reshape(-1, 3),
    )


def _count_within(centre: np.ndarray, points: np.ndarray, radius: float) -> int:
    """How many ``points`` lie within ``radius`` of ``centre`` (0 is a measurement)."""
    if len(points) == 0:
        return 0
    return int((np.linalg.norm(points - centre, axis=1) <= radius).sum())


def _ordinal_in(residues: list[Residue], target: Residue) -> int:
    """The 1-BASED position of ``target`` in ``residues``, by object identity.

    The numbering atomium/ProteinMPNN use (every chain renumbered 1..L). The caller
    passes ``binder_res[idx1 - 1]``, so this is always the supplied position: it ties
    the value to the residue the measurement ran on but is NOT a check. The His-name
    check and ``--expect-binder-len-col`` stand behind idx1."""
    for i, res in enumerate(residues):
        if res is target:
            return i + 1
    raise ValueError(
        f"residue {target.label()} is not in the binder chain's residue list; the "
        f"centre position cannot be established"
    )


def _resnum_of(residue: Residue) -> int | None:
    """The residue's auth_seq_id as an integer, or None when the file gave none.

    An insertion code is dropped (``112A`` -> 112). None is NA, never 0, and forces
    centres_idx1_eq_resnum to 0: an agreement that cannot be shown is not one.
    """
    text = residue.seq.strip()
    try:
        return int(text)
    except ValueError:
        return None


def _fmt_dist(value: float | None) -> str:
    return "NA" if value is None else f"{value:.2f}"


def measure(
    path: Path,
    binder_chain: str,
    target_chains: list[str],
    centres: list[int],
    min_dist: float,
    max_dist: float,
    max_angle: float,
    env_radius: float,
    expect_len: int | None,
    expect_col: str,
    metrics: dict,
) -> dict:
    """One structure + one centre list -> metrics.

    Raises on anything that makes the measurement or the MAPPING untrustworthy; main()
    records it as the status with every science column NA. ``metrics`` is filled in
    place, so an error row keeps the facts known before the failure (chains found,
    binder length)."""
    found = read_structure(path)
    metrics["chain_ids_found"] = ",".join(found.chains_present)
    metrics["structure_has_hydrogens"] = found.has_hydrogens

    if binder_chain not in found.residues:
        raise ValueError(
            f"binder chain {binder_chain} absent from {path.name} "
            f"(have: {','.join(found.chains_present) or 'none'})"
        )
    absent = [c for c in target_chains if c not in found.residues]
    if absent:
        raise ValueError(
            f"target chain(s) {','.join(absent)} absent from {path.name} "
            f"(have: {','.join(found.chains_present) or 'none'})"
        )

    # Binder residues in file order (what a centre index indexes), minus het groups.
    binder_res = [
        r for r in found.residues[binder_chain] if not _excluded_residue(r.resname)
    ]
    if not binder_res:
        raise ValueError(
            f"binder chain {binder_chain} of {path.name} has no amino-acid residues"
        )
    metrics["binder_chain_n_res"] = len(binder_res)

    if expect_len is not None:
        metrics["expect_binder_len"] = expect_len
        if expect_len != len(binder_res):
            raise ValueError(
                f"binder chain {binder_chain} of {path.name} has {len(binder_res)} "
                f"residues but {expect_col} says {expect_len}; the centre indices "
                f"would land on the wrong residues"
            )

    # Same het filter as the binder; target_n_het_atoms makes an ignored ligand visible.
    target_all = [r for c in target_chains for r in found.residues[c]]
    target_res = [r for r in target_all if not _excluded_residue(r.resname)]
    metrics["target_chain_n_res"] = len(target_res)
    metrics["target_n_het_atoms"] = sum(
        len(r.names) for r in target_all if _excluded_residue(r.resname)
    )
    target_xyz, target_labels, target_classes = _polar_atoms(target_res)
    if len(target_xyz) == 0:
        raise ValueError(
            f"target chain(s) {','.join(target_chains)} of {path.name} hold no N/O "
            f"heavy atoms in an amino-acid residue; this is not a protein target"
        )

    # pLDDT scale from the file; an all-zero B column is 'unknown', not '0-1'.
    if found.max_bfac is None or found.max_bfac <= 0.0:
        scale = "unknown"
    elif found.max_bfac <= PLDDT_UNIT_MAX:
        scale = "0-1"
    else:
        scale = "0-100"
    metrics["plddt_scale"] = scale

    # Charged TARGET atoms only: the binder's own charges are also in the free binder.
    env_acidic_xyz, env_basic_xyz = _charged_atoms(target_res)
    # Every H/D in the file, for the tautomer reading only (never for a contact).
    h_xyz = np.asarray(found.h_xyz, dtype=float).reshape(-1, 3)

    per_centre: list[str] = []
    # Per nitrogen (target side): centre count, minimum distance and partner cells.
    n_hb_target_n = {n: 0 for n in RING_NEIGHBOURS}
    min_target_n: dict[str, float | None] = {n: None for n in RING_NEIGHBOURS}
    partners_target_n: dict[str, list[str]] = {n: [] for n in RING_NEIGHBOURS}
    # The same partner cells for the intra-binder side, naming the binder partners.
    partners_binder_n: dict[str, list[str]] = {n: [] for n in RING_NEIGHBOURS}
    n_hb_target = 0
    n_hb_binder = 0
    n_near_target = 0
    n_close_target = 0
    n_close_binder = 0
    min_target: float | None = None
    min_binder: float | None = None
    plddt_min: float | None = None
    n_his_ok = 0
    # (idx1, resnum) per reported centre; resnum is None without a usable integer.
    centre_positions: list[tuple[int, int | None]] = []
    # (tautomer, n_acidic, n_basic) per reported centre; None tautomer = no hydrogens.
    centre_extras: list[tuple[str | None, int, int]] = []
    # The role tally, over all centres (not capped): these feed row-level counts.
    role_entries: list[str] = []
    n_acceptor_role = 0
    n_donor_role = 0
    n_switch_candidate = 0
    n_switch_candidate_clean = 0
    n_incoherent = 0
    n_ambiguous = 0
    n_unclassified = 0
    # The same tally for the intra-binder side.
    role_entries_binder: list[str] = []
    n_acceptor_role_binder = 0
    n_donor_role_binder = 0
    n_switch_candidate_binder = 0
    n_switch_candidate_clean_binder = 0
    n_incoherent_binder = 0
    n_ambiguous_binder = 0
    n_unclassified_binder = 0
    # Binder entries over all centres: the denominator of the per-centre role counts.
    n_contacts_binder_entries = 0
    # Across both sides: acceptor-role somewhere, no donor-role or incoherent anywhere.
    n_switch_candidate_clean_any = 0
    n_with_tautomer = 0
    n_no_ring_h = 0
    n_hip = 0
    n_foreign_h = 0

    # Designed centres as 0-based list offsets, for the '*' partner mark only.
    designed_idx = frozenset(p - 1 for p in centres)
    for idx1 in centres:
        idx = idx1 - 1  # 1-based position -> offset into binder_res (file order)
        if idx >= len(binder_res):
            raise ValueError(
                f"centre position {idx1} is past the end of binder chain "
                f"{binder_chain}, which has {len(binder_res)} residues in {path.name}"
            )
        his = binder_res[idx]
        if his.resname.strip().upper() not in HIS_NAMES:
            raise ValueError(
                f"centre position {idx1} of binder chain {binder_chain} is "
                f"{his.resname} {his.label()}, not a histidine -- either the positions "
                f"do not map onto this structure, or they are not His centres at all "
                f"(an acid_p / all_p design pins Asp/Glu and is out of scope here)"
            )
        n_his_ok += 1
        # Updated as we go: an error row shows how many centres were confirmed His.
        metrics["centers_his_ok"] = n_his_ok

        # Position of `his`, the verified residue (_ordinal_in is not itself a check).
        report_this_centre = len(centre_positions) < MAX_CENTRE_COLUMNS
        if report_this_centre:
            centre_positions.append((_ordinal_in(binder_res, his), _resnum_of(his)))

        # The binder side drops the centre residue and only the BACKBONE of i-1 / i+1.
        binder_xyz, binder_labels, binder_classes = _binder_partner_atoms(
            binder_res, idx, designed_idx
        )

        target_n = _centre_nitrogens(his, target_xyz, min_dist, max_dist, max_angle)
        binder_n = _centre_nitrogens(his, binder_xyz, min_dist, max_dist, max_angle)
        hb_t, near_t, close_t, d_t = _merge_nitrogens(target_n, len(target_xyz))
        hb_b, near_b, close_b, d_b = _merge_nitrogens(binder_n, len(binder_xyz))

        # The per-nitrogen view of the SAME masks the merge above consumed.
        for n_name, contacts in target_n.items():
            if bool(contacts.hb.any()):
                n_hb_target_n[n_name] += 1
            if contacts.min_d is not None:
                current = min_target_n[n_name]
                min_target_n[n_name] = (
                    contacts.min_d if current is None else min(current, contacts.min_d)
                )
            partners_target_n[n_name].append(
                _partner_entries(idx1, contacts, target_labels)
            )
        for n_name, contacts in binder_n.items():
            partners_binder_n[n_name].append(
                _partner_entries(idx1, contacts, binder_labels)
            )

        # Roles of the target contacts (tautomer if the file has H, else partner class).
        carries: dict[str, bool] | None = None
        tautomer: str | None = None
        if found.has_hydrogens:
            tautomer, carries, foreign_h = _tautomer_of(
                his, h_xyz, found.h_names, found.h_keys
            )
            if carries is None:
                # 'none': broken reading, falls back to the partner class.
                n_no_ring_h += 1
            else:
                n_with_tautomer += 1
            if tautomer == "HD1+HE2":
                n_hip += 1
            if foreign_h:
                n_foreign_h += 1
        roles = _contact_roles(target_n, target_classes, carries)
        role_entries.append(
            f"{idx1}:" + ":".join(f"{f}={roles[f]}" for f in ROLE_FIELDS)
        )
        n_acceptor_role += 1 if roles["acc"] > 0 else 0
        n_donor_role += 1 if roles["don"] > 0 else 0
        # Candidate: >= 1 acceptor-role, no donor-role; amb/incoh/unclass only reported.
        n_switch_candidate += 1 if (roles["acc"] > 0 and roles["don"] == 0) else 0
        # Clean: also no incoherent contact (it could invert the verdict); gate on it.
        n_switch_candidate_clean += (
            1 if (roles["acc"] > 0 and roles["don"] == 0 and roles["incoh"] == 0) else 0
        )
        n_incoherent += roles["incoh"]
        n_ambiguous += roles["amb"]
        n_unclassified += roles["unclass"]

        # Intra-binder roles: same function and tautomer reading, binder-side masks.
        roles_b = _contact_roles(binder_n, binder_classes, carries)
        role_entries_binder.append(
            f"{idx1}:" + ":".join(f"{f}={roles_b[f]}" for f in ROLE_FIELDS)
        )
        n_contacts_binder_entries += sum(roles_b[f] for f in ROLE_FIELDS)
        n_acceptor_role_binder += 1 if roles_b["acc"] > 0 else 0
        n_donor_role_binder += 1 if roles_b["don"] > 0 else 0
        n_switch_candidate_binder += (
            1 if (roles_b["acc"] > 0 and roles_b["don"] == 0) else 0
        )
        n_switch_candidate_clean_binder += (
            1
            if (roles_b["acc"] > 0 and roles_b["don"] == 0 and roles_b["incoh"] == 0)
            else 0
        )
        n_incoherent_binder += roles_b["incoh"]
        n_ambiguous_binder += roles_b["amb"]
        n_unclassified_binder += roles_b["unclass"]
        # Across both sides, per centre: neither the AND nor the OR of the two verdicts.
        n_switch_candidate_clean_any += (
            1
            if (
                (roles["acc"] + roles_b["acc"]) > 0
                and (roles["don"] + roles_b["don"]) == 0
                and (roles["incoh"] + roles_b["incoh"]) == 0
            )
            else 0
        )

        if report_this_centre:
            # pH context around the ring centroid; 0 is a measurement, never NA.
            centroid = _ring_centroid(his)
            centre_extras.append(
                (
                    tautomer,
                    _count_within(centroid, env_acidic_xyz, env_radius),
                    _count_within(centroid, env_basic_xyz, env_radius),
                )
            )

        n_hb_target += 1 if hb_t > 0 else 0
        n_hb_binder += 1 if hb_b > 0 else 0
        n_near_target += 1 if near_t > 0 else 0
        n_close_target += 1 if close_t > 0 else 0
        n_close_binder += 1 if close_b > 0 else 0
        if d_t is not None:
            min_target = d_t if min_target is None else min(min_target, d_t)
        if d_b is not None:
            min_binder = d_b if min_binder is None else min(min_binder, d_b)

        for i, name in enumerate(his.names):
            if name in BACKBONE_ATOMS:
                continue
            b = his.bfac[i]
            if b is not None:
                plddt_min = b if plddt_min is None else min(plddt_min, b)

        per_centre.append(
            f"{idx1}:{his.resname.strip().upper()}"
            f":target_min={_fmt_dist(d_t)}:n_partner_atoms={hb_t}"
            f":binder_min={_fmt_dist(d_b)}:n_partner_atoms_binder={hb_b}"
            f":n_near_target={near_t}:n_near_binder={near_b}"
            f":n_too_close_target={close_t}:n_too_close_binder={close_b}"
        )

    if len(centres) > MAX_CENTRE_COLUMNS:
        # Extra centres are counted everywhere but get no position column: say so.
        print(
            f"WARNING: {path.name} has {len(centres)} centres and only the first "
            f"{MAX_CENTRE_COLUMNS} get centre<k>_idx1 / centre<k>_resnum / "
            f"centre<k>_tautomer / centre<k>_n_acidic_env / centre<k>_n_basic_env "
            f"columns (compare hisbond centres_reported with centers_n, and read "
            f"centers_spec for the full list); every other column covers all of them",
            file=sys.stderr,
        )

    metrics.update(
        {
            "centers_his_ok": n_his_ok,
            "per_centre": ";".join(per_centre),
            "n_centres_hb_target": n_hb_target,
            "any_target_hb_proxy": n_hb_target > 0,
            "min_dist_target": min_target,
            # Per nitrogen, target side: each <= merged <= their sum.
            **{
                f"n_centres_hb_target_{suffix}": n_hb_target_n[n_name]
                for n_name, suffix in NITROGEN_SUFFIX.items()
            },
            **{
                f"any_target_hb_proxy_{suffix}": n_hb_target_n[n_name] > 0
                for n_name, suffix in NITROGEN_SUFFIX.items()
            },
            **{
                f"min_dist_target_{suffix}": min_target_n[n_name]
                for n_name, suffix in NITROGEN_SUFFIX.items()
            },
            **{
                f"per_centre_{suffix}_target": ";".join(partners_target_n[n_name])
                for n_name, suffix in NITROGEN_SUFFIX.items()
            },
            # Their intra-binder twins ('*' marks a designed-centre partner).
            **{
                f"per_centre_{suffix}_binder": ";".join(partners_binder_n[n_name])
                for n_name, suffix in NITROGEN_SUFFIX.items()
            },
            # What those contacts are: five counts per centre partitioning its entries.
            "per_centre_role_target": ";".join(role_entries),
            "n_centres_acceptor_role_target": n_acceptor_role,
            "n_centres_donor_role_target": n_donor_role,
            "n_centres_switch_candidate": n_switch_candidate,
            "n_centres_switch_candidate_clean": n_switch_candidate_clean,
            "any_switch_candidate": n_switch_candidate > 0,
            # Trust counts: contacts with no role (incoherent, ambiguous, unclassified).
            "n_contacts_incoherent_target": n_incoherent,
            "n_contacts_ambiguous_target": n_ambiguous,
            "n_contacts_unclassified_target": n_unclassified,
            # The same reading of the intra-binder contacts n_centres_hb_binder counts.
            "per_centre_role_binder": ";".join(role_entries_binder),
            "n_centres_acceptor_role_binder": n_acceptor_role_binder,
            "n_centres_donor_role_binder": n_donor_role_binder,
            "n_centres_switch_candidate_binder": n_switch_candidate_binder,
            "n_centres_switch_candidate_clean_binder": n_switch_candidate_clean_binder,
            "any_switch_candidate_binder": n_switch_candidate_binder > 0,
            "n_contacts_incoherent_binder": n_incoherent_binder,
            "n_contacts_ambiguous_binder": n_ambiguous_binder,
            "n_contacts_unclassified_binder": n_unclassified_binder,
            "n_contacts_binder_entries": n_contacts_binder_entries,
            # Both sides at once, per centre: something to break and nothing to keep.
            "n_centres_switch_candidate_clean_any": n_switch_candidate_clean_any,
            "any_switch_candidate_clean_any": n_switch_candidate_clean_any > 0,
            # Reading defects as row totals (cover centres past the cap); expected 0.
            "n_centres_no_ring_h": n_no_ring_h,
            "n_centres_hip": n_hip,
            "n_centres_foreign_h": n_foreign_h,
            # Evidence behind the roles; do not pool 'mixed' rows with either mode.
            "role_basis": (
                "tautomer"
                if n_with_tautomer == len(centres)
                else "partner_class"
                if n_with_tautomer == 0
                else "mixed"
            ),
            # Positions in the centres column's order; NA past the row's own centres.
            "centres_reported": len(centre_positions),
            "centres_idx1_eq_resnum": int(
                bool(centre_positions)
                and all(rn is not None and idx1 == rn for idx1, rn in centre_positions)
            ),
            **{
                f"centre{k}_{field}": value
                for k, (idx1, resnum) in enumerate(centre_positions, start=1)
                for field, value in (("idx1", idx1), ("resnum", resnum))
            },
            # Tautomer and environment of the same centres (NA tautomer: no hydrogens).
            **{
                f"centre{k}_{field}": value
                for k, (tautomer, n_acid, n_base) in enumerate(centre_extras, start=1)
                for field, value in (
                    ("tautomer", tautomer),
                    ("n_acidic_env", n_acid),
                    ("n_basic_env", n_base),
                )
            },
            "n_centres_hb_binder": n_hb_binder,
            "min_dist_binder": min_binder,
            "n_centres_polar_dist_only_target": n_near_target,
            "n_centres_too_close_target": n_close_target,
            "n_centres_too_close_binder": n_close_binder,
            # NA, not 0.0, when there was no usable B-factor ('unknown' scale).
            "his_plddt_min": None if scale == "unknown" else plddt_min,
        }
    )
    return metrics


# --- output ---


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def write_result(
    out_dir: Path, name: str, status: str, path: str, metrics: dict
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / f"{name}.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(RESULT_COLUMNS)
        w.writerow(
            [name, status, path, *(_fmt(metrics.get(c)) for c in METRIC_COLUMNS)]
        )


def _parse_expect(raw: str, column: str) -> int | None:
    """An expected binder length from a table cell: an integer, or a sequence's length.

    Anything else raises: a silently ignored check is worse than no check.
    """
    text = raw.strip()
    if not text:
        return None
    if text == UNRESOLVED_EXPECT:
        raise ValueError(
            f"the binder-length check was asked for but {column} resolves nowhere "
            f"for this design (neither on its own row nor on any ancestor); a check "
            f"that cannot be made is an error, not a skipped check"
        )
    try:
        return int(text)
    except ValueError:
        pass
    # '60.0' is the same length: pandas floats an int column once any row is NA.
    try:
        value = float(text)
    except ValueError:
        pass
    else:
        if value.is_integer():
            return int(value)
        raise ValueError(
            f"{column} holds {raw!r}, which is not a whole number of residues; a "
            f"binder length cannot be fractional"
        )
    seq = text.replace(" ", "").upper()
    if seq.isalpha():
        return len(seq)
    raise ValueError(
        f"{column} holds {raw!r}, which is neither an integer length nor a plain "
        f"one-letter sequence; the binder-length check cannot be made"
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Do the designed His centres hydrogen-bond the target?"
    )
    ap.add_argument(
        "--task-file",
        type=Path,
        required=True,
        help="Tab-separated 'name<TAB>structure<TAB>centres<TAB>centres_depth<TAB>"
        "expect_col<TAB>expect_value<TAB>centres_error' rows, one per design, SEVEN "
        "fields always. 'centres' is the resolved comma-separated 1-based position "
        "list and 'centres_depth' how far up the lineage it was resolved (0 = the "
        "design's own table); 'centres_error' carries the submit-time reason it "
        "could not be resolved and is re-raised as this design's status (it and "
        "'centres' are mutually exclusive).",
    )
    ap.add_argument(
        "--binder-chain",
        required=True,
        help="The binder chain ID, in the auth_asym_id namespace.",
    )
    ap.add_argument(
        "--target-chains", required=True, help="Comma-joined target chains."
    )
    ap.add_argument(
        "--min-dist",
        type=float,
        required=True,
        help="Distance floor (A): a partner closer than this is a clash, "
        "not a hydrogen bond, and is excluded from both counts. Required, not "
        "defaulted: hisbond.sh always passes it and the value is echoed into a "
        "column, so a silent fallback would mislabel the window.",
    )
    ap.add_argument("--max-dist", type=float, required=True)
    ap.add_argument("--max-angle", type=float, required=True)
    ap.add_argument(
        "--env-radius",
        type=float,
        required=True,
        help="Radius (A) around a centre's imidazole ring centroid in which TARGET "
        "acidic (ASP OD1/OD2, GLU OE1/OE2) and basic (LYS NZ, ARG NE/NH1/NH2) ATOMS "
        "are counted as centre<k>_n_acidic_env / centre<k>_n_basic_env. Required "
        f"here; the tool's default ({DEFAULT_ENV_RADIUS}) lives on run_hisbond.py's "
        "--env-radius. Must be > 0. Echoed as the env_radius column.",
    )
    ap.add_argument(
        "--centers-arg",
        default="",
        help="The --centers template verbatim, echoed as the hisbond_centers_arg "
        "column: the provenance of every position (e.g. "
        "'{atomiumph_centre1_idx1},{atomiumph_centre2_idx1}'). Never parsed here -- "
        "the task file carries the already-resolved positions.",
    )
    ap.add_argument(
        "--input-column",
        default="the input structure column",
        help="Only used to name the column in the 'could not resolve' error.",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(namespace=Args())

    target_chains = split_list(args.target_chains)
    if not target_chains:
        raise ValueError("--target-chains is empty")
    binder_chain = args.binder_chain.strip()
    if not binder_chain:
        raise ValueError("--binder-chain is empty")
    if binder_chain in target_chains:
        raise ValueError(
            f"--binder-chain {binder_chain} is also in --target-chains; the binder "
            f"cannot be its own binding partner"
        )
    if args.max_dist <= 0:
        raise ValueError("--max-dist must be > 0")
    if args.min_dist < 0:
        raise ValueError("--min-dist must be >= 0")
    if args.min_dist >= args.max_dist:
        raise ValueError(
            f"--min-dist ({args.min_dist}) must be below --max-dist "
            f"({args.max_dist}); otherwise the distance window is empty and no "
            f"contact can ever qualify"
        )
    if not 0 < args.max_angle <= 180:
        raise ValueError("--max-angle must be in (0, 180]")
    if args.env_radius <= 0:
        raise ValueError(
            f"--env-radius must be > 0 A, got {args.env_radius}; a non-positive "
            f"radius counts nothing and the environment columns would read as an "
            f"empty neighbourhood rather than as a refused question"
        )

    echo = {
        "binder_chain": binder_chain,
        "target_chains": ",".join(target_chains),
        "centers_arg": args.centers_arg,
        "min_dist": args.min_dist,
        "max_dist": args.max_dist,
        "max_angle": args.max_angle,
        "env_radius": args.env_radius,
    }

    with open(args.task_file) as f:
        members = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    for member in members:
        if len(member) != 7:
            raise ValueError(
                f"malformed task row {member!r}: expected exactly 7 tab-separated "
                f"fields 'name<TAB>structure<TAB>centres<TAB>centres_depth<TAB>"
                f"expect_col<TAB>expect_value<TAB>centres_error', got {len(member)}"
            )
        (
            name,
            src,
            centres_cell,
            centres_depth,
            expect_col,
            expect_raw,
            centres_error,
        ) = member
        expect_col = expect_col or "the expected-binder-length column"
        metrics = dict(echo)
        try:
            # THE CENTRES FIRST, and before the structure: they do not depend on it,
            # and a reader interrogating a structure-less error row needs to see WHICH
            # positions it was asked about (centers_spec / centers_n / the depth).
            if centres_error.strip():
                # Resolved at submit time and failed there; the reason rides through.
                raise ValueError(
                    f"could not resolve --centers {args.centers_arg!r} for this "
                    f"design: {centres_error.strip()}"
                )
            if not centres_cell.strip():
                raise ValueError(
                    f"no centres reached this design from --centers "
                    f"{args.centers_arg!r}; there is nothing to measure"
                )
            centres = parse_centres(centres_cell)
            # Every supplied position, past MAX_CENTRE_COLUMNS too (trust column).
            metrics["centers_spec"] = ",".join(str(p) for p in centres)
            metrics["centers_n"] = len(centres)
            if centres_depth.strip():
                metrics["centers_lineage_depth"] = int(centres_depth)
            if not src.strip():
                raise ValueError(
                    f"could not resolve a structure path from column "
                    f"{args.input_column!r} for this design -- neither on its own row "
                    f"nor on any ancestor (parent_name/parent_table)"
                )
            expect_len = _parse_expect(expect_raw, expect_col)
            measure(
                Path(src),
                binder_chain,
                target_chains,
                centres,
                args.min_dist,
                args.max_dist,
                args.max_angle,
                args.env_radius,
                expect_len,
                expect_col,
                metrics,
            )
            write_result(args.out_dir, name, "OK", src, metrics)
            print(
                f"{name}: centres={metrics['centers_n']} "
                f"at={metrics['centers_spec']} "
                f"depth={metrics.get('centers_lineage_depth')} "
                f"his_ok={metrics['centers_his_ok']} "
                f"any_target_hb_proxy={metrics['any_target_hb_proxy']} "
                f"n_centres_hb_target={metrics['n_centres_hb_target']} "
                f"(nd1={metrics['n_centres_hb_target_nd1']} "
                f"ne2={metrics['n_centres_hb_target_ne2']}) "
                f"min_dist_target={_fmt_dist(metrics['min_dist_target'])} "
                f"too_close={metrics['n_centres_too_close_target']} "
                f"role_basis={metrics['role_basis']} "
                f"switch={metrics['n_centres_switch_candidate']} "
                f"clean={metrics['n_centres_switch_candidate_clean']} "
                f"incoh={metrics['n_contacts_incoherent_target']} "
                f"amb={metrics['n_contacts_ambiguous_target']} "
                f"unclass={metrics['n_contacts_unclassified_target']} "
                f"no_ring_h={metrics['n_centres_no_ring_h']} "
                f"hip={metrics['n_centres_hip']} "
                f"foreign_h={metrics['n_centres_foreign_h']} "
                f"| {metrics['per_centre']} "
                f"| roles {metrics['per_centre_role_target']} "
                f"| ND1 {metrics['per_centre_nd1_target']} "
                f"| NE2 {metrics['per_centre_ne2_target']}"
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            msg = " ".join(str(e).split())  # one line: it lands in a TSV cell
            print(f"{name}: ERROR {msg}", file=sys.stderr)
            # Echo and trust facts survive; every science column is empty (NA, never 0).
            kept = {
                k: v
                for k, v in metrics.items()
                if k in ECHO_COLUMNS or k in TRUST_COLUMNS
            }
            write_result(args.out_dir, name, f"error: {msg}", "", kept)


if __name__ == "__main__":
    main()
