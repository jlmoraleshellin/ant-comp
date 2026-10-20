#!/usr/bin/env python3
"""
Per-design worker: what pKa does PROPKA3 predict for each designed histidine of the
binder, in the complex and in the free binder, and how big is the shift? How to read
the columns and the traps: ``.claude/skills/propka/SKILL.md``.

Two PROPKA3 runs per design, on temporary copies of ONE input file:

  complex   every chain of the file, hydrogens stripped
  free      ONLY the binder chains' records of the SAME file, SAME COORDINATES,
            hydrogens stripped: a RIGID unbinding, recorded as
            free_state_model = 'chain_only_same_coordinates'

A file with no chain outside the binder chains is an error: the two runs would be the
same molecule and every dpka 0. ``--free-only`` is the explicit single-run mode: only
the free run is made, every complex-dependent column is NA by construction and
``free_only`` is 1.

dpka = pka_complex - pka_free; NEGATIVE is the wanted sign for bind-at-neutral,
release-in-acid (the complex disfavours the protonated His more than the free binder
does). For one residue dG_bind(pH) - dG_bind(inf. pH) = RT ln[(1 + 10^(pKa_free - pH))
/ (1 + 10^(pKa_complex - pH))]; the F, switch and release columns are arithmetic on the
two pKa values, with no third run. The stab_* block asks a separate question on the
free state only: does protonating the designed histidines destabilise the binder's own
fold? It is non-fatal: a failure leaves stab_* NA with stab_error set and the row OK.

Centres are histidines named by (chain, residue number) in auth_seq_id numbering: a spec
(a literal list; the builder resolves any {expr} islands) or 'auto' (every His of the
binder chains). Each must be a histidine, uniquely numbered in its chain, or the design
is an error row.

PROPKA is empirical (about one pK unit of error) and is run on predicted or relaxed
coordinates. It does not assign the tautomer (hisbond's centre<k>_tautomer does).

THE INPUT FILE IS NEVER MODIFIED. PROPKA reads temporary copies (under TMPDIR, kept
with --keep-propka-input) that carry four edits: (1) every H/D atom dropped (element
columns 77-78; on ATOM records only, the first letter of the name when the element is
blank); (2) HID/HIE/HIP/HSD/HSE/HSP renamed HIS; (3) a bare TER at every chain change,
so the binder's N-terminus is searched afresh; (4) only the first model. HETATM records
are copied as they are (PROPKA ignores HOH and a few buffers, but any other ligand
becomes a group and can shift a pKa). Residue counting and His lookups use ATOM records.

Error-as-data: every failure is written as 'error: ...' with the trust columns kept and
every science column EMPTY. NA is not 0.
"""

import argparse
import csv
import gzip
import logging
import sys
import tempfile
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Tolerated spellings of a histidine; the temp files PROPKA reads carry 'HIS' for all.
HIS_NAMES = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}

# The first MAX_CENTRE_COLUMNS centres get per-centre columns (centres_reported says how
# many did); every other column covers all centres. THIS MODULE IS THE ONE DEFINITION:
# run_propka.py and collect_propka.py import it, and it must import only the standard
# library at module level because it also runs in the tool image, which has no prosapia.
MAX_CENTRE_COLUMNS = 4
CENTRE_FIELDS = (
    "chain",
    "resnum",
    "pka_complex",
    "pka_free",
    "dpka",
    "buried_complex",
    "buried_free",
    "n_coupled_complex",
    "n_coupled_free",
)
CENTRE_STR_FIELDS = frozenset({"chain"})
CENTRE_INT_FIELDS = frozenset({"resnum", "n_coupled_complex", "n_coupled_free"})
# Fields that need the complex run: NA by construction under --free-only.
CENTRE_COMPLEX_FIELDS = frozenset(
    {"pka_complex", "dpka", "buried_complex", "n_coupled_complex"}
)

# One-letter codes for binder_chains_sequence; anything else reads 'X'.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

# The centre spec that means "every histidine of the binder chains".
AUTO_CENTRES = "auto"

# How the free state was made, written as a column so nobody has to infer it.
FREE_STATE_MODEL = "chain_only_same_coordinates"

# Written by run_propka.py into the sub-manifest when --expect-binder-len-col resolves
# nowhere for a design: a check asked for and not made is an error row.
UNRESOLVED_EXPECT = "__unresolved__"

# ---------------------------------------------------------------------------
# The pH-dependence model (arithmetic on the two runs' pKa values; no third run)
# ---------------------------------------------------------------------------
# For ONE histidine, the two-state independent-sites factor by which binding is
# WEAKENED at pH relative to the high-pH limit:
#
#     F_i(pH) = (1 + 10^(pKa_free_i - pH)) / (1 + 10^(pKa_complex_i - pH))
#
# Independent sites multiply: F = prod_i F_i. F -> 1 at high pH and -> 10^(-dpka) at low
# pH, and since the denominator is >= 1, F_i(pH) <= 1 + 10^(pKa_free_i - pH) (the
# CEILING), a bound set by the free binder alone. The SWITCH RATIO F(pH) / F(pH_ref) is
# K_D(pH) / K_D(pH_ref) under the same model, and is what a pH switch is judged on.
F_MODEL = "two_state_independent_sites"

# Defaults of the pH layer. A pH that leaves the list stops being produced; its old
# columns go stale in a table that holds them.
DEFAULT_PH_VALUES = (7.4, 6.5, 6.0, 5.5)
DEFAULT_RELEASE_FACTORS = (3.0, 10.0)
# The denominator of the switch ratio: the pH the binder is meant to HOLD at.
DEFAULT_REFERENCE_PH = 7.4
MAX_PH_VALUES = 6
MAX_RELEASE_FACTORS = 4
SITE_SETS = ("designed", "all")

# The titration midpoint is found by bisection on [0, 14] (strictly decreasing in pH).
MIDPOINT_PH_LO = 0.0
MIDPOINT_PH_HI = 14.0
MIDPOINT_TOL = 1e-6

# The release pH is a downward grid scan, never a root: F is not monotonic in pH.
RELEASE_PH_HI = 9.0
RELEASE_PH_LO = 3.0
RELEASE_PH_STEP = 0.01

# ---------------------------------------------------------------------------
# FOLD STABILITY vs pH (the stab_* columns): FREE STATE ONLY
# ---------------------------------------------------------------------------
# PROPKA's own folding free energy (reference 'neutral'): ConformationContainer.
# calculate_folding_energy, the per-group sum of Group.calculate_folding_energy, evaluated
# directly at each pH. SIGN: a LOWER dG_fold is MORE stable (propka's optimum is the
# MINIMUM of get_folding_profile), so stab_loss_<tag> = dG_fold(pH) - dG_fold(reference)
# is POSITIVE when the binder is less stable at the lower pH. The split is designed
# centres vs EVERYTHING ELSE TITRATABLE in the binder (any non-designed His included).
STAB_PREFIX = "stab_"
STAB_REFERENCE = "neutral"
# PROPKA's own default profile grid, passed explicitly so stab_opt_* cannot drift.
STAB_PROFILE_GRID = (0.0, 14.0, 0.1)
# The integer pHs the .pka file prints: proof the direct evaluation is the same function.
STAB_CHECK_PHS = (5.0, 6.0, 7.0, 8.0, 9.0)
STAB_GRID_TOL = 0.051


def num_token(value: float, keep_decimal: bool) -> str:
    """A number as a canonical, collision-free column token.

    pH values keep at least one decimal (6.0 -> '6.0'), so ph6p0 reads alongside
    ph6p5; release factors drop a bare '.0' (3.0 -> '3', 2.5 -> '2.5'). Two decimals
    is the resolution of both -- finer than the pKa values themselves are worth.
    """
    text = f"{value:.2f}".rstrip("0")
    if text.endswith("."):
        text = text + "0" if keep_decimal else text[:-1]
    return text


def ph_tag(value: float) -> str:
    """``6.5`` -> ``'ph6p5'``."""
    return "ph" + num_token(value, keep_decimal=True).replace(".", "p")


def factor_tag(value: float) -> str:
    """``3.0`` -> ``'x3'``; ``2.5`` -> ``'x2p5'``."""
    return "x" + num_token(value, keep_decimal=False).replace(".", "p")


def parse_ph_values(text: str) -> list[float]:
    """``'7.4,6.5'`` -> ``[7.4, 6.5]``, with every rule the flag promises.

    Validated HERE as well as at submit time: the worker is runnable standalone, and
    a pH list that silently lost a value would produce a row whose columns do not
    mean what their names say.
    """
    values: list[float] = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = float(token)
        except ValueError as e:
            raise ValueError(f"--ph-values: {token!r} is not a number") from e
        if not 0.0 <= value <= 14.0:
            raise ValueError(
                f"--ph-values: {token!r} is outside pH 0-14; a pH outside that range "
                f"is not a question about water"
            )
        # Rounded to the two decimals the column TAG carries, so the ph_values echo is
        # exactly the question the numbers answer: '6.449' would be tagged ph6p45 and
        # then evaluated at 6.449, and no column would say so.
        values.append(round(value, 2))
    if not values:
        raise ValueError("--ph-values must name at least one pH")
    if len(values) > MAX_PH_VALUES:
        raise ValueError(
            f"--ph-values names {len(values)} pH values; at most {MAX_PH_VALUES} are "
            f"allowed (each one adds four columns to the table)"
        )
    tags = [ph_tag(v) for v in values]
    if len(set(tags)) != len(tags):
        raise ValueError(
            f"--ph-values lists the same pH twice ({','.join(tags)}); two identical "
            f"columns would differ only by name"
        )
    return values


def parse_release_factors(text: str) -> list[float]:
    """``'3,10'`` -> ``[3.0, 10.0]``. Every factor must be > 1."""
    values: list[float] = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = float(token)
        except ValueError as e:
            raise ValueError(f"--release-factors: {token!r} is not a number") from e
        if not value > 1.0:
            raise ValueError(
                f"--release-factors: {token!r} must be > 1; a factor of 1 or less is "
                f"not a weakening and the release pH would be meaningless"
            )
        # Rounded like the pH values, for the same reason: the release_factors echo is
        # the factor the grid scan actually compared against.
        values.append(round(value, 2))
    if not values:
        raise ValueError("--release-factors must name at least one factor")
    if len(values) > MAX_RELEASE_FACTORS:
        raise ValueError(
            f"--release-factors names {len(values)} factors; at most "
            f"{MAX_RELEASE_FACTORS} are allowed (each one adds two columns)"
        )
    tags = [factor_tag(v) for v in values]
    if len(set(tags)) != len(tags):
        raise ValueError(
            f"--release-factors lists the same factor twice ({','.join(tags)})"
        )
    return values


def parse_reference_ph(text: str, ph_values: list[float]) -> float:
    """The denominator pH of the switch ratio. MUST be one of ``ph_values``.

    Compared on the canonical TAG, never on float equality: '6.50' and 6.5 are the
    same column and must be the same reference. A reference outside the list is
    refused rather than silently added, because the ratio's denominator has to be a
    pH the row actually reports -- otherwise no column would let a reader check it.
    """
    token = text.strip()
    if not token:
        raise ValueError("--reference-ph must be a pH value")
    try:
        value = float(token)
    except ValueError as e:
        raise ValueError(f"--reference-ph: {token!r} is not a number") from e
    if not 0.0 <= value <= 14.0:
        raise ValueError(f"--reference-ph: {token!r} is outside pH 0-14")
    value = round(value, 2)
    tags = [ph_tag(v) for v in ph_values]
    if ph_tag(value) not in tags:
        raise ValueError(
            f"--reference-ph {token!r} is not in --ph-values "
            f"({','.join(num_token(v, True) for v in ph_values)}); the switch ratio "
            f"divides by F at the reference pH, so that pH must be one this run "
            f"reports -- add it to --ph-values or pick one of them"
        )
    return value


def derived_columns(
    ph_values: list[float], factors: list[float], reference_ph: float
) -> list[str]:
    """The pH-dependence columns of one run, in reading order.

    Names are ``<quantity>[_<pH or factor tag>][_<site set>]``: the site set (designed,
    all) always comes LAST. They depend on --ph-values / --release-factors, so the
    collector rebuilds this list from each row's own echo.
    """
    columns: list[str] = []
    reference_tag = ph_tag(reference_ph)
    for value in ph_values:
        tag = ph_tag(value)
        for site_set in SITE_SETS:
            columns.append(f"f_{tag}_{site_set}")
        for site_set in SITE_SETS:
            columns.append(f"f_ceiling_{tag}_{site_set}")
        # No ratio at the reference pH: it is 1 by construction and invites a filter.
        if tag != reference_tag:
            for site_set in SITE_SETS:
                columns.append(f"switch_{tag}_{site_set}")
    # Fold stability, free state only: dG_fold at every pH, losses against the reference.
    for value in ph_values:
        tag = ph_tag(value)
        columns.append(f"stab_dg_fold_{tag}")
        if tag != reference_tag:
            columns.append(f"stab_loss_{tag}")
            columns.append(f"stab_loss_{tag}_designed")
            columns.append(f"stab_loss_{tag}_other")
    for site_set in SITE_SETS:
        columns.append(f"ph_half_free_{site_set}")
        columns.append(f"ph_half_complex_{site_set}")
    for factor in factors:
        tag = factor_tag(factor)
        for site_set in SITE_SETS:
            columns.append(f"ph_release_{tag}_{site_set}")
        for site_set in SITE_SETS:
            columns.append(f"f_dips_below_{tag}_{site_set}")
    columns.append("ph_release_saturated")
    return columns


# The echo of the question and the facts established about the file: KEPT ON ERROR
# ROWS, where they diagnose a wrong chain letter or a wrong binder length.
TRUST_COLUMNS = [
    "propka_version",
    "binder_chains",
    "binder_chains_n_res",
    "binder_chains_n_het_res",
    "binder_chains_sequence",
    "expect_binder_len",
    "chain_ids_found",
    "centres_source",
    "free_only",
    "centres_n",
    "centres_reported",
    "centres_his_ok",
    "centres_pka_ok",
    "n_h_stripped",
    "n_h_stripped_guessed",
    "n_his_renamed",
    # Split per run: the complex run prechecks the whole assembly, so one truncated
    # residue in the target raises its count for every design. GATE ON _free.
    "n_propka_warnings_complex",
    "n_propka_warnings_free",
    # Complete before the centre loop (the likeliest error), so kept on error rows.
    "binder_term_groups_complex",
    "binder_term_groups_free",
    "term_groups_match",
    "free_state_model",
    # The pH question's echo, and the fact that decides whether the 'all' set exists.
    "ph_values",
    "reference_ph",
    "release_factors",
    "f_model",
    "n_binder_his_no_pka",
]

# The MEASUREMENT: every one is EMPTY on an error row.
SCIENCE_COLUMNS = [
    *[
        f"centre{k}_{field}"
        for k in range(1, MAX_CENTRE_COLUMNS + 1)
        for field in CENTRE_FIELDS
    ],
    "n_centres_coupled",
    "per_centre_pka_detail",
    "min_dpka",
    "n_binder_his",
    "binder_his_pka",
    # The non-pH half of the stability block: optional on an OK row, because a failure
    # there must not turn sound pKa measurements into an error row.
    "stab_n_groups_titrated",
    "stab_n_his_groups",
    "stab_profile_check_maxdiff",
    "stab_share_check",
    "stab_opt_ph",
    "stab_opt_dg",
    "stab_error",
    # The pH columns for the DEFAULT flags (other flags write other names; see
    # metric_columns_for).
    *derived_columns(
        list(DEFAULT_PH_VALUES), list(DEFAULT_RELEASE_FACTORS), DEFAULT_REFERENCE_PH
    ),
]

METRIC_COLUMNS = TRUST_COLUMNS + SCIENCE_COLUMNS
RESULT_COLUMNS = ["name", "status", "path", *METRIC_COLUMNS]


_DEFAULT_DERIVED = set(
    derived_columns(
        list(DEFAULT_PH_VALUES), list(DEFAULT_RELEASE_FACTORS), DEFAULT_REFERENCE_PH
    )
)


def metric_columns_for(
    ph_values: list[float], factors: list[float], reference_ph: float
) -> list[str]:
    """The metric columns a run with THESE flags writes, in reading order."""
    base = [c for c in SCIENCE_COLUMNS if c not in _DEFAULT_DERIVED]
    return TRUST_COLUMNS + base + derived_columns(ph_values, factors, reference_ph)


class StructureError(ValueError):
    """The file was read but holds nothing measurable, or could not be parsed."""


class PropkaError(ValueError):
    """PROPKA ran but produced nothing usable for this design."""


# ---------------------------------------------------------------------------
# Reading atom records
#
# Plain text, never a structure library: the image holds propka and nothing else, and
# PROPKA itself reads fixed-column PDB, so the file the tool writes must be the same
# bytes it read apart from the four documented edits.
# ---------------------------------------------------------------------------


class Record:
    """One copied atom line, plus the fields this tool indexes on."""

    __slots__ = ("line", "chain", "resname", "resnum", "icode", "is_atom")

    def __init__(
        self,
        line: str,
        chain: str,
        resname: str,
        resnum: int | None,
        icode: str,
        is_atom: bool,
    ) -> None:
        self.line = line
        self.chain = chain
        self.resname = resname
        self.resnum = resnum
        self.icode = icode
        self.is_atom = is_atom

    @property
    def key(self) -> tuple[str, str, int | None, str]:
        """Residue identity as the file gives it."""
        return (self.chain, self.resname, self.resnum, self.icode)


class Structure:
    """Every heavy-atom record of the first model, in file order."""

    def __init__(self) -> None:
        self.records: list[Record] = []
        self.chain_ids: list[str] = []
        # Chains that carry ATOM records: a chain of waters or ions is not a partner.
        self.atom_chain_ids: list[str] = []
        self.n_h_stripped = 0
        self.n_h_stripped_guessed = 0
        self.n_his_renamed = 0


def _is_hydrogen(name: str, element: str, is_atom: bool) -> tuple[bool, bool]:
    """``(is an H/D atom, the NAME fallback decided it)``.

    The element columns (77-78) win wherever they are filled. When they are blank the
    only guess available is the first alphabetic character of the atom name -- and
    that guess is made ON ATOM RECORDS ONLY. On a HETATM it would delete a mercury
    ('HG'), a holmium ('HO'), a hafnium ('HF'), a helium ('HE') or a sulphur-named
    ligand atom ('HS') from both PROPKA inputs, silently and with no column saying so.
    A HETATM with a blank element field is therefore KEPT whatever it is named.
    """
    if element:
        return element.upper() in ("H", "D"), False
    if not is_atom:
        return False, False
    for ch in name.strip():
        if ch.isalpha():
            return ch.upper() in ("H", "D"), ch.upper() in ("H", "D")
    return False, False


def _read_text(path: Path) -> str:
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return f.read()
    return path.read_text()


def read_structure(path: Path) -> Structure:
    """Parse a PDB file into heavy-atom records, applying edits 1, 2 and 4.

    CIF is NOT read here: PROPKA 3.5.1 refuses any extension but '.pdb'
    (propka/input.py: "Unknown input file type"), so the manifest builder converts a
    CIF to PDB with prosapia's ensure_pdb -- the same helper the pyrosetta tool uses --
    before the path ever reaches this worker.
    """
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    if path.name.lower().endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        raise StructureError(
            f"{path.name} is a CIF; PROPKA reads PDB only. The manifest builder "
            f"converts CIF inputs with ensure_pdb, so a CIF arriving here means the "
            f"worker was called by hand on one -- convert it first."
        )

    out = Structure()
    n_atoms = 0
    renamed: set[tuple[str, int | None, str, str]] = set()
    try:
        text = _read_text(path)
    except Exception as e:  # noqa: BLE001 - any read failure is "unreadable"
        raise StructureError(f"unreadable structure {path.name}: {e}") from e

    for line in text.splitlines():
        record = line[:6]
        if record == "ENDMDL" and n_atoms:
            break  # edit 4: first model only
        if record not in ("ATOM  ", "HETATM"):
            continue
        if len(line) < 54:
            raise StructureError(f"truncated {record.strip()} record: {line!r}")
        n_atoms += 1

        name = line[12:16].strip()
        element = line[76:78].strip() if len(line) >= 78 else ""
        is_h, by_name = _is_hydrogen(name, element, record == "ATOM  ")
        if is_h:
            out.n_h_stripped += 1  # edit 1
            if by_name:
                out.n_h_stripped_guessed += 1
            continue

        # PROPKA reads the chain from COLUMN 22 ONLY (line[21]). A two-character
        # chain ID spilling into column 21 would be invisible to it, so the whole
        # design is refused rather than measured against a chain PROPKA thinks is
        # blank.
        if line[20].strip():
            raise StructureError(
                f"{path.name} carries a chain ID in column 21 ({line[20:22]!r}); "
                f"PROPKA reads column 22 only, so a two-character chain ID would be "
                f"silently truncated"
            )
        chain = line[21].strip()
        resname = line[17:20].strip().upper()
        icode = line[26].strip()
        raw_num = line[22:26].strip()
        try:
            resnum: int | None = int(raw_num)
        except ValueError:
            # Hybrid-36 or a junk field. PROPKA can decode hybrid-36; this worker
            # cannot, so the residue stays in the file (PROPKA still sees it) but
            # carries no number here and can never be matched to a centre. A BINDER
            # record like that is refused in measure(), where the binder chain is
            # known.
            resnum = None

        written = line
        if resname in HIS_NAMES and resname != "HIS":
            written = f"{line[:17]}HIS{line[20:]}"  # edit 2
            # RESIDUES, not atom records: the same residue contributes many lines.
            renamed.add((chain, resnum, icode, resname))

        if chain not in out.chain_ids:
            out.chain_ids.append(chain)
        if record == "ATOM  " and chain not in out.atom_chain_ids:
            out.atom_chain_ids.append(chain)
        out.records.append(
            Record(written, chain, "HIS" if resname in HIS_NAMES else resname,
                   resnum, icode, record == "ATOM  ")
        )

    out.n_his_renamed = len(renamed)

    if n_atoms == 0:
        raise StructureError(f"no ATOM/HETATM records in {path.name}")
    if not out.records:
        raise StructureError(f"{path.name} holds hydrogens only, no heavy atoms")
    return out


def write_pdb(records: list[Record], dst: Path) -> None:
    """Write the records out, with a TER at every chain change (edit 3) and an END."""
    lines: list[str] = []
    previous: str | None = None
    for rec in records:
        if previous is not None and rec.chain != previous:
            lines.append("TER   ")
        lines.append(rec.line)
        previous = rec.chain
    if lines:
        lines.append("TER   ")
    lines.append("END   ")
    dst.write_text("\n".join(lines) + "\n")


def residues_of(
    records: list[Record], chains: frozenset[str], atom_only: bool
) -> list[Record]:
    """The FIRST record of each residue of ``chains``, in file order."""
    out: list[Record] = []
    last: tuple[str, str, int | None, str] | None = None
    for rec in records:
        if rec.chain not in chains:
            continue
        if atom_only and not rec.is_atom:
            continue
        if rec.key != last:
            last = rec.key
            out.append(rec)
    return out


# ---------------------------------------------------------------------------
# PROPKA
# ---------------------------------------------------------------------------


class _WarningCounter(logging.Handler):
    """Counts WARNING-and-above records PROPKA emits during one run, AND PRINTS THEM.

    Counting alone would be a dead end: ``propka.run.single`` installs no handler of
    its own, so without this the text of every warning would be lost and the
    n_propka_warnings_* columns would send a reader to a task .err that says nothing.
    Each line carries the design and the run it belongs to, because one task holds
    several designs and each design has two runs.
    """

    def __init__(self, tag: str) -> None:
        super().__init__(level=logging.WARNING)
        self.tag = tag
        self.n = 0

    def emit(self, record: logging.LogRecord) -> None:
        if record.levelno < logging.WARNING:
            return
        self.n += 1
        try:
            message = record.getMessage()
        except Exception:  # noqa: BLE001 - a broken format string must not kill a task
            message = repr(record.msg)
        print(f"{self.tag} propka {record.levelname}: {message}", file=sys.stderr)


class GroupResult:
    """One PROPKA histidine group, as this tool needs it."""

    __slots__ = (
        "pka",
        "model_pka",
        "buried_pct",
        "desolv_reg",
        "desolv_re",
        "sc_hb",
        "bb_hb",
        "coulomb",
        "n_coupled",
        "reason",
    )

    def __init__(self) -> None:
        self.pka: float | None = None
        # How many groups PROPKA found NON-COVALENTLY coupled to this one -- the
        # groups it marks with '*' in its own summary. Non-zero means PROPKA is
        # telling you its pairwise-independent model is breaking down for this group,
        # which for a dpka (a DIFFERENCE of two such pKa values) is the warning that
        # matters most. 0 is the ordinary case.
        self.n_coupled: int = 0
        self.model_pka: float | None = None
        self.buried_pct: float | None = None
        self.desolv_reg: float | None = None
        self.desolv_re: float | None = None
        self.sc_hb: float | None = None
        self.bb_hb: float | None = None
        self.coulomb: float | None = None
        # Why there is no pKa, when there is none. Never an empty string on a miss:
        # the caller puts it in the error message.
        self.reason: str = ""


def propka_version() -> str:
    import propka

    return str(getattr(propka, "__version__", "unknown"))


def run_propka(
    pdb_path: Path, label: str, design: str, binder_chains: frozenset[str]
) -> tuple[dict[tuple[str, int], GroupResult], list[str], int, Any]:
    """Run PROPKA3 on one file.

    Returns its HISTIDINE groups by (chain, resnum), the terminal groups of the binder
    chains, the number of warnings it emitted, and the MolecularContainer (which the
    stab_* block reads from the free run, so propka is not re-run).

    Values come from the Python API (propka.run.single) and its conformation-average
    container 'AVR', the object PROPKA's own .pka file is written from; nothing is
    parsed out of PROPKA's text output.
    """
    import propka.run

    logger = logging.getLogger("propka")
    counter = _WarningCounter(f"{design} [{label}]")
    previous_level = logger.level
    logger.addHandler(counter)
    # PROPKA logs its whole determinant and summary table at INFO on every run. The
    # tool reads the objects instead, so the text is pure noise in a task log.
    logger.setLevel(logging.WARNING)
    try:
        mol = propka.run.single(str(pdb_path), optargs=("-q",), write_pka=False)
    except Exception as e:  # noqa: BLE001 - recorded as data by the caller
        raise PropkaError(f"propka failed on the {label}: {e}") from e
    finally:
        logger.removeHandler(counter)
        logger.setLevel(previous_level)

    conformation = getattr(mol, "conformations", {}).get("AVR")
    if conformation is None:
        raise PropkaError(
            f"propka produced no 'AVR' conformation for the {label}; the file held no "
            f"usable molecular conformation"
        )

    groups: dict[tuple[str, int], GroupResult] = {}
    term_groups: list[str] = []
    for group in conformation.groups:
        residue_type = str(getattr(group, "residue_type", "")).strip()
        if residue_type in ("N+", "C-"):
            # Terminal groups PROPKA built on the binder chains. N+ depends on TER/OXT/
            # MODEL placement, so it can differ between runs: term_groups_match.
            atom = group.atom
            chain_id = str(getattr(atom, "chain_id", "")).strip()
            if chain_id in binder_chains:
                try:
                    number = int(atom.res_num)
                except (TypeError, ValueError):
                    # A trust extra: an unreadable number must not fail the row, and is
                    # kept RAW so two different residues cannot compare equal.
                    term_groups.append(
                        f"{residue_type}<unreadable {atom.res_num!r}>{chain_id}"
                    )
                else:
                    term_groups.append(f"{residue_type}{number}{chain_id}")
            continue
        if not _is_his_group(group):
            continue
        atom = group.atom
        try:
            key = (str(atom.chain_id).strip(), int(atom.res_num))
        except (TypeError, ValueError):
            continue
        result = GroupResult()
        if key in groups:
            # Two HIS groups at one (chain, resnum): insertion codes, which PROPKA
            # drops from its residue number. Which one a pKa belongs to is unknowable,
            # so neither is reported.
            groups[key] = _ambiguous(
                f"two histidine groups share chain {key[0]} residue {key[1]} in the "
                f"{label} (insertion codes?), so neither pKa can be attributed"
            )
            continue
        # COVALENT coupling only. PROPKA sets coupled_titrating_group in
        # find_covalently_coupled_groups, which penalises all but one member of a
        # covalently coupled system and then omits them from its own summary. For an
        # ordinary protein histidine this branch is close to dead code -- the
        # mechanism that actually fires is NON-covalent coupling, which does not
        # suppress a pKa and is reported as n_coupled below.
        coupled = getattr(group, "coupled_titrating_group", None)
        if coupled is not None:
            result.reason = (
                f"propka discarded this histidine in the {label}: it is COVALENTLY "
                f"coupled with {str(getattr(coupled, 'label', '?')).strip()} and was "
                f"penalised, so propka reports no pKa for it"
            )
            groups[key] = result
            continue
        try:
            pka = float(group.pka_value)
        except (TypeError, ValueError):
            pka = float("nan")
        if pka != pka or pka in (float("inf"), float("-inf")):  # NaN / inf
            result.reason = f"propka returned no finite pKa in the {label}"
            groups[key] = result
            continue
        result.pka = pka
        # NON-covalent coupling (what PROPKA marks with '*'): it does NOT suppress the
        # pKa, it says the pKa is less trustworthy; two designed histidines near each
        # other are the canonical trigger. A missing attribute must not read as "not
        # coupled" (a silent 0 on a gate column), so it is an error: re-check the names
        # 'non_covalently_coupled_groups' and 'residue_type' before moving the pin.
        if not hasattr(group, "non_covalently_coupled_groups"):
            raise PropkaError(
                f"propka's histidine group has no 'non_covalently_coupled_groups' in the "
                f"{label} (a propka release other than the pinned 3.5.1?): the coupling "
                f"counts would read 0 where the honest answer is unknown"
            )
        result.n_coupled = len(group.non_covalently_coupled_groups or ())
        result.model_pka = _maybe_float(getattr(group, "model_pka", None))
        buried = _maybe_float(getattr(group, "buried", None))
        result.buried_pct = None if buried is None else 100.0 * buried
        result.desolv_reg = _maybe_float(getattr(group, "energy_volume", None))
        result.desolv_re = _maybe_float(getattr(group, "energy_local", None))
        determinants = getattr(group, "determinants", {}) or {}
        result.sc_hb = _determinant_sum(determinants, "sidechain")
        result.bb_hb = _determinant_sum(determinants, "backbone")
        result.coulomb = _determinant_sum(determinants, "coulomb")
        groups[key] = result

    return groups, sorted(term_groups), counter.n, mol


def _is_his_group(group: Any) -> bool:
    """Is this propka group a HISTIDINE group?

    THE ONE RULE, used by BOTH the histidine dictionary of ``run_propka`` and the
    designed set of ``folding_stability``, so the two cannot drift. Deliberately an
    exact match on the stripped ``residue_type`` and NOT case-folded: that is what
    ``run_propka`` has always done, and a second, laxer rule here would let
    ``folding_stability`` call a group designed that the pKa columns never saw.
    """
    return str(getattr(group, "residue_type", "")).strip() == "HIS"


def _group_key(group: Any) -> tuple[str, int] | None:
    """``(chain, resnum)`` of a propka group, or None when either is unreadable.

    The same ``(chain, resnum)`` key the histidine dictionary of ``run_propka`` is
    built with -- but the key ALONE does not identify a histidine: PROPKA also builds
    an N+ group on residue 1's N atom, a C- group on the last residue, and ion/ligand
    groups, any of which can share a centre's (chain, resnum). Callers that mean "the
    designed histidine" must pair this with ``_is_his_group``, exactly as
    ``run_propka`` does.
    """
    atom = getattr(group, "atom", None)
    if atom is None:
        return None
    try:
        return (str(atom.chain_id).strip(), int(atom.res_num))
    except (TypeError, ValueError, AttributeError):
        return None


def folding_stability(
    mol: Any,
    designed_keys: set[tuple[str, int]],
    ph_values: list[float],
    reference_ph: float,
) -> dict[str, Any]:
    """The stab_* block for ONE free-binder run. Read the SIGN note at the top.

    ``mol`` is the MolecularContainer of the FREE run (binder chain only, same
    coordinates). Everything here is one extra pass over its 'AVR' groups -- propka is
    not re-run.

    THE PARTITION. ``total`` is propka's own
    ``ConformationContainer.calculate_folding_energy``, which sums
    ``Group.calculate_folding_energy`` over EVERY group (non-titratable ones return
    0.00). ``his`` is that same per-group term summed over the DESIGNED CENTRES only,
    and ``other`` is summed INDEPENDENTLY over every remaining group -- not as
    ``total - his``, which would make ``stab_share_check`` identically zero and
    prove nothing. The check then verifies that the two parts really do reconstruct
    propka's own total at every reported pH.

    WHAT THE CHECK DOES NOT CATCH: if ``designed_keys`` matched no group at all, ``his``
    is 0, ``other`` is the total, and the check still passes. ``stab_n_his_groups`` is
    the column that catches that -- it must equal ``centres_n``. (On this code path the
    centre loop has already raised if the free run produced no histidine group for a
    centre, so a mismatch here means the key convention drifted, not the chemistry.)
    """
    conformation = mol.conformations["AVR"]
    parameters = conformation.parameters
    groups = list(conformation.groups)
    # HISTIDINE GROUPS ONLY, and then by key. The key alone is not enough: PROPKA
    # builds an N+ group on residue 1's N atom and a C- group on the last residue, and
    # either can carry a centre's (chain, resnum) -- on a binder whose centre happens
    # to be residue 1 or the last residue, a key-only match would move the TERMINUS's
    # folding term into the designed share and report it as the histidine's.
    designed = [
        g for g in groups if _is_his_group(g) and _group_key(g) in designed_keys
    ]
    if len(designed) != len(designed_keys):
        raise PropkaError(
            f"{len(designed)} histidine group(s) of the free run match the "
            f"{len(designed_keys)} designed centre(s): the designed share would read 0 "
            f"and the 'other' share absorb the whole loss"
        )
    # The complement BY OBJECT IDENTITY, never by re-testing the predicate: that is
    # what makes `designed` and `others` an exact partition of `groups` whatever the
    # predicate does, so stab_share_check stays a real reconstruction of propka's
    # own total instead of a tautology that two overlapping or two leaky filters would
    # both pass.
    _chosen = {id(g) for g in designed}
    others = [g for g in groups if id(g) not in _chosen]

    out: dict[str, Any] = {
        "stab_n_groups_titrated": sum(
            1 for g in groups if bool(getattr(g, "titratable", False))
        ),
        "stab_n_his_groups": len(designed),
    }

    def total_at(ph: float) -> float:
        return float(
            conformation.calculate_folding_energy(ph=ph, reference=STAB_REFERENCE)
        )

    def part_at(subset: list[Any], ph: float) -> float:
        return float(
            sum(
                g.calculate_folding_energy(parameters, ph=ph, reference=STAB_REFERENCE)
                for g in subset
            )
        )

    reference_dg = total_at(reference_ph)
    reference_his = part_at(designed, reference_ph)
    reference_other = part_at(others, reference_ph)
    reference_tag = ph_tag(reference_ph)

    share_diffs = [abs(reference_his + reference_other - reference_dg)]
    for ph in ph_values:
        tag = ph_tag(ph)
        total = total_at(ph)
        his = part_at(designed, ph)
        other = part_at(others, ph)
        share_diffs.append(abs(his + other - total))
        out[f"stab_dg_fold_{tag}"] = total
        if tag == reference_tag:
            continue
        # POSITIVE = LESS STABLE AT THE LOWER pH. See the SIGN note at the top of the
        # module: propka's optimum is the MINIMUM of this function.
        out[f"stab_loss_{tag}"] = total - reference_dg
        out[f"stab_loss_{tag}_designed"] = his - reference_his
        out[f"stab_loss_{tag}_other"] = other - reference_other
    out["stab_share_check"] = max(share_diffs)

    # propka's own profile: the optimum, and the proof that the direct evaluations
    # above call the same function as the "Free energy of folding" block of a .pka
    # file. The grid is passed explicitly so neither depends on propka's options.
    profile, opt, _range80, _stable = mol.get_folding_profile(
        conformation="AVR", reference=STAB_REFERENCE, grid=STAB_PROFILE_GRID
    )
    opt_ph, opt_dg = (opt + (None, None))[:2] if isinstance(opt, tuple) else (None, None)
    out["stab_opt_ph"] = None if opt_ph is None else float(opt_ph)
    # dg is only meaningful together with the pH it belongs to: if propka could not
    # determine the optimum pH, its 1e6 sentinel must not be published as an energy.
    out["stab_opt_dg"] = None if opt_ph is None or opt_dg is None else float(opt_dg)

    maxdiff = None
    if profile:
        diffs: list[float] = []
        for target in STAB_CHECK_PHS:
            nearest = min(profile, key=lambda point: abs(float(point[0]) - target))
            if abs(float(nearest[0]) - target) > STAB_GRID_TOL:
                continue
            diffs.append(abs(float(nearest[1]) - total_at(target)))
        maxdiff = max(diffs) if diffs else None
    out["stab_profile_check_maxdiff"] = maxdiff
    return out


def _ambiguous(reason: str) -> GroupResult:
    result = GroupResult()
    result.reason = reason
    return result


def _maybe_float(value) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return None if out != out else out


def _determinant_sum(determinants: dict, kind: str) -> float:
    """The total of PROPKA's determinants of one kind, in pK units."""
    total = 0.0
    for determinant in determinants.get(kind, []) or []:
        value = _maybe_float(getattr(determinant, "value", None))
        if value is not None:
            total += value
    return total


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------


def _pow10(exponent: float) -> float:
    """``10 ** exponent``, clamped so a pathological pKa cannot raise OverflowError.

    Real pKa values sit in roughly -5..25 and pH in 0..14, so the clamp is never
    reached by anything this tool measures; it exists so that a corrupt number ends
    as a huge float a reader can see rather than as a crashed task.
    """
    if exponent > 300.0:
        return float("1e300")
    if exponent < -300.0:
        return 0.0
    return 10.0**exponent


def f_factor(pka_complex: float, pka_free: float, ph: float) -> float:
    """F_i(pH) for one site: the factor by which binding is weakened at this pH."""
    return (1.0 + _pow10(pka_free - ph)) / (1.0 + _pow10(pka_complex - ph))


def f_product(sites: list[tuple[float, float]], ph: float) -> float:
    """F(pH) = prod_i F_i, over (pka_complex, pka_free) pairs. Empty set -> 1.0."""
    total = 1.0
    for pka_complex, pka_free in sites:
        total *= f_factor(pka_complex, pka_free, ph)
    return total


def f_ceiling(sites: list[tuple[float, float]], ph: float) -> float:
    """prod_i (1 + 10^(pKa_free_i - pH)): the bound F cannot exceed at this pH,
    whatever the complex pKa values are. Set by the FREE binder alone."""
    total = 1.0
    for _pka_complex, pka_free in sites:
        total *= 1.0 + _pow10(pka_free - ph)
    return total


def _protonated_sum(pkas: list[float], ph: float) -> float:
    """sum_i 1/(1 + 10^(pH - pKa_i)): how many of the sites are protonated at pH."""
    return sum(1.0 / (1.0 + _pow10(ph - pka)) for pka in pkas)


def titration_midpoint(pkas: list[float]) -> float | None:
    """The pH at which HALF the sites are protonated: sum_i theta_i(pH) = n/2.

    For ONE site this is exactly its pKa. For several it is a single SUMMARY number
    that hides the individual values -- two sites at 5.0 and 8.0 give the same
    midpoint as two sites at 6.5. The summed fraction is strictly decreasing in pH,
    so bisection on [0, 14] is sound; None when the half-point lies outside that
    window (every pKa above 14 or below 0), which must never be read as a pH of 0.
    """
    if not pkas:
        return None
    target = len(pkas) / 2.0

    def g(ph: float) -> float:
        return _protonated_sum(pkas, ph) - target

    lo, hi = MIDPOINT_PH_LO, MIDPOINT_PH_HI
    if g(lo) < 0.0 or g(hi) > 0.0:
        return None
    while hi - lo > MIDPOINT_TOL:
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def release_ph(sites: list[tuple[float, float]], factor: float) -> tuple[float | None, bool]:
    """``(highest pH in [3, 9] where F >= factor, already saturated at pH 9)``.

    Scanned DOWNWARDS on a 0.01 grid and NOT solved as a root: F is NOT guaranteed
    monotonic in pH, so the first hit from the top is the only well-defined answer.

    READ IT AS "the highest pH at which F reaches X", AND NOTHING MORE. It does NOT
    mean F >= X everywhere below that pH. A site whose complex pKa is ABOVE its free
    pKa (a positive dpka) contributes a factor below 1 that keeps falling as the pH
    drops, so a mixed binder can cross X on the way down and cross back: sites
    (free 6.5, complex 5.0) and (free 4.0, complex 6.0) peak near pH 5.5 and have
    collapsed well below X by pH 4. dips_below() reports exactly that case per row,
    and f_*_ph<tag> at the pH of interest is the number to gate on.

    Two different NA cases, which the caller must keep apart:
      * (None, False) -- F never reaches the factor anywhere in [3, 9].
      * (None, True)  -- F is ALREADY at or above the factor at pH 9.0, so the answer
                         lies above the scanned window. Reporting 9.0 would read as a
                         release pH of 9, which is a different and wrong claim.
    """
    if not sites:
        return None, False
    if f_product(sites, RELEASE_PH_HI) >= factor:
        return None, True
    steps = int(round((RELEASE_PH_HI - RELEASE_PH_LO) / RELEASE_PH_STEP))
    for step in range(steps + 1):
        # Built from an integer step so the grid cannot drift: 9.00, 8.99, ... 3.00.
        ph = round(RELEASE_PH_HI - step * RELEASE_PH_STEP, 2)
        if f_product(sites, ph) >= factor:
            return ph, False
    return None, False


def dips_below(sites: list[tuple[float, float]], factor: float, start_ph: float) -> bool:
    """Does F fall below ``factor`` anywhere between ``start_ph`` and pH 3.0?

    The test of whether "weakened at least X-fold at and below this pH" -- the reading
    a scientist will reach for -- is TRUE for this row. Same 0.01 grid as release_ph,
    started at the release pH it found. True means the reading is false for this row:
    F recrosses X on the way down, so the complex is NOT weakened X-fold at every pH
    below the release pH.
    """
    steps = int(round((start_ph - RELEASE_PH_LO) / RELEASE_PH_STEP))
    for step in range(steps + 1):
        ph = round(start_ph - step * RELEASE_PH_STEP, 2)
        if f_product(sites, ph) < factor:
            return True
    return False


def _fmt_pka(value: float | None) -> str:
    """A pKa for a text cell. 'NA' is an ABSENCE; 0.00 would be a measurement."""
    return "NA" if value is None else f"{value:.2f}"


def _parse_expect(raw: str, column: str) -> int | None:
    """An expected binder length from a table cell: an integer, or a sequence's length.

    A sequence may carry '/' chain breaks (not residues). Anything else raises: a
    silently ignored check is worse than no check.
    """
    text = raw.strip()
    if not text:
        return None
    if text == UNRESOLVED_EXPECT:
        raise ValueError(
            f"the binder-length check was asked for but {column} resolves nowhere for "
            f"this design (neither on its own row nor on any ancestor); a check that "
            f"cannot be made is an error, not a skipped check"
        )
    try:
        return int(text)
    except ValueError:
        pass
    # An integer column round-trips as '60.0' once any row of it is NA (pandas floats it).
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
    seq = text.replace(" ", "").replace("/", "").upper()
    if seq.isalpha():
        return len(seq)
    raise ValueError(
        f"{column} holds {raw!r}, which is neither an integer length nor a plain "
        f"one-letter sequence; the binder-length check cannot be made"
    )


def _as_resnum(text: str, position: int, source: str) -> int:
    """A residue number from a token: an integer, or an integral float string."""
    try:
        return int(text)
    except ValueError:
        pass
    try:
        value = float(text)
    except ValueError as e:
        raise ValueError(
            f"centre {position}: {text!r} is not a residue number (from {source})"
        ) from e
    if not value.is_integer():
        raise ValueError(
            f"centre {position}: {text!r} is not a whole residue number (from {source})"
        )
    return int(value)


def parse_centres(
    cell: str, chains: list[str], source: str
) -> list[tuple[str, int]] | None:
    """A centre spec -> ``[(chain, resnum)]``, or ``None`` for 'auto' (every histidine).

    Tokens are ``<chain>:<resnum>`` or a bare ``<resnum>``, the latter only when exactly
    one binder chain is analysed. An empty token is refused. ``source`` says where the
    spec came from, for messages.
    """
    text = cell.strip()
    if text.lower() == AUTO_CENTRES:
        return None
    if not text:
        raise ValueError(
            f"no centre residue number could be resolved from {source} for this "
            f"design (the spec resolved to nothing); there is nothing to measure"
        )
    values: list[tuple[str, int]] = []
    for position, raw_token in enumerate(text.split(","), start=1):
        token = raw_token.strip()
        if not token:
            raise ValueError(
                f"{source} has an empty slot (slot {position}): remove the stray comma"
            )
        chain_text, sep, number_text = token.rpartition(":")
        if sep:
            chain = chain_text.strip()
            if chain not in chains:
                raise ValueError(
                    f"centre {position}: chain {chain!r} in {token!r} is not one of "
                    f"the binder chains {','.join(chains)} (from {source})"
                )
        elif len(chains) == 1:
            chain = chains[0]
        else:
            raise ValueError(
                f"centre {position}: {token!r} names no chain but {len(chains)} binder "
                f"chains ({','.join(chains)}) are analysed; write <chain>:<resnum> "
                f"(from {source})"
            )
        values.append((chain, _as_resnum(number_text.strip(), position, source)))
    if len(set(values)) != len(values):
        raise ValueError(f"{source} names the same centre twice: {values}")
    return values


def measure(
    src: Path,
    binder_chains: list[str],
    centres_cell: str,
    centres_source: str,
    free_only: bool,
    expect_len: int | None,
    expect_col: str,
    keep_dir: Path | None,
    metrics: dict,
    ph_values: list[float],
    factors: list[float],
    reference_ph: float,
    design: str = "design",
) -> None:
    """Fill ``metrics`` for one design. Raises on anything that makes a number a lie."""
    chain_set = frozenset(binder_chains)
    structure = read_structure(src)
    metrics["chain_ids_found"] = ",".join(structure.chain_ids)
    metrics["n_h_stripped"] = structure.n_h_stripped
    metrics["n_h_stripped_guessed"] = structure.n_h_stripped_guessed
    metrics["n_his_renamed"] = structure.n_his_renamed

    absent = [c for c in binder_chains if c not in structure.chain_ids]
    if absent:
        raise StructureError(
            f"binder chain(s) {','.join(absent)} not in {src.name}; chains found: "
            f"{','.join(structure.chain_ids) or '(none)'}"
        )
    partners = [c for c in structure.atom_chain_ids if c not in chain_set]
    if not free_only and not partners:
        raise StructureError(
            f"{src.name} holds no chain outside the binder chains "
            f"{','.join(binder_chains)}: the complex and the free state would be the "
            f"same molecule and every dpka 0 by construction. Pass --free-only to run "
            f"the free state alone."
        )

    binder_res = residues_of(structure.records, chain_set, atom_only=True)
    binder_het = residues_of(structure.records, chain_set, atom_only=False)
    metrics["binder_chains_n_res"] = len(binder_res)
    metrics["binder_chains_n_het_res"] = len(binder_het) - len(binder_res)
    # The only post-hoc evidence of WHICH molecule the chain letters named, whether or
    # not a length was expected ('/' between chains).
    metrics["binder_chains_sequence"] = "/".join(
        "".join(THREE_TO_ONE.get(r.resname, "X") for r in binder_res if r.chain == c)
        for c in binder_chains
    )
    if not binder_res:
        raise StructureError(
            f"chain(s) {','.join(binder_chains)} hold no ATOM residues in {src.name} "
            f"(only HETATM?)"
        )
    unnumbered = [r for r in binder_res if r.resnum is None]
    if unnumbered:
        raise StructureError(
            f"{len(unnumbered)} residue(s) of the binder chains in {src.name} have an "
            f"unreadable residue number (hybrid-36 or junk in columns 23-26); a centre "
            f"could not be matched reliably"
        )

    if expect_len is not None:
        metrics["expect_binder_len"] = expect_len
        if len(binder_res) != expect_len:
            raise StructureError(
                f"binder chains {','.join(binder_chains)} have {len(binder_res)} "
                f"residues but {expect_col} says {expect_len}; the chain letters may "
                f"name the wrong molecule, and every pKa below would be the wrong "
                f"residue's"
            )

    # --- the His check, before any pKa is computed -------------------------------
    his_all = [r for r in binder_res if r.resname == "HIS"]
    named = parse_centres(centres_cell, binder_chains, centres_source)
    if named is None:
        if not his_all:
            raise StructureError(
                f"centres 'auto': chain(s) {','.join(binder_chains)} of {src.name} "
                f"hold no histidine"
            )
        named = [(r.chain, int(r.resnum)) for r in his_all]
    metrics["centres_n"] = len(named)
    metrics["centres_reported"] = min(len(named), MAX_CENTRE_COLUMNS)
    if len(named) > MAX_CENTRE_COLUMNS:
        print(
            f"WARNING: {design}: {len(named)} centres, only the first "
            f"{MAX_CENTRE_COLUMNS} get centre<k>_* columns (centres_reported); every "
            f"other column covers all of them",
            file=sys.stderr,
        )
    by_key: dict[tuple[str, int], list[Record]] = {}
    for residue in binder_res:
        by_key.setdefault((residue.chain, int(residue.resnum)), []).append(residue)

    centres: list[Record] = []
    his_ok = 0
    for k, (chain, number) in enumerate(named, start=1):
        found = by_key.get((chain, number), [])
        if not found:
            raise StructureError(
                f"centre {k}: chain {chain} of {src.name} has no residue numbered "
                f"{number} (from {centres_source})"
            )
        if len(found) > 1:
            raise StructureError(
                f"centre {k}: residue number {number} occurs {len(found)} times in "
                f"chain {chain} of {src.name} (insertion codes or repeated numbering); "
                f"which residue the centre means is ambiguous"
            )
        residue = found[0]
        if residue.resname != "HIS":
            raise StructureError(
                f"centre {k}: residue {number} of chain {chain} in {src.name} is "
                f"{residue.resname}, not a histidine; the number names the wrong "
                f"residue"
            )
        centres.append(residue)
        his_ok += 1
        metrics["centres_his_ok"] = his_ok
    metrics["centres_his_ok"] = his_ok

    # --- the PROPKA runs ---------------------------------------------------------
    # Recorded before the runs: a row without propka_version failed on the file.
    metrics["propka_version"] = propka_version()
    complex_groups: dict[tuple[str, int], GroupResult] = {}
    complex_terms: list[str] = []
    with tempfile.TemporaryDirectory(prefix="propka_") as tmp:
        tmp_dir = Path(tmp)
        complex_pdb = tmp_dir / "complex.pdb"
        free_pdb = tmp_dir / "free.pdb"
        write_pdb([r for r in structure.records if r.chain in chain_set], free_pdb)
        if not free_only:
            write_pdb(structure.records, complex_pdb)
        if keep_dir is not None:
            keep_dir.mkdir(parents=True, exist_ok=True)
            (keep_dir / "free.pdb").write_text(free_pdb.read_text())
            if not free_only:
                (keep_dir / "complex.pdb").write_text(complex_pdb.read_text())

        # Each count is recorded as its run returns, so a failure of the second run
        # still leaves the first run's count on the error row.
        if not free_only:
            complex_groups, complex_terms, n_warn_complex, _complex_mol = run_propka(
                complex_pdb, "complex", design, chain_set
            )
            metrics["n_propka_warnings_complex"] = n_warn_complex
        free_groups, free_terms, n_warn_free, free_mol = run_propka(
            free_pdb, "free", design, chain_set
        )
        metrics["n_propka_warnings_free"] = n_warn_free

    # The binder's terminal groups per run SHOULD be identical (C- from an OXT atom,
    # N+ from the first nitrogen); a difference is a file artefact that would put a
    # spurious shift on an N-terminal-proximal centre. Written before anything else
    # that can raise, so they survive on the likeliest error row.
    metrics["binder_term_groups_free"] = ",".join(free_terms)
    if not free_only:
        metrics["binder_term_groups_complex"] = ",".join(complex_terms)
        metrics["term_groups_match"] = 1 if complex_terms == free_terms else 0

    # --- every histidine of the binder chains ------------------------------------
    metrics["n_binder_his"] = len(his_all)
    entries: list[str] = []
    # The 'all' site set at FULL PRECISION, from the group objects, never re-read from
    # the two-decimal binder_his_pka text. Under --free-only the complex slot is NaN
    # (the free-only columns never read it).
    all_sites: list[tuple[float, float]] = []
    n_no_pka = 0
    for residue in his_all:
        number = int(residue.resnum)
        in_complex = complex_groups.get((residue.chain, number))
        in_free = free_groups.get((residue.chain, number))
        pka_complex = None if in_complex is None else in_complex.pka
        pka_free = None if in_free is None else in_free.pka
        entries.append(
            f"{residue.chain}:{number}:{_fmt_pka(pka_complex)}:{_fmt_pka(pka_free)}"
        )
        if pka_free is None or (not free_only and pka_complex is None):
            # One histidine without a pKa makes the product over 'all' undefined: the
            # product over the rest would be a different quantity under the same name.
            n_no_pka += 1
            continue
        all_sites.append((float("nan") if free_only else pka_complex, pka_free))
    metrics["binder_his_pka"] = ";".join(entries)
    metrics["n_binder_his_no_pka"] = n_no_pka

    # --- the designed centres ----------------------------------------------------
    details: list[str] = []
    dpkas: list[float] = []
    designed_sites: list[tuple[float, float]] = []
    pka_ok = 0
    n_coupled_centres = 0
    for k, residue in enumerate(centres, start=1):
        number = int(residue.resnum)
        in_complex = complex_groups.get((residue.chain, number))
        in_free = free_groups.get((residue.chain, number))
        for state, group in (("complex", in_complex), ("free binder", in_free)):
            if state == "complex" and free_only:
                continue
            if group is None:
                raise PropkaError(
                    f"centre {k}: propka produced no histidine group for residue "
                    f"{number} of chain {residue.chain} in the {state} run -- the "
                    f"residue is a HIS in the file, so its side chain is probably "
                    f"incomplete"
                )
            if group.pka is None:
                raise PropkaError(f"centre {k} (residue {number}): {group.reason}")
        assert in_free is not None and in_free.pka is not None  # for type checkers
        # Full precision, straight off the group objects: the F columns must not be a
        # recomputation from the rounded per-centre cells.
        if free_only:
            designed_sites.append((float("nan"), in_free.pka))
        else:
            assert in_complex is not None and in_complex.pka is not None
            dpkas.append(in_complex.pka - in_free.pka)
            designed_sites.append((in_complex.pka, in_free.pka))
        pka_ok += 1
        # Written inside the loop, so an error on centre 2 still records centre 1's.
        metrics["centres_pka_ok"] = pka_ok
        if k <= MAX_CENTRE_COLUMNS:
            metrics[f"centre{k}_chain"] = residue.chain
            metrics[f"centre{k}_resnum"] = number
            metrics[f"centre{k}_pka_free"] = in_free.pka
            metrics[f"centre{k}_buried_free"] = in_free.buried_pct
            metrics[f"centre{k}_n_coupled_free"] = in_free.n_coupled
            if not free_only:
                metrics[f"centre{k}_pka_complex"] = in_complex.pka
                metrics[f"centre{k}_dpka"] = dpkas[-1]
                metrics[f"centre{k}_buried_complex"] = in_complex.buried_pct
                metrics[f"centre{k}_n_coupled_complex"] = in_complex.n_coupled
        if not free_only:
            if in_complex.n_coupled or in_free.n_coupled:
                n_coupled_centres += 1
            details.append(
                f"{residue.chain}:{number}:model={_fmt_pka(in_complex.model_pka)}"
                f":buried={_fmt_pka(in_complex.buried_pct)}"
                f":desolv_reg={_fmt_pka(in_complex.desolv_reg)}"
                f":desolv_re={_fmt_pka(in_complex.desolv_re)}"
                f":sc_hb={_fmt_pka(in_complex.sc_hb)}"
                f":bb_hb={_fmt_pka(in_complex.bb_hb)}"
                f":coulomb={_fmt_pka(in_complex.coulomb)}"
            )
    metrics["centres_pka_ok"] = pka_ok
    metrics["free_state_model"] = FREE_STATE_MODEL
    if not free_only:
        # Centres PROPKA found non-covalently coupled in EITHER run: one clause to filter.
        metrics["n_centres_coupled"] = n_coupled_centres
        metrics["per_centre_pka_detail"] = ";".join(details)
        metrics["min_dpka"] = min(dpkas)

    # --- fold stability vs pH, FREE STATE ONLY (no extra PROPKA run) ---------------
    # Caught, not raised: these columns are an addition, and a failure here must leave
    # stab_* NA and the row OK. stab_error says why: an NA there is not a zero.
    try:
        metrics.update(
            folding_stability(
                free_mol,
                {(r.chain, int(r.resnum)) for r in centres},
                ph_values,
                reference_ph,
            )
        )
        metrics["stab_error"] = ""
    except Exception as e:  # noqa: BLE001 - recorded as data, never fatal here
        metrics["stab_error"] = " ".join(str(e).split())[:200]

    # --- the pH dependence (arithmetic only) -------------------------------------
    # DESIGNED is always available (a centre without both pKa values has raised). ALL
    # is every binder histidine and is refused outright when any lacks a pKa, so an
    # 'all' column is never a product over a silently smaller set. Under --free-only
    # only what needs the free pKa alone is written: f_ceiling_* and ph_half_free_*.
    site_sets: dict[str, list[tuple[float, float]] | None] = {
        "designed": designed_sites,
        "all": all_sites if n_no_pka == 0 else None,
    }
    reference_tag = ph_tag(reference_ph)
    for ph in ph_values:
        tag = ph_tag(ph)
        for name, sites in site_sets.items():
            if sites is None:
                continue
            if not free_only:
                metrics[f"f_{tag}_{name}"] = f_product(sites, ph)
            metrics[f"f_ceiling_{tag}_{name}"] = f_ceiling(sites, ph)
    for name, sites in site_sets.items():
        if sites is None:
            continue
        metrics[f"ph_half_free_{name}"] = titration_midpoint([free for _c, free in sites])
        if not free_only:
            metrics[f"ph_half_complex_{name}"] = titration_midpoint(
                [complexed for complexed, _f in sites]
            )
    if free_only:
        return
    # The switch ratio F(pH) / F(reference), from the same f_product calls as the F
    # columns, so dividing the two columns by hand reproduces it exactly.
    for name, sites in site_sets.items():
        if sites is None:
            continue
        denominator = f_product(sites, reference_ph)
        for ph in ph_values:
            tag = ph_tag(ph)
            if tag != reference_tag:
                metrics[f"switch_{tag}_{name}"] = f_product(sites, ph) / denominator
    saturated: list[str] = []
    for factor in factors:
        tag = factor_tag(factor)
        for name, sites in site_sets.items():
            if sites is None:
                continue
            value, at_ceiling = release_ph(sites, factor)
            metrics[f"ph_release_{tag}_{name}"] = value
            if value is not None:
                # NA exactly when the release pH is NA: nothing for the claim to be about.
                metrics[f"f_dips_below_{tag}_{name}"] = int(
                    dips_below(sites, factor, value)
                )
            if at_ceiling:
                # NA because the answer lies ABOVE the scanned window, not because the
                # factor is never reached; named, not reported as a pH of 9.
                saturated.append(f"{tag}_{name}")
    metrics["ph_release_saturated"] = ";".join(saturated)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def write_result(
    out_dir: Path,
    name: str,
    status: str,
    path: str,
    metrics: dict,
    metric_columns: list[str] | None = None,
) -> None:
    """Write one design's TSV.

    ``metric_columns`` is the run's OWN column list (the pH-dependence names depend on
    --ph-values / --release-factors); it defaults to the list for the default flags,
    which is what the standalone and the fallback paths want.
    """
    columns = METRIC_COLUMNS if metric_columns is None else metric_columns
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / f"{name}.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "status", "path", *columns])
        w.writerow([name, status, path, *(_fmt(metrics.get(c)) for c in columns)])


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


class Args(argparse.Namespace):
    task_file: Path | None
    structure: Path | None
    name: str
    resnums: str
    binder_chains: str
    expect_binder_len: str
    expect_col: str
    centres_source: str
    free_only: bool
    input_column: str
    keep_propka_input: bool
    ph_values: str
    reference_ph: str
    release_factors: str
    out_dir: Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "PROPKA3 pKa of each designed histidine, in the complex and in the free "
            "binder (the binder chains alone, same coordinates)."
        )
    )
    ap.add_argument(
        "--task-file",
        type=Path,
        default=None,
        help="PIPELINE MODE. Tab-separated 'name<TAB>structure<TAB>centres<TAB>"
        "expect_col<TAB>expect_value' rows, one per design; centres is a comma-joined "
        "spec ('45,52', 'A:45,B:45' or 'auto'), empty slots allowed at the END only. "
        "The last two fields may be empty. Mutually exclusive with --structure.",
    )
    ap.add_argument(
        "--structure",
        type=Path,
        default=None,
        help="STANDALONE MODE. One PDB file to measure; use with --resnums. Writes "
        "the same <name>.tsv into --out-dir as a pipeline task would.",
    )
    ap.add_argument("--name", default="standalone", help="Standalone design name.")
    ap.add_argument(
        "--resnums",
        default="",
        help="STANDALONE MODE: the centre spec, e.g. '45,52', 'A:45,B:45' or 'auto' "
        "(auth_seq_id numbering).",
    )
    ap.add_argument(
        "--binder-chains",
        required=True,
        help="Comma-joined binder chain IDs (one or more) in the auth_asym_id "
        "namespace; together they are the free state.",
    )
    ap.add_argument(
        "--free-only",
        action="store_true",
        help="Run the free state only (no partner chain needed): complex-dependent "
        "columns are NA by construction and free_only is 1.",
    )
    ap.add_argument(
        "--expect-binder-len",
        default="",
        help="STANDALONE MODE: expected total residues of the binder chains (an integer "
        "or a one-letter sequence). Empty means the check is not made.",
    )
    ap.add_argument(
        "--expect-col",
        default="--expect-binder-len",
        help="Only used to name the column in the binder-length error message.",
    )
    ap.add_argument(
        "--centres-source",
        default="--resnums",
        help="Only echoed into the centres_source column and error messages.",
    )
    ap.add_argument(
        "--input-column",
        default="the input structure column",
        help="Only used to name the column in the 'could not resolve' error.",
    )
    ap.add_argument(
        "--keep-propka-input",
        action="store_true",
        help="Keep the temporary PDB files PROPKA read, as <out-dir>/<name>_propka/"
        "complex.pdb and free.pdb, so the numbers can be reproduced with `propka3`.",
    )
    ap.add_argument(
        "--ph-values",
        default=",".join(num_token(v, keep_decimal=True) for v in DEFAULT_PH_VALUES),
        help="Comma-separated pH values at which the weakening factor F is reported. "
        f"Default '{','.join(num_token(v, True) for v in DEFAULT_PH_VALUES)}'. Each "
        f"between 0 and 14, at most {MAX_PH_VALUES}, no duplicates.",
    )
    ap.add_argument(
        "--reference-ph",
        default=num_token(DEFAULT_REFERENCE_PH, keep_decimal=True),
        help="The pH the switch ratio divides by: the pH the binder is meant to HOLD "
        f"at. Default {num_token(DEFAULT_REFERENCE_PH, True)}. MUST be one of "
        "--ph-values.",
    )
    ap.add_argument(
        "--release-factors",
        default=",".join(
            num_token(v, keep_decimal=False) for v in DEFAULT_RELEASE_FACTORS
        ),
        help="Comma-separated weakening factors X for the release-pH columns. "
        f"Default "
        f"'{','.join(num_token(v, False) for v in DEFAULT_RELEASE_FACTORS)}'. Each "
        f"> 1, at most {MAX_RELEASE_FACTORS}, no duplicates.",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(namespace=Args())

    # Validated here as well as at submit time: the worker is runnable standalone.
    ph_values = parse_ph_values(args.ph_values)
    factors = parse_release_factors(args.release_factors)
    reference_ph = parse_reference_ph(args.reference_ph, ph_values)
    run_columns = metric_columns_for(ph_values, factors, reference_ph)

    binder_chains = [c.strip() for c in args.binder_chains.split(",") if c.strip()]
    if not binder_chains or any(len(c) != 1 for c in binder_chains):
        raise ValueError(
            f"--binder-chains must name one or more single-character chain IDs, got "
            f"{args.binder_chains!r}; PROPKA reads a single chain column"
        )
    if len(set(binder_chains)) != len(binder_chains):
        raise ValueError(f"--binder-chains lists a chain twice: {binder_chains}")
    if (args.task_file is None) == (args.structure is None):
        raise ValueError(
            "pass exactly one of --task-file (pipeline) or --structure (standalone)"
        )

    if args.task_file is not None:
        with open(args.task_file) as f:
            members = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    else:
        members = [
            [
                args.name,
                str(args.structure),
                args.resnums,
                args.expect_col,
                args.expect_binder_len,
            ]
        ]

    failures = 0
    for member in members:
        if len(member) < 3:
            raise ValueError(
                f"malformed task row {member!r}: expected "
                f"'name<TAB>structure<TAB>centres[<TAB>expect_col<TAB>expect_value]'"
            )
        member = member + [""] * (5 - len(member))
        name, src, centres_cell, expect_col, expect_raw = member[:5]
        expect_col = expect_col or "the expected-binder-length column"
        # The echo of the question survives every failure.
        metrics: dict = {
            "binder_chains": ",".join(binder_chains),
            "centres_source": args.centres_source,
            "free_only": int(args.free_only),
            "free_state_model": FREE_STATE_MODEL,
            # In the CANONICAL token form the column names are built from.
            "ph_values": ",".join(num_token(v, True) for v in ph_values),
            "reference_ph": num_token(reference_ph, True),
            "release_factors": ",".join(num_token(v, False) for v in factors),
            "f_model": F_MODEL,
        }
        try:
            if not src.strip():
                raise ValueError(
                    f"could not resolve a structure path from column "
                    f"{args.input_column!r} for this design -- neither on its own row "
                    f"nor on any ancestor (parent_name/parent_table)"
                )
            expect_len = _parse_expect(expect_raw, expect_col)
            keep_dir = (
                args.out_dir / f"{name}_propka" if args.keep_propka_input else None
            )
            measure(
                Path(src),
                binder_chains,
                centres_cell,
                args.centres_source,
                args.free_only,
                expect_len,
                expect_col,
                keep_dir,
                metrics,
                ph_values,
                factors,
                reference_ph,
                design=name,
            )
            write_result(args.out_dir, name, "OK", src, metrics, run_columns)
            print(_summary(name, metrics, ph_values, args.free_only))
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            failures += 1
            msg = " ".join(str(e).split())  # one line: it lands in a TSV cell
            print(f"{name}: ERROR {msg}", file=sys.stderr)
            # The echo and the trust facts survive; EVERY SCIENCE COLUMN IS EMPTY (NA).
            kept = {k: v for k, v in metrics.items() if k in TRUST_COLUMNS}
            write_result(args.out_dir, name, f"error: {msg}", "", kept, run_columns)

    if failures:
        print(
            f"propka: {failures}/{len(members)} design(s) in this task recorded an "
            f"error row",
            file=sys.stderr,
        )


def _summary(name: str, metrics: dict, ph_values: list[float], free_only: bool) -> str:
    """One log line per OK design."""
    head = (
        f"{name}: centres={metrics['centres_n']} his_ok={metrics['centres_his_ok']} "
        f"pka_ok={metrics['centres_pka_ok']} n_binder_his={metrics['n_binder_his']} "
        f"h_stripped={metrics['n_h_stripped']} "
        f"warnings={metrics.get('n_propka_warnings_complex', 'NA')}"
        f"/{metrics['n_propka_warnings_free']} no_pka={metrics['n_binder_his_no_pka']}"
    )
    if free_only:
        return f"{head} free_only | {metrics['binder_his_pka']}"
    f_summary = " ".join(
        f"{ph_tag(ph)}={metrics[f'f_{ph_tag(ph)}_designed']:.2f}" for ph in ph_values
    )
    return (
        f"{head} min_dpka={metrics['min_dpka']:.2f} "
        f"coupled={metrics['n_centres_coupled']} "
        f"terms={metrics['binder_term_groups_complex']}"
        f"|{metrics['binder_term_groups_free']} match={metrics['term_groups_match']} "
        f"| F {f_summary} | {metrics['binder_his_pka']} "
        f"| {metrics['per_centre_pka_detail']}"
    )


if __name__ == "__main__":
    main()
