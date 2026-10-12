"""Collect ``design_ph.py`` outputs into a child table -- one row per designed sequence.

Lineage: the submitter symlinks each input complex as ``grp_<g>/inputs/<design_name>.pdb``
and ``design_ph.py`` names every output file after the input stem, so

    <out_dir>/**/<parent_name>.metrics.csv
    <out_dir>/**/<parent_name>.fa
    <out_dir>/**/<parent_name>.ext.fa
    <out_dir>/**/<parent_name>.centers.json
    <out_dir>/**/<parent_name>.couplings.csv      (--redesign-scope coupling only)

and the stem IS the parent-table row. The search is a RECURSIVE glob on purpose: the
subdirectory design_ph.py writes into under ``--out_folder`` is not documented, so this
collector finds the files wherever they land rather than hardcoding a layout. A stem
whose files are in two places at once is an error, not a silent first-match.

THE JOIN IS THE RISK
--------------------
A row's sequence comes from ``<stem>.fa`` and its numbers from ``<stem>.metrics.csv``;
if those two orders ever diverge, every design silently inherits another design's
``ddE_mean``. So the join uses every scrap of evidence the two files share, and reports
what it did:

  1. **key, on the DESIGN ID** -- the primary ladder, and the one the real program
     supports. A design header's first comma-token is its id (``<stem>_lam0_s0``) and
     the metrics file repeats it in its ``name`` column. If those ids are unique on
     both sides and the two SETS are identical, rows are paired value by value.
     ``join_mode = "key"``, ``join_key = "name"``.
     If the ids are NOT unique -- which a multi-lambda sweep may cause, since whether
     ``0.25`` survives into the id is unverified -- the composite ``(id, lambda)`` is
     tried, both halves being present in both files; ``join_key = "name,lambda"``.
     Neither unique ⇒ **failed**, never positional.
     Then ``lambda`` and the header's ``mutations`` (against ``n_mutations``) are
     compared on every pair, and ``join_checked_fields`` lists ONLY the fields that
     actually had a value on both sides. The rounded measurement echoes
     (``dE_complex``, ``ddE``) are deliberately never compared: the header prints 3
     decimals and the CSV 4, and the real data contains the tie that proves it.
  2. **key, on shared ``key=value`` header fields** -- the legacy ladder, for a header
     format with no design id. Same rules.
  3. **positional** -- reached only when neither key exists at all: rows matched by
     order, after the counts are asserted equal and every shared field is compared on
     every pair. ``join_checked_fields`` empty here means nothing could corroborate it.
  4. **failed** -- counts differ, a key exists but the two value SETS disagree, a
     cross-checked field disagrees on some pair, a required file is missing, the stem is
     ambiguous, a metrics row has more fields than the header, a centre does not match
     the designed sequence (``centers_seq_check``), or ``centers.json`` and a row's own
     ``centers`` cell name different centres (``centers_match_metrics``). The collector
     then mints ONE row for that parent with an ``error`` status and the evidence. It
     never emits numbers it cannot place, and never writes 0 where it means "no
     measurement".

A unique key whose value sets disagree is NOT "no key available" -- it is positive
evidence that the two files are not in correspondence, so it fails rather than quietly
falling back to order.

THE METRICS FILE IS NOT ALL NUMBERS
-----------------------------------
READ FROM SOURCE at 7ffb756 (design_ph.py:1859-1877). ``dE_complex`` and ``dE_apo`` are
``;``-joined PER-CENTRE 4-dp strings -- a single number with one centre, ``"a;b"`` with
two, ``"nan;nan"`` on the apo side when there are no target chains. ``dE_complex_mean``
and ``dE_apo_mean`` are the numeric means, and
``ddE_mean = mean(dE_complex) - mean(dE_apo)``.

A column that is a float on some rows and a string on others cannot be filtered, so the
two per-centre columns are SPLIT: the numeric column keeps the single value (NA when
there are several centres) and ``<name>_centres`` carries the verbatim string.
``ddE_centres`` and ``ddE_min`` (the minimum over centres) are derived from them.

Two contradictions fail the parent: a per-centre count that disagrees with the row's own
``centers`` cell or with ``centers.json``, and a ``mean(per-centre ddE)`` more than 2e-3
from the file's ``ddE_mean``.

What it does NOT collect
------------------------
``summary.csv`` is ignored: it ranks designs across structures, which sorting the
prosapia table does natively and with the lineage intact. ``<stem>.couplings.csv`` stays
a FILE (``couplings_path``, ``""`` when absent): a per-position influence map is not one
value per design.
"""

import csv
import json
import re
from pathlib import Path
from typing import Any, Iterable, NamedTuple

import pandas as pd

from prosapia.core import Collected, CollectCtx, CollectEach, DesignCtx

# Sidecar written by the run (see run_atomiumph.py). Chains are recorded there because
# they are run-wide: without them the table cannot say WHICH molecule was redesigned.
META_FILENAME = ".meta.json"
EXPECTED_KEY = "atomiumph_expected_per_structure"
BINDER_CHAINS_KEY = "atomiumph_binder_chains"
TARGET_CHAINS_KEY = "atomiumph_target_chains"
CENTER_MODE_KEY = "atomiumph_center_mode"
OUTPUT_SCOPE_KEY = "atomiumph_output_scope"
# VERIFIED 2026-10-09: the model-noise level decides which centres design_ph picks
# (two runs over the same structures chose different centres at n03 and n05), and the
# dir-label was until now its ONLY record. Written by the run from here on; for a leaf
# written before that, recovered from the sub-manifest -- see ``noise_and_lam_from_tasks``.
MODEL_NOISE_KEY = "atomiumph_model_noise"
LAM_SWEEP_KEY = "atomiumph_lam_sweep"
# Same reason as the noise level: it changes which positions were allowed to move, so a
# table that mixes two values of it is not ranking like against like.
PROTECT_TITRATABLE_KEY = "atomiumph_protect_titratable"
# Modelled only because the upstream client defaults it to 0.5 where design_ph.py
# defaults it to 0.0. Above 0 the apo graph enters sampling, so two values of it are
# two different experiments.
APO_WEIGHT_KEY = "atomiumph_apo_weight"

# Header fields that IDENTIFY a draw. Only these may key or cross-check a join; a
# shared field that is a measurement (``seq_rec=0.43`` against a full-precision CSV
# value) would contradict a correct pairing. See ``is_identity_field``.
IDENTITY_FIELDS = frozenset(
    {"lambda", "lam", "sample", "seed", "index", "design", "id", "name"}
)

# Written per subgroup by atomiumph.sh from the image's baked /opt/atomium_ph/COMMIT.
COMMIT_FILENAME = "image_commit.txt"

# Metrics columns that are ``;``-joined PER-CENTRE strings, not numbers.
# READ FROM SOURCE at 7ffb756 (design_ph.py:1859-1860):
#     "dE_complex": ";".join(f"{v:.4f}" for v in de_complex),
#     "dE_apo":     ";".join(f"{v:.4f}" for v in de_apo),
# With one centre the cell is a single 4-dp number and reads as one; with several it is
# "0.1234;0.5678" and would land in the same column as a STRING. A column that is
# sometimes a float and sometimes a string cannot be filtered, so each of these is split:
# the numeric column keeps the single value (NA when there are several) and a sibling
# ``<name>_centres`` column carries the verbatim string.
PER_CENTRE_COLUMNS = ("dE_complex", "dE_apo")

# Marker csv.DictReader files surplus fields under. A row with more fields than the
# header has shifted values; it is an error, not something to trim.
OVERFLOW_KEY = "__overflow__"

# Suffixes design_ph.py appends to the input stem. Order matters: ".ext.fa" must be
# tested before ".fa" or every ext FASTA would be filed under the stem "<name>.ext".
SUFFIXES = (
    ("metrics", ".metrics.csv"),
    ("couplings", ".couplings.csv"),
    ("centers", ".centers.json"),
    ("ext_fa", ".ext.fa"),
    ("fa", ".fa"),
)

# Directories the run creates that never hold design_ph output.
SKIP_DIRS = ("atomiumph_logs", "atomiumph_tasks", "inputs")

# Headers a FASTA preamble entry may use for the echoed input sequence. Matched
# EXACTLY (case-insensitively, after stripping '>' and whitespace) rather than by a
# substring: a loose match would silently drop a real design.
#
# VERIFIED 2026-10-09 against a real run: design_ph.py writes the echo as
# ``><stem>_native`` -- the plain name QUALIFIED BY THE INPUT STEM. ``drop_preamble``
# therefore also accepts ``<stem>_<p>`` for every ``p`` below, with the stem passed in
# by the caller. Still an exact match on the whole first comma-token, never a
# substring: a real design header that merely CONTAINS "native" survives, and so does
# one belonging to another stem.
PREAMBLE_HEADERS = {
    "_init_sequence",
    "init",
    "init_sequence",
    "native",
    "native_sequence",
    "input",
    "input_sequence",
    "seed",
    "seed_sequence",
    "wt",
}

# Column names this collector owns. A metrics column colliding with one of these is
# emitted as ``metrics_<name>`` instead of overwriting it -- a trust column must never
# be shadowed by the data it vouches for.
RESERVED = {
    "join_mode",
    "join_key",
    "join_checked_fields",
    "ext_join_mode",
    "ext_checked_fields",
    "ext_status",
    "center_mode",
    "output_scope",
    "n_metrics_rows",
    "n_fa_entries",
    "n_ext_entries",
    "n_expected",
    "n_preamble_dropped",
    "entry_index",
    "binder_chains",
    "target_chains",
    "image_commit",
    "sequence",
    "sequence_ext",
    "seq_len",
    "centers",
    "n_centers",
    "centers_source",
    "centers_resnum_key",
    "placement_score_min",
    "charge_field_min",
    "charge_field_max",
    "n_redesign",
    "centers_seq_check",
    "centers_match_metrics",
    "dE_complex_centres",
    "dE_apo_centres",
    "ddE_centres",
    "ddE_min",
    "placement_metric",
    "redesign_scope",
    "output_predates_meta",
    "native_sequence",
    "n_diff_native",
    "model_noise",
    "lam_sweep",
    "protect_titratable",
    # design_ph.py echoes apo_weight into the metrics CSV too (design_ph.py:1869).
    # The run-recorded value keeps the plain name -- it is what the SUBMIT asked for --
    # and the program's echo becomes metrics_apo_weight, so the two can be compared
    # instead of one silently overwriting the other.
    "apo_weight",
    "metrics_path",
    "couplings_path",
    "header",
    # VERIFIED 2026-10-09: <stem>.metrics.csv really does carry a ``name`` column (the
    # design id) and a ``centers`` column ("10:HIS-S"). Both are reserved so they are
    # renamed to ``metrics_name`` / ``metrics_centers`` instead of landing on, or being
    # overwritten by, this collector's own columns. See EXPECTED_COLLISIONS.
    "name",
}

# Metrics columns KNOWN to collide, so their rename is normal and is not announced. An
# unexpected collision still prints, because it means upstream grew a column that
# shadows a trust column and a reader needs to know the table layout moved.
EXPECTED_COLLISIONS = {"name", "centers", "apo_weight"}


# --------------------------------------------------------------------------- #
# file discovery
# --------------------------------------------------------------------------- #


def index_outputs(out_dir: Path) -> dict[str, dict[str, list[Path]]]:
    """Map ``stem -> kind -> [paths]`` for every design_ph output under ``out_dir``.

    Recursive, because the subdirectory layout under ``--out_folder`` is undocumented.
    Lists, not single paths, so an ambiguous stem (the same name written twice) can be
    reported as an error instead of resolved by luck.
    """
    index: dict[str, dict[str, list[Path]]] = {}
    for path in sorted(out_dir.rglob("*")):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.relative_to(out_dir).parts[:-1]):
            continue
        name = path.name
        for kind, suffix in SUFFIXES:
            if name.endswith(suffix) and len(name) > len(suffix):
                stem = name[: -len(suffix)]
                index.setdefault(stem, {}).setdefault(kind, []).append(path)
                break
    return index


def commit_for(path: Path, out_dir: Path) -> str:
    """The image commit recorded beside this output, walking up to ``out_dir``.

    atomiumph.sh writes ``image_commit.txt`` into each subgroup dir, so a table can
    distinguish two runs built from different pins of a private, moving branch.
    """
    here = path.parent
    while True:
        marker = here / COMMIT_FILENAME
        if marker.is_file():
            try:
                return marker.read_text().strip()
            except OSError:
                return ""
        if here == out_dir or here.parent == here:
            return ""
        here = here.parent


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #


def parse_fasta(path: Path, join: str = "") -> list[tuple[str, str]]:
    """``[(header, sequence)]`` in file order; the header excludes the leading '>'.

    ``join`` is what a WRAPPED sequence's continuation lines are joined with. ``""`` is
    right for a canonical FASTA, where the sequence is one letter per character, and
    WRONG for ``.ext.fa``, whose sequence is space-separated tokens: joining
    ``"… GLU-D"`` to ``"HIS-S …"`` with nothing fuses two residues into the single
    nonsense token ``GLU-DHIS-S``. The real files are unwrapped, so this never bit --
    which is exactly why it had to be fixed before a longer binder makes them wrapped.
    Call it with ``join=" "`` for the extended FASTA.
    """
    entries: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                entries.append((header, join.join(seq)))
            header, seq = line[1:].strip(), []
        elif line and header is not None:
            seq.append(line)
    if header is not None:
        entries.append((header, join.join(seq)))
    return entries


def parse_header(header: str) -> dict[str, str]:
    """``key=value`` fields of a FASTA header, comma-separated.

    ``gen_sample_0, lambda=0.25, ddE=1.2`` -> ``{'gen_sample': '0', 'lambda': '0.25',
    'ddE': '1.2'}``. A bare trailing ``_<int>`` token (AtomiUM's own style) becomes a
    field so it can serve as a join key; any other bare token is ignored.
    """
    fields: dict[str, str] = {}
    for part in header.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            fields[key.strip()] = value.strip()
            continue
        head, _, tail = part.rpartition("_")
        if head and tail.isdigit():
            fields[head] = tail
    return fields


def read_metrics(path: Path) -> tuple[list[str], list[dict[str, str]], int]:
    """``(column order, rows, n_overflow)`` of a metrics CSV, every value kept as text.

    Text on purpose: the coercion to numbers happens once, visibly, in ``_as_number`` --
    pandas' own inference would turn an empty cell into NaN and a '0' into 0.0 without
    distinguishing them.

    ``n_overflow`` counts rows carrying MORE fields than the header. Those rows' values
    are shifted, so the caller must treat them as an error; a short row (fewer fields)
    is fine and its missing values become NA.
    """
    n_overflow = 0
    with open(path, newline="") as f:
        reader = csv.DictReader(f, restkey=OVERFLOW_KEY)
        columns = [c for c in (reader.fieldnames or []) if c]
        rows: list[dict[str, str]] = []
        for raw in reader:
            if raw.get(OVERFLOW_KEY):
                n_overflow += 1
            rows.append(
                {
                    k: (v if v is not None else "")
                    for k, v in raw.items()
                    if k and k != OVERFLOW_KEY
                }
            )
    return columns, rows, n_overflow


def _as_number(text: str) -> Any:
    """'3' -> 3, '0.25' -> 0.25, 'abc' -> 'abc', '' -> NA. Never 0 for a missing value."""
    text = text.strip()
    if text == "":
        return pd.NA
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def _norm(text: str) -> str:
    """Join-key normalisation: numbers compare as numbers ('0.250' == '0.25')."""
    text = text.strip()
    try:
        return repr(float(text))
    except ValueError:
        return text.casefold()


# --------------------------------------------------------------------------- #
# the join
# --------------------------------------------------------------------------- #


class Join(NamedTuple):
    """The outcome of pairing FASTA entries with metrics rows.

    ``mode``   "key" | "positional" | "failed".
    ``keys``   the fields rows were MATCHED on (key mode only).
    ``checked`` the shared fields that were COMPARED value-by-value. In key mode these
                are the keys; in positional mode a non-empty list means the order was
                corroborated, an empty one means nothing could corroborate it.
    ``pairs``  (entry, metrics row) in FASTA order; empty when mode is "failed".
    ``note``   why it failed, for the error row's status.
    """

    mode: str
    keys: list[str]
    checked: list[str]
    pairs: list[tuple[tuple[str, str], dict[str, str]]]
    note: str = ""


def is_identity_field(name: str, bare: set[str]) -> bool:
    """Is this header field an IDENTIFIER of a draw, rather than a measurement?

    Only identifiers may key or cross-check a join. A header field that is a *metric*
    is worse than useless here: the sibling tool's FASTA embeds values like
    ``seq_rec=0.43``, and comparing that rounded text against the CSV's full-precision
    ``0.431872`` would contradict a perfectly correct pairing and fail every parent.
    """
    return name.casefold() in IDENTITY_FIELDS or name in bare


def bare_field_names(header: str) -> set[str]:
    """Fields ``parse_header`` synthesised from a bare ``<word>_<int>`` token.

    Those are positional identifiers by construction (``gen_sample_0``), so they are
    identity-like however they are named.
    """
    names: set[str] = set()
    for part in header.split(","):
        part = part.strip()
        if not part or "=" in part:
            continue
        head, _, tail = part.rpartition("_")
        if head and tail.isdigit():
            names.add(head)
    return names


def shared_fields(
    headers: list[dict[str, str]], columns: list[str], bare: set[str]
) -> list[str]:
    """Identity fields that EVERY header carries and a metrics column also names."""
    if not headers:
        return []
    lower_to_col = {c.casefold(): c for c in columns}
    return [
        k
        for k in headers[0]
        if k.casefold() in lower_to_col
        and all(k in h for h in headers)
        and is_identity_field(k, bare)
    ]


def first_token(header: str) -> str:
    """The first comma-separated token of a FASTA header.

    VERIFIED 2026-10-09: for a design entry this is the design id and it is exactly the
    ``name`` column of ``<stem>.metrics.csv`` (``tnf_15_a2_model_0_lam0_s0``); for the
    preamble entry it is ``<stem>_native``.
    """
    return header.split(",")[0].strip()


# The stage-1 seeds design_ph.py writes BEFORE the designs under --seed-mode generate.
# READ FROM SOURCE at 7ffb756 (design_ph.py:1772-1779)::
#
#     if args.seed_mode == "generate":
#         for copy_index in range(seeds.shape[0]):
#             fasta_lines.append(f">{name}_stage1_seed_s{copy_index}\n" + ...)
#
# One per --num-seq, after the ``_native`` echo. They are NOT designs: they are what
# stage 1 generated before the selective infill ran. Unrecognised, they would make
# ``n_fa_entries`` exceed ``n_metrics_rows`` by num_seq and fail EVERY parent -- so a
# whole seed mode would be unusable. Matched by an anchored pattern on the full first
# comma-token, with the known stem, so a design can never match it.
_STAGE1_SEED_RE = r"stage1_seed_s\d+"


def is_preamble_header(header: str, stem: str | None) -> bool:
    """Is this header a non-design entry (the input echo, or a stage-1 seed)?

    EXACT match on the whole first comma-token, against the bare names and -- this is
    the part the real program needed -- against ``<stem>_<name>``, which is what
    design_ph.py actually writes (``>tnf_15_a2_model_0_native``), plus
    ``<stem>_stage1_seed_s<int>`` under ``--seed-mode generate``. Never a substring
    test: ``>native_like_lam0_s0`` is a design and must survive, and so must a design
    whose id merely contains "seed".
    """
    token = first_token(header).casefold()
    if token in PREAMBLE_HEADERS:
        return True
    if not stem:
        return False
    prefix = f"{stem.casefold()}_"
    if not token.startswith(prefix):
        return False
    rest = token[len(prefix) :]
    return rest in PREAMBLE_HEADERS or re.fullmatch(_STAGE1_SEED_RE, rest) is not None


def drop_preamble(
    entries: list[tuple[str, str]], stem: str | None = None
) -> tuple[list[tuple[str, str]], int]:
    """Strip a leading echo of the input sequence, if the header says that is what it is."""
    dropped = 0
    while entries and is_preamble_header(entries[0][0], stem):
        entries = entries[1:]
        dropped += 1
    return entries, dropped


# Metrics columns that may hold the design id, i.e. the FASTA header's first token.
# VERIFIED 2026-10-09: design_ph.py writes ``name``. Kept as a tuple so a rename
# upstream is a one-line change, but deliberately short -- an unverified guess here
# would key the join on something that is not an identifier.
DESIGN_ID_COLUMNS = ("name",)

# Header field -> (metrics column, comparison) for the fields cross-checked AFTER rows
# have been paired on the design id. An allowlist, not "everything shared", because:
#
#   the design header prints the measurements ROUNDED TO 3 DECIMALS
#   (``dE_complex=0.720``) while the CSV keeps 4 (``0.7195``) and a full-precision mean
#   besides. Comparing those would FAIL A CORRECT PAIRING -- and the real data already
#   contains the exact tie at the rounding boundary (0.720 vs 0.7195) that proves it.
#
# Only identifiers are compared: ``lambda`` numerically (header ``0`` == csv ``0.0``)
# and the header's ``mutations`` against the csv's ``n_mutations`` as integers.
CROSS_CHECKS: tuple[tuple[str, str, str], ...] = (
    ("lambda", "lambda", "number"),
    ("mutations", "n_mutations", "int"),
)


def _agree(header_text: str, metrics_text: str, kind: str) -> bool:
    """Do a header field and a metrics cell carry the same identifier value?"""
    a, b = header_text.strip(), metrics_text.strip()
    if a == "" or b == "":
        return True  # nothing to contradict; absence is not disagreement
    try:
        if kind == "int":
            return int(round(float(a))) == int(round(float(b)))
        return float(a) == float(b)
    except ValueError:
        return _norm(a) == _norm(b)


def _design_id_join(
    entries: list[tuple[str, str]],
    headers: list[dict[str, str]],
    columns: list[str],
    rows: list[dict[str, str]],
) -> Join | None:
    """Pair on the design id -- the header's first token and the metrics ``name`` cell.

    Returns ``None`` when this evidence is not available (no id column, or the header
    ids are not unique), so the caller falls through to the shared-field logic. Returns
    a FAILED join when the id IS available but the two files name different draws or a
    cross-checked identifier disagrees: that is a contradiction, and a contradiction is
    never answered with a positional guess.

    This is the join the real program supports. It is also the only one that can tell
    two draws apart at ``lambda 0`` with ``--temperature 0.1``, where both draws can
    come back with the SAME sequence and the SAME ddE to 6 decimals and differ in
    nothing but ``_s0`` / ``_s1``.
    """
    if not entries or not rows:
        return None
    lower_to_col = {c.casefold(): c for c in columns}
    id_col = next(
        (lower_to_col[c] for c in DESIGN_ID_COLUMNS if c in lower_to_col), None
    )
    if id_col is None:
        return None

    lam_col = lower_to_col.get("lambda")
    h_ids: list[tuple[str, ...]] = [(first_token(h),) for h, _ in entries]
    m_ids: list[tuple[str, ...]] = [(str(r.get(id_col, "")).strip(),) for r in rows]
    key_names = [id_col]

    # THE COMPOSITE FALLBACK. Whether design_ph renders lambda 0.25 into the design id
    # as `lam0.25` or truncates it to `lam0` is UNVERIFIED -- the only real output seen
    # had a single lambda, so a multi-lambda sweep could produce the SAME first token
    # for several draws. If it does, the id alone is not a key, but (id, lambda) still
    # is: both are present in both files. Falling back to POSITION here instead would
    # pair draws that differ in nothing the collector can see, which is exactly the
    # silent mis-mapping this tool exists to prevent.
    if not all(i[0] for i in h_ids) or len(set(h_ids)) != len(h_ids):
        if lam_col is None or not all("lambda" in h for h in headers):
            return None  # genuinely no key; the shared-field logic may still find one
        h_ids = [(first_token(h), _norm(hdr["lambda"]))
                 for (h, _), hdr in zip(entries, headers)]
        m_ids = [(str(r.get(id_col, "")).strip(), _norm(str(r.get(lam_col, ""))))
                 for r in rows]
        key_names = [id_col, lam_col]
        if not all(i[0] for i in h_ids) or len(set(h_ids)) != len(h_ids):
            return Join(
                "failed", [], key_names, [],
                f"the FASTA design ids are not unique, and nor is "
                f"({', '.join(key_names)}): {len(set(h_ids))} distinct key(s) over "
                f"{len(h_ids)} entries, so no two draws can be told apart",
            )

    label = ",".join(key_names)
    if not all(i[0] for i in m_ids) or len(set(m_ids)) != len(m_ids):
        return Join(
            "failed", [], key_names, [],
            f"the metrics key ({label}) is not unique ({len(set(m_ids))} distinct "
            f"value(s) over {len(m_ids)} row(s)), so it cannot identify a draw",
        )
    if set(h_ids) != set(m_ids):
        missing = sorted(":".join(k) for k in set(h_ids) - set(m_ids))
        shown = ", ".join(missing[:3]) + (" …" if len(missing) > 3 else "")
        return Join(
            "failed", [], key_names, [],
            f"FASTA headers and the metrics ({label}) name different draws "
            f"({len(set(h_ids) - set(m_ids))} FASTA id(s) absent from the metrics "
            f"file: {shown})",
        )

    by_id = {m: r for m, r in zip(m_ids, rows)}
    pairs = [(entry, by_id[i]) for entry, i in zip(entries, h_ids)]

    checked = list(key_names)
    for field, column, kind in CROSS_CHECKS:
        col = lower_to_col.get(column.casefold())
        if col is None or not any(field in h for h in headers):
            continue
        compared = 0
        for i, (header, row) in enumerate(zip(headers, [r for _, r in pairs])):
            cell = str(row.get(col, "")).strip()
            # F6: only count a field as CHECKED where both sides actually had a value.
            # A field present in the header but blank in every metrics cell compared
            # nothing, and must not appear in join_checked_fields as if it had.
            if field not in header or header[field].strip() == "" or cell == "":
                continue
            compared += 1
            if not _agree(header[field], cell, kind):
                return Join(
                    "failed", [], checked + [col], [],
                    f"design id pairing contradicted at entry {i} "
                    f"({':'.join(h_ids[i])!r}): header says {field}={header[field]!r}, "
                    f"metrics column {col!r} says {cell!r}",
                )
        if compared and col not in checked:
            checked.append(col)
    return Join("key", key_names, checked, pairs)


def align_entries(
    entries: list[tuple[str, str]], columns: list[str], rows: list[dict[str, str]]
) -> Join:
    """Pair FASTA entries with metrics rows, using every shared IDENTITY field as
    evidence (see ``is_identity_field`` for why a shared metric field is excluded)."""
    raw = [h for h, _ in entries]
    headers = [parse_header(h) for h in raw]
    bare = set().union(*(bare_field_names(h) for h in raw)) if raw else set()

    # Strongest evidence first: the design id, which the real program writes on both
    # sides. Falls through (None) only when it is genuinely unavailable.
    id_join = _design_id_join(entries, headers, columns, rows)
    if id_join is not None:
        return id_join

    shared = shared_fields(headers, columns, bare)
    lower_to_col = {c.casefold(): c for c in columns}

    def metric(row: dict[str, str], field: str) -> str:
        return row.get(lower_to_col[field.casefold()], "")

    if shared and rows:
        h_keys = [tuple(_norm(h[k]) for k in shared) for h in headers]
        m_keys = [tuple(_norm(metric(r, k)) for k in shared) for r in rows]
        unique = len(set(h_keys)) == len(h_keys) and len(set(m_keys)) == len(m_keys)
        if unique:
            if set(h_keys) != set(m_keys):
                # A usable key whose two value sets disagree is POSITIVE EVIDENCE that
                # the files are not in correspondence. Falling back to order here would
                # be answering a contradiction with a guess.
                return Join(
                    "failed", [], shared, [],
                    f"shared field(s) {', '.join(shared)} are unique in both files but "
                    f"name different draws ({len(set(h_keys) - set(m_keys))} FASTA "
                    f"value(s) absent from the metrics file)",
                )
            by_key = {tuple(_norm(metric(r, k)) for k in shared): r for r in rows}
            pairs = [
                (entry, by_key[tuple(_norm(h[k]) for k in shared)])
                for entry, h in zip(entries, headers)
            ]
            return Join("key", shared, shared, pairs)

    if len(entries) != len(rows):
        # `checked` is [] here on purpose: nothing was compared against anything, so
        # naming fields would claim corroboration that never happened.
        return Join(
            "failed", [], [], [],
            f"{len(entries)} FASTA entries vs {len(rows)} metrics rows, and no shared "
            f"join key in the headers",
        )

    pairs = list(zip(entries, rows))
    # A shared field that is not jointly unique (``lambda`` is exactly that, and is the
    # field most likely to be present in both) still CORROBORATES a positional pairing:
    # compare it on every pair rather than discarding it.
    for field in shared:
        for i, ((header, _), row) in enumerate(pairs):
            want, got = _norm(headers[i][field]), _norm(metric(row, field))
            if want != got:
                return Join(
                    "failed", [], shared, [],
                    f"positional pairing contradicted at entry {i}: header says "
                    f"{field}={headers[i][field]!r}, metrics row says "
                    f"{metric(row, field)!r} ({header!r})",
                )
    return Join("positional", [], shared, pairs)


class ExtJoin(NamedTuple):
    """The outcome of pairing the extended FASTA onto the canonical one.

    Same four states and the same evidence rules as ``Join``: ``mode`` is "key",
    "positional", "missing" or "failed", ``checked`` names the identity fields compared
    value-by-value (empty with "positional" means nothing could corroborate the order),
    ``seqs`` are the extended sequences in the order of ``entries``.
    """

    mode: str
    checked: list[str]
    seqs: list[str]
    note: str = ""


def align_ext(
    entries: list[tuple[str, str]], ext_entries: list[tuple[str, str]]
) -> ExtJoin:
    """Extended sequences in the order of ``entries``.

    A full mirror of ``align_entries``, because ``sequence_ext`` is as mappable-wrong as
    ``ddE_mean`` is: match on shared identity fields when they are jointly unique; FAIL
    when such a key exists but the two files name different draws; otherwise pair by
    position with the counts asserted equal AND every shared identity field compared on
    every pair. The caller records the outcome in ``ext_status``/``ext_checked_fields``
    rather than letting a reader assume ``sequence_ext`` was verified.
    """
    if not ext_entries:
        return ExtJoin("missing", [], [])
    raw = [h for h, _ in entries]
    ext_raw = [h for h, _ in ext_entries]
    headers = [parse_header(h) for h in raw]
    ext_headers = [parse_header(h) for h in ext_raw]
    bare = set().union(*(bare_field_names(h) for h in raw + ext_raw))

    # Strongest evidence first, exactly as in _design_id_join: the header's first
    # token. VERIFIED 2026-10-09 -- .ext.fa repeats the design headers verbatim (minus
    # the native echo), so this is the join the real program supports.
    h_ids: list[tuple[str, ...]] = [(first_token(h),) for h in raw]
    e_ids: list[tuple[str, ...]] = [(first_token(h),) for h in ext_raw]
    key_names = ["design_id"]
    usable = (
        all(i[0] for i in h_ids)
        and all(i[0] for i in e_ids)
        and len(set(h_ids)) == len(h_ids)
        and len(set(e_ids)) == len(e_ids)
    )
    # Same composite fallback as _design_id_join: if the id alone repeats (which a
    # multi-lambda sweep may do -- see there), (id, lambda) is still a key, and pairing
    # by POSITION instead would be a guess dressed as a mapping.
    if not usable and all("lambda" in h for h in headers + ext_headers):
        h_ids = [(first_token(r), _norm(h["lambda"])) for r, h in zip(raw, headers)]
        e_ids = [(first_token(r), _norm(h["lambda"]))
                 for r, h in zip(ext_raw, ext_headers)]
        key_names = ["design_id", "lambda"]
        usable = (
            all(i[0] for i in h_ids)
            and all(i[0] for i in e_ids)
            and len(set(h_ids)) == len(h_ids)
            and len(set(e_ids)) == len(e_ids)
        )
        if not usable:
            return ExtJoin(
                "failed", key_names, [],
                f"neither the design id nor ({', '.join(key_names)}) is unique across "
                f"the two FASTAs, so no two draws can be told apart",
            )
    if usable:
        if set(h_ids) != set(e_ids):
            return ExtJoin(
                "failed", key_names, [],
                f"the canonical and extended FASTAs name different draws "
                f"({len(set(h_ids) - set(e_ids))} canonical id(s) absent from "
                f"the .ext.fa)",
            )
        by_id = {i: h for i, h in zip(e_ids, ext_headers)}
        seq_by_id = {i: seq for i, (_, seq) in zip(e_ids, ext_entries)}
        checked = list(key_names)
        for field, _column, kind in CROSS_CHECKS:
            if field in checked or not any(field in h for h in headers):
                continue
            compared = 0
            for i, (h, ident) in enumerate(zip(headers, h_ids)):
                other = by_id[ident]
                # F6: both sides must actually carry a value for this to count.
                if (field not in h or field not in other
                        or h[field].strip() == "" or other[field].strip() == ""):
                    continue
                compared += 1
                if not _agree(h[field], other[field], kind):
                    return ExtJoin(
                        "failed", checked + [field], [],
                        f"design id pairing contradicted at entry {i} "
                        f"({':'.join(ident)!r}): canonical header says "
                        f"{field}={h[field]!r}, extended header says "
                        f"{other[field]!r}",
                    )
            if compared:
                checked.append(field)
        return ExtJoin("key", checked, [seq_by_id[i] for i in h_ids])

    shared = (
        [
            k
            for k in headers[0]
            if all(k in h for h in headers + ext_headers) and is_identity_field(k, bare)
        ]
        if headers
        else []
    )

    if shared:
        h_keys = [tuple(_norm(h[k]) for k in shared) for h in headers]
        e_keys = [tuple(_norm(h[k]) for k in shared) for h in ext_headers]
        if len(set(h_keys)) == len(h_keys) and len(set(e_keys)) == len(e_keys):
            if set(h_keys) != set(e_keys):
                return ExtJoin(
                    "failed", shared, [],
                    f"shared field(s) {', '.join(shared)} are unique in both FASTAs but "
                    f"name different draws",
                )
            by_key = {k: seq for k, (_, seq) in zip(e_keys, ext_entries)}
            return ExtJoin("key", shared, [by_key[k] for k in h_keys])

    if len(entries) != len(ext_entries):
        return ExtJoin(
            "failed", [], [],
            f"{len(entries)} canonical sequences vs {len(ext_entries)} extended ones",
        )

    for field in shared:
        for i, (h, e) in enumerate(zip(headers, ext_headers)):
            if _norm(h[field]) != _norm(e[field]):
                return ExtJoin(
                    "failed", shared, [],
                    f"positional pairing contradicted at entry {i}: canonical header "
                    f"says {field}={h[field]!r}, extended header says {e[field]!r}",
                )
    return ExtJoin("positional", shared, [seq for _, seq in ext_entries])


# --------------------------------------------------------------------------- #
# centers.json
# --------------------------------------------------------------------------- #

_CENTER_LIST_KEYS = ("centers", "selected_centers", "chosen_centers", "pinned_centers")
_CHAIN_KEYS = ("chain", "chain_id", "chain_name", "binder_chain", "asym_id")
# ORDERED BY NUMBERING CONVENTION, most explicit first. Which one matched is reported
# in ``centers_resnum_key``: this is the only place in the tool that infers a numbering
# convention. Read that column before mapping a centre onto a structure.
#
# VERIFIED 2026-10-09 against a real run and the input CIF: design_ph.py writes
# ``position``, and it is a ZERO-BASED INDEX INTO THE BINDER CHAIN, not a PDB residue
# number. (tnf: position 10 is auth residue 11, MET in the input, and the designed
# sequence carries H at index 10. f3: position 0 is auth residue 1, SER in the input,
# and the designed sequence carries H at index 0.) An index-like key is therefore
# reported with an explicit ``_index0`` suffix, and the centre is checked against the
# designed sequence -- see ``check_centers_against_sequence``.
_RESNUM_KEYS = ("resnum", "res_num", "res_id", "resi", "residue", "seqpos", "position", "pos")
# Keys whose value this tool treats as a 0-based index into the designed sequence. Only
# these get the ``_index0`` label and only these are sequence-checked; an explicit
# ``resnum``-style key is left unlabelled and unchecked, because nothing here has
# verified ITS base.
_INDEX0_KEYS = frozenset({"position", "pos", "seqpos", "index", "idx"})
# Keys that may carry the centre's protonation TOKEN. ``token`` is the verified one.
# The others are accepted for the label, but only a value that LOOKS like a full token
# (``HIS-S``, or a bare 3-letter code) is ever read as a residue -- see
# ``_center_token``. A ``state`` of "S" is a protonation state, not serine.
_STATE_KEYS = ("token", "state", "target_state", "microstate", "residue_state", "aa")
# How much the model WANTED this centre where it put it. VERIFIED 2026-10-09: the real
# file records ``placement_score`` per centre (tnf 0.708, f3 0.093 -- the f3 His site
# was barely wanted at all). Collected as the MINIMUM over the centres, because the
# weakest centre is what decides whether the placement was any good.
#
# ONE key, not a list of plausible aliases: a bare ``score`` in some other schema could
# be any quantity at all, and a mis-sourced placement score reads exactly like a real
# one. If upstream renames it this is a one-line change, and until then the column is
# NA -- which is the honest answer.
_PLACEMENT_KEYS = ("placement_score",)

# Keys that record which positions were allowed to move. ``free_positions`` is the one
# the real program writes (VERIFIED 2026-10-09: len 12 for tnf and 6 for f3, matching
# "Stage 2: redesign positions = 12 / 6" in the task .out).
_REDESIGN_KEYS = (
    "n_redesign",
    "free_positions",
    "redesign_set",
    "redesign_positions",
    "redesign",
)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
ONE_LETTERS = frozenset(THREE_TO_ONE.values())


def residue_letter(token: Any) -> str | None:
    """The 1-letter residue a protonation TOKEN names, or ``None`` if unrecognised.

    DERIVED FROM THE TOKEN, never from ``--center-mode``: ``HIS-S``/``HIS-P`` -> ``H``,
    ``ASP-D``/``ASP-P`` -> ``D``, ``GLU-D`` -> ``E``. A bare 3-letter code works too.
    ``None`` means "cannot check this", which is NA -- not a failure.

    **A bare 1-letter value is deliberately NOT accepted.** In several schemas a
    ``state`` field holds a protonation state letter, and reading ``"state": "S"`` as
    serine would invent a residue and then fail a perfectly correct centre on it. Only
    a full token (``XXX`` or ``XXX-Y``) is read; see ``_center_token``, which is what
    decides which JSON field is allowed to be a token at all.
    """
    head = str(token).split("-")[0].strip().upper()
    if len(head) != 3:
        return None
    return THREE_TO_ONE.get(head)


def _center_token(value: Any) -> str | None:
    """``value`` if it looks like a protonation token, else ``None``.

    A token is ``XXX`` or ``XXX-Y`` with a 3-letter head. Anything shorter is a state
    flag or a letter code whose meaning this tool has not verified, and reading it as a
    residue would invent a check rather than perform one.
    """
    if value is None:
        return None
    text = str(value).strip()
    head = text.split("-")[0]
    return text if len(head) == 3 and head.isalpha() else None


def _parse_metrics_centers(cell: Any) -> set[tuple[int, str]] | None:
    """``{(index, token)}`` parsed from a metrics ``centers`` cell, or ``None``.

    The verified single-centre form is ``10:HIS-S``. **How several centres are joined
    is UNVERIFIED** -- the only real output had one -- so every plausible separator is
    accepted (``;`` ``,`` ``|`` and whitespace). ``None`` means "could not be read as
    (index, token) pairs", which the caller turns into NA, never into a failure: a
    separator nobody has seen must not condemn a correct run.
    """
    if cell is None:
        return None
    text = str(cell).strip()
    if not text:
        return None
    out: set[tuple[int, str]] = set()
    for part in re.split(r"[;,|\s]+", text):
        part = part.strip()
        if not part:
            continue
        index_text, sep, token = part.partition(":")
        if not sep:
            return None
        index = _as_index(index_text)
        tok = _center_token(token)
        if index is None or tok is None:
            return None
        out.add((index, tok.upper()))
    return out or None


class Centers(NamedTuple):
    """``columns`` to put on every child row, ``checks`` to verify against each design.

    ``checks`` is ``[(index0, letter, token)]`` for the centres whose numbering key is
    index-like AND whose token names a residue. ``pairs`` is ``{(index0, token)}`` for
    every centre that had both, used to cross-check each metrics row's own ``centers``
    cell. Both are kept out of ``columns`` on purpose: they are evidence for a check,
    not columns.
    """

    columns: dict[str, Any]
    checks: list[tuple[int, str, str]]
    pairs: set[tuple[int, str]] | None = None


def _first(mapping: dict, keys: tuple[str, ...]) -> tuple[Any, str]:
    """``(value, key)`` of the first present key, or ``(None, "")``."""
    for key in keys:
        if key in mapping and mapping[key] not in (None, ""):
            return mapping[key], key
    return None, ""


def _as_index(value: Any) -> int | None:
    """``value`` as an integer index, or ``None`` when it is not one."""
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    try:
        text = str(value).strip()
        return int(text) if text == str(int(text)) else None
    except (TypeError, ValueError):
        return None


def parse_centers(path: Path | None) -> Centers:
    """Audit columns from ``<stem>.centers.json``, plus the centres to sequence-check.

    The real schema (VERIFIED 2026-10-09) is::

        {"centers": [{"position": 10, "token": "HIS-S", "target_idx": 22,
                      "comparison_idx": 21, "placement_score": 0.708}],
         "free_positions": [9, 11, …], "center_mode": "his_s",
         "placement_region": "interface"}

    There is NO chain field, NO residue number and NO charge-field range in it, so
    ``charge_field_min``/``charge_field_max`` stay NA (the range is printed in the task
    ``.out`` only). ``placement_score`` IS there per centre and is collected as the
    MINIMUM over the centres. Reading stays defensive for other schemas: ``centers_source``
    records which key the centres came from and ``unparsed`` means nothing recognisable
    was there. Unparsed means the CENTRE columns are NA -- it says nothing about the
    metrics columns, which come from a different file.
    """
    out: dict[str, Any] = {
        "centers": "",
        "n_centers": pd.NA,
        "centers_source": "missing",
        "centers_resnum_key": "",
        "placement_metric": "",
        "redesign_scope": "",
        "placement_score_min": pd.NA,
        "charge_field_min": pd.NA,
        "charge_field_max": pd.NA,
        "n_redesign": pd.NA,
    }
    if path is None or not path.is_file():
        return Centers(out, [], None)
    try:
        obj = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        out["centers_source"] = "unreadable"
        return Centers(out, [], None)

    entries: list | None = None
    source = "unparsed"
    if isinstance(obj, list):
        entries, source = obj, "<root list>"
    elif isinstance(obj, dict):
        for key in _CENTER_LIST_KEYS:
            value = obj.get(key)
            if isinstance(value, list):
                entries, source = value, key
                break

    checks: list[tuple[int, str, str]] = []
    pairs: set[tuple[int, str]] = set()
    pairs_complete = True
    if entries is not None:
        labels, resnum_keys = [], []
        placements: list[float] = []
        for entry in entries:
            if isinstance(entry, dict):
                placement, _ = _first(entry, _PLACEMENT_KEYS)
                score = _num_or_na(placement) if placement is not None else pd.NA
                if score is not pd.NA:
                    placements.append(float(score)) # type: ignore
                chain, _ = _first(entry, _CHAIN_KEYS)
                resnum, resnum_key = _first(entry, _RESNUM_KEYS)
                state, _ = _first(entry, _STATE_KEYS)
                index = _as_index(resnum)
                token = _center_token(state)
                if index is not None and token is not None:
                    pairs.add((index, token.upper()))
                else:
                    # One centre this tool cannot express as (index, token) makes the
                    # SET incomparable, so the metrics cross-check must say NA rather
                    # than report a difference that is only its own blindness.
                    pairs_complete = False
                if resnum_key:
                    # Label the numbering convention IN the column, so a reader cannot
                    # mistake an index for a residue number by reading `centers` alone.
                    resnum_keys.append(
                        f"{resnum_key}_index0"
                        if resnum_key in _INDEX0_KEYS
                        else resnum_key
                    )
                    letter = residue_letter(token) if token is not None else None
                    if resnum_key in _INDEX0_KEYS and index is not None and letter:
                        checks.append((index, letter, str(token)))
                parts = [str(x) for x in (chain, resnum, state) if x is not None]
                labels.append(":".join(parts) if parts else json.dumps(entry, sort_keys=True))
            else:
                labels.append(str(entry))
                pairs_complete = False
        out["centers"] = ";".join(labels)
        out["n_centers"] = len(entries)
        out["centers_resnum_key"] = ",".join(sorted(set(resnum_keys)))
        # The MINIMUM, not the mean: one barely-wanted centre is enough to make the
        # whole placement doubtful, and a strong second centre must not average it
        # away. NA -- never 0 -- when no centre recorded the key: 0 is a real score
        # (a centre the model was indifferent to) and must stay distinguishable from
        # "this file did not say".
        out["placement_score_min"] = min(placements) if placements else pd.NA
    out["centers_source"] = source

    if isinstance(obj, dict):
        # READ FROM SOURCE at 7ffb756 (design_ph.py:2001-2005). placement_metric
        # decides WHICH positions become centres, and its default changed between
        # commits -- ddE at 7ffb756, dE_complex at dee5a4d, where the key does not
        # exist at all. The same structure therefore gets different centres under the
        # two builds, so the table has to say which. Empty (NA) on older output.
        for key in ("placement_metric", "redesign_scope"):
            value = obj.get(key)
            if isinstance(value, str) and value:
                out[key] = value
        lo, hi = _charge_field_range(obj)
        out["charge_field_min"], out["charge_field_max"] = lo, hi
        for key in _REDESIGN_KEYS:
            value = obj.get(key)
            if isinstance(value, bool):
                continue
            if isinstance(value, int):
                out["n_redesign"] = value
                break
            if isinstance(value, list):
                out["n_redesign"] = len(value)
                break
    return Centers(out, checks, pairs if (pairs and pairs_complete) else None)


def check_centers_against_sequence(
    checks: list[tuple[int, str, str]], sequence: str
) -> tuple[Any, str]:
    """Does the designed sequence carry each centre's residue at the centre's index?

    This is the ONE trust check that can catch a wrong index base, and it is why the
    tool may print a centre position at all. A disagreement is never a column value to
    read past: the caller fails the parent.

    Returns ``(True, "")`` when every checkable centre agreed, ``(False, why)`` on a
    disagreement, and ``(pd.NA, "")`` when the check could not run -- no centres, no
    recognisable residue in any token, or a multi-chain ``sequence`` (``--output-scope
    all`` joins chains with '/', so a binder-chain index no longer indexes it).
    """
    if not checks or not sequence or "/" in sequence:
        return pd.NA, ""
    for index, letter, token in checks:
        if not 0 <= index < len(sequence):
            return False, (
                f"centre {token!r} is at index {index}, outside the "
                f"{len(sequence)}-residue designed sequence"
            )
        if sequence[index] != letter:
            return False, (
                f"centre {token!r} is at index {index}, where the designed sequence "
                f"has {sequence[index]!r} and not the {letter!r} the token names"
            )
    return True, ""


def check_centers_against_metrics(
    pairs: set[tuple[int, str]] | None, cell: Any
) -> tuple[Any, str]:
    """Does this row's own ``centers`` cell name the same centres as ``centers.json``?

    ``centers.json`` is written ONCE per structure and stamped on every child row;
    ``metrics_centers`` is written per DESIGN. Nothing checked they agreed until now,
    and with ``--centers-per-type 2`` there is a real set to disagree about.

    Returns ``(True, "")`` on an exact set match, ``(False, why)`` on a clear
    disagreement, and ``(pd.NA, "")`` whenever either side could not be read as
    ``(index, token)`` pairs -- an unverified multi-centre separator must produce NA,
    never a failure.
    """
    if not pairs:
        return pd.NA, ""
    other = _parse_metrics_centers(cell)
    if other is None:
        return pd.NA, ""
    if other == pairs:
        return True, ""
    only_json = sorted(f"{i}:{t}" for i, t in pairs - other)
    only_csv = sorted(f"{i}:{t}" for i, t in other - pairs)
    return False, (
        f"centers.json names {{{', '.join(sorted(f'{i}:{t}' for i, t in pairs))}}} but "
        f"this row's metrics centers cell names "
        f"{{{', '.join(sorted(f'{i}:{t}' for i, t in other))}}}"
        + (f"; only in centers.json: {', '.join(only_json)}" if only_json else "")
        + (f"; only in the metrics row: {', '.join(only_csv)}" if only_csv else "")
    )


def check_sequence_ext(ext_seq: str, sequence: str) -> tuple[bool, str]:
    """Is ``sequence_ext`` the same molecule as ``sequence``, token for token?

    ``sequence_ext`` is the design with explicit protonation tokens, and it is as
    mappable-wrong as any metric: a correct-looking token string belonging to another
    draw would be read as this design's protonation pattern. Two assertions, both cheap:
    one token per residue, and every token that resolves to a residue must BE that
    residue. A token this tool cannot resolve (an unknown 3-letter head) is skipped --
    unrecognised is not disagreement.
    """
    tokens = ext_seq.split()
    if not tokens:
        return False, "the extended sequence is empty"
    if len(tokens) != len(sequence):
        return False, (
            f"{len(tokens)} extended token(s) against {len(sequence)} residue(s) in the "
            f"canonical sequence"
        )
    for i, token in enumerate(tokens):
        letter = residue_letter(token)
        if letter is None:
            letter = token.strip().upper() if len(token.strip()) == 1 else None
        if letter is not None and letter != sequence[i]:
            return False, (
                f"token {i} is {token!r} ({letter}) but the canonical sequence has "
                f"{sequence[i]!r} there"
            )
    return True, ""


def split_per_centre(cell: Any) -> list[str]:
    """The per-centre fields of a ``;``-joined metrics cell, in file order."""
    text = "" if cell is None else str(cell).strip()
    return [x.strip() for x in text.split(";")] if text else []


def _floats(parts: list[str]) -> list[float] | None:
    """``parts`` as floats, or ``None`` if any is not a finite number.

    ``nan`` is deliberately rejected: design_ph writes ``nan`` per centre for the apo
    side when there are no target chains, and a mean over it is not a measurement.
    """
    out: list[float] = []
    for part in parts:
        try:
            value = float(part)
        except (TypeError, ValueError):
            return None
        if value != value or value in (float("inf"), float("-inf")):
            return None
        out.append(value)
    return out


class PerCentre(NamedTuple):
    """Per-centre energies for ONE design row, and whether they are self-consistent.

    ``columns`` are the values to stamp; ``note`` is non-empty when the row contradicts
    itself, in which case the caller fails the parent rather than publishing a mean
    whose parts do not add up.
    """

    columns: dict[str, Any]
    note: str = ""


def per_centre_energies(
    row: dict[str, str], lower_to_col: dict[str, str], n_centers: Any
) -> PerCentre:
    """Split ``dE_complex`` / ``dE_apo`` into numeric + verbatim + derived columns.

    READ FROM SOURCE at 7ffb756 (design_ph.py:1859-1868): both cells are ``;``-joined
    4-dp per-centre values, ``dE_complex_mean``/``dE_apo_mean`` are the numeric means,
    and ``ddE_mean = mean(dE_complex) - mean(dE_apo)``, which is also the mean of the
    per-centre ddE. Not yet observed in a run -- the only real output had one centre.

    ``ddE_min`` is the point of this function: ``ddE_mean`` cannot distinguish two
    histidines that each carry half the switch from one that carries all of it and one
    that carries none. The minimum over centres can.

    Two contradictions FAIL the parent (``note``):
      * the per-centre counts disagree with each other or with the row's ``centers``
        cell -- the values cannot be attributed to centres at all;
      * ``mean(per-centre ddE)`` differs from the file's ``ddE_mean`` by more than
        2e-3, which 4-dp rounding cannot explain.
    """
    out: dict[str, Any] = {
        "dE_complex_centres": "",
        "dE_apo_centres": "",
        "ddE_centres": "",
        "ddE_min": pd.NA,
    }
    raw: dict[str, list[str]] = {}
    for name in PER_CENTRE_COLUMNS:
        col = lower_to_col.get(name.casefold())
        parts = split_per_centre(row.get(col, "")) if col else []
        raw[name] = parts
        out[f"{name}_centres"] = ";".join(parts)

    complex_parts, apo_parts = raw["dE_complex"], raw["dE_apo"]
    if not complex_parts or not apo_parts:
        return PerCentre(out)

    if len(complex_parts) != len(apo_parts):
        return PerCentre(out, (
            f"{len(complex_parts)} per-centre dE_complex value(s) against "
            f"{len(apo_parts)} dE_apo value(s)"
        ))
    # The row's OWN centres cell, not centers.json: this is a per-row consistency
    # check, and centers.json is compared separately by centers_match_metrics.
    centres_col = lower_to_col.get("centers")
    centres_cell = row.get(centres_col, "") if centres_col else ""
    cell_centres = split_per_centre(centres_cell)
    if cell_centres and len(cell_centres) != len(complex_parts):
        return PerCentre(out, (
            f"{len(complex_parts)} per-centre energy value(s) against "
            f"{len(cell_centres)} centre(s) in this row's centers cell "
            f"({centres_cell!r})"
        ))
    if not isinstance(n_centers, int) or n_centers <= 0:
        pass
    elif len(complex_parts) != n_centers:
        return PerCentre(out, (
            f"{len(complex_parts)} per-centre energy value(s) against "
            f"{n_centers} centre(s) in centers.json"
        ))

    complex_values = _floats(complex_parts)
    apo_values = _floats(apo_parts)
    if complex_values is None or apo_values is None:
        # nan apo (no target chains) or an unparseable cell: the verbatim strings are
        # still collected, the derived numbers stay NA, and nothing is contradicted.
        return PerCentre(out)

    ddes = [c - a for c, a in zip(complex_values, apo_values)]
    out["ddE_centres"] = ";".join(f"{v:.4f}" for v in ddes)
    out["ddE_min"] = min(ddes)

    mean_col = lower_to_col.get("dde_mean")
    stated = row.get(mean_col, "").strip() if mean_col else ""
    if stated:
        try:
            value = float(stated)
        except ValueError:
            return PerCentre(out)
        if value == value:  # not nan
            got = sum(ddes) / len(ddes)
            # 2e-3: the parts are rounded to 4 dp before we see them, the mean is not.
            if abs(got - value) > 2e-3:
                return PerCentre(out, (
                    f"the per-centre ddE values average {got:.4f} but the file's "
                    f"ddE_mean is {value:.4f}; the columns do not describe the same "
                    f"design"
                ))
    return PerCentre(out)


def hamming(a: str, b: str) -> Any:
    """Per-position differences between two equal-length strings, else NA.

    NA on a length mismatch rather than a truncated count: two sequences of different
    length are not two versions of one design, and a number here would be read as one.
    """
    if not a or not b or len(a) != len(b):
        return pd.NA
    return sum(1 for x, y in zip(a, b) if x != y)


def _charge_field_range(obj: dict) -> tuple[Any, Any]:
    """The charge-field range, in whichever of three shapes the file uses, else NA/NA."""
    if "charge_field_min" in obj or "charge_field_max" in obj:
        return _num_or_na(obj.get("charge_field_min")), _num_or_na(obj.get("charge_field_max"))
    for key in ("charge_field", "charge_field_range", "charge_range"):
        value = obj.get(key)
        if isinstance(value, (list, tuple)) and len(value) == 2:
            return _num_or_na(value[0]), _num_or_na(value[1])
        if isinstance(value, dict):
            return _num_or_na(value.get("min")), _num_or_na(value.get("max"))
    return pd.NA, pd.NA


def _num_or_na(value: Any) -> Any:
    if isinstance(value, bool) or value is None:
        return pd.NA
    if isinstance(value, (int, float)):
        return value
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return pd.NA


# --------------------------------------------------------------------------- #
# the collector
# --------------------------------------------------------------------------- #


def read_meta(out_dir: Path) -> dict:
    """The run's sidecar, or ``{}`` for a run_dir written without one."""
    path = out_dir / META_FILENAME
    if not path.is_file():
        return {}
    try:
        meta = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return meta if isinstance(meta, dict) else {}


# Flags recoverable from the sub-manifest's argv tail, keyed by the column they fill.
# ``--lam`` is not here: it travels as its own TSV field, not in the tail.
_TAIL_FLAGS = {"model_noise": "--checkpoint_model_noise",
               "protect_titratable": "--protect_titratable",
               "apo_weight": "--apo_weight"}


def noise_and_lam_from_tasks(out_dir: Path) -> tuple[str, str]:
    """``(model_noise, lam_sweep)`` recovered from the sub-manifests, else ``("", "")``.

    Thin wrapper over ``run_params_from_tasks`` kept for readability at the call site.

    The fallback for a leaf whose ``.meta.json`` predates those keys. ``run_atomiumph``
    writes one TSV row per subgroup into ``atomiumph_tasks/task_<t>.tsv``::

        <grp_dir>\\t<binder chains>\\t<target chains>\\t<lam …>\\t<design_ph argv tail>

    and the argv tail carries ``--checkpoint_model_noise <level>``. Values that are not
    identical across every row are reported as unknown (``""``) rather than attributed
    to rows they may not describe -- the leaf was reused by two different submits.
    Never a guessed default: the whole point of the column is that n03 and n05 pick
    DIFFERENT centres, so a wrong label is worse than none.
    """
    noises: set[str] = set()
    lams: set[str] = set()
    for task_file in sorted((out_dir / "atomiumph_tasks").glob("task_*.tsv")):
        try:
            text = task_file.read_text()
        except OSError:
            continue
        for line in text.splitlines():
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            lams.add(fields[3].strip())
            tail = fields[4].split()
            for i, token in enumerate(tail):
                if token == "--checkpoint_model_noise" and i + 1 < len(tail):
                    noises.add(tail[i + 1])
                elif token.startswith("--checkpoint_model_noise="):
                    noises.add(token.split("=", 1)[1])
    noise = noises.pop() if len(noises) == 1 else ""
    lam = lams.pop() if len(lams) == 1 else ""
    return noise, lam


def run_params_from_tasks(out_dir: Path) -> dict[str, str]:
    """Every recoverable run parameter from the sub-manifests: ``{column: value}``.

    The fallback for a leaf whose ``.meta.json`` predates a key. A value that is not
    IDENTICAL across every sub-manifest row is omitted rather than attributed to rows it
    may not describe. Never a guessed default -- an absent key must read as NA, because
    the whole point of these columns is that two values of them are not comparable.
    """
    found: dict[str, set[str]] = {k: set() for k in _TAIL_FLAGS}
    lams: set[str] = set()
    for task_file in sorted((out_dir / "atomiumph_tasks").glob("task_*.tsv")):
        try:
            text = task_file.read_text()
        except OSError:
            continue
        for line in text.splitlines():
            fields = line.split("\t")
            if len(fields) < 5 or not line.strip():
                continue
            lams.add(fields[3].strip())
            tail = fields[4].split()
            for column, flag in _TAIL_FLAGS.items():
                for i, token in enumerate(tail):
                    if token == flag and i + 1 < len(tail):
                        found[column].add(tail[i + 1])
                    elif token.startswith(f"{flag}="):
                        found[column].add(token.split("=", 1)[1])
    out = {k: (v.pop() if len(v) == 1 else "") for k, v in found.items()}
    out["lam_sweep"] = lams.pop() if len(lams) == 1 else ""
    return out


def _error_row(parent: str, message: str, data: dict[str, Any]) -> Collected:
    """One child row recording that this parent's output could not be trusted.

    ``path=""`` on purpose: a <leaf>_path pointing at a .metrics.csv would be handed to
    a sequence tool by a later ``-i atomiumph_path``.
    """
    return Collected(name=f"{parent}_p0", parent=parent, path="",
                     status=f"error: {message}", data=data)


def collect_atomiumph(ctx: CollectCtx) -> CollectEach:
    """Per-parent collector. atomiumph is a create tool: the framework iterates the
    ready parents and this mints one child row (``<parent>_p<i>``) per designed
    sequence, carrying ``parent`` for lineage. Rows are rebuilt from disk, so
    re-running is idempotent."""
    index = index_outputs(ctx.out_dir)
    meta = read_meta(ctx.out_dir)
    expected = meta.get(EXPECTED_KEY)
    expected = expected if isinstance(expected, int) else pd.NA
    run_columns = {
        "binder_chains": str(meta.get(BINDER_CHAINS_KEY, "")),
        "target_chains": str(meta.get(TARGET_CHAINS_KEY, "")),
        # Run-wide, and both decide how the row must be READ: center_mode fixes what a
        # high ddE_mean means (which pH binds), output_scope 'all' makes `sequence` a
        # slash-joined multi-chain string so seq_len stops being a binder length.
        "center_mode": str(meta.get(CENTER_MODE_KEY, "")),
        "output_scope": str(meta.get(OUTPUT_SCOPE_KEY, "")),
    }
    # The noise level decides WHICH CENTRES design_ph picks, so two runs over the same
    # structures at n03 and n05 are not comparable. The run records it; a leaf written
    # before that falls back to the sub-manifest, and to "" (NA) if even that is gone.
    sidecar = {
        "model_noise": str(meta.get(MODEL_NOISE_KEY, "") or ""),
        "lam_sweep": str(meta.get(LAM_SWEEP_KEY, "") or ""),
        "protect_titratable": str(meta.get(PROTECT_TITRATABLE_KEY, "") or ""),
        "apo_weight": str(meta.get(APO_WEIGHT_KEY, "") or ""),
    }
    if not all(sidecar.values()):
        from_tasks = run_params_from_tasks(ctx.out_dir)
        for key, value in sidecar.items():
            sidecar[key] = value or from_tasks.get(key, "")
    run_columns.update(sidecar)
    warned: set[str] = set()
    stale_warned: set[str] = set()

    # F4, THE COLLECTOR HALF -- ADVISORY, NOT FATAL. See the long note below.
    #
    # The hazard is real: `grp_<gid>` restarts at 0 on every submit, nothing prunes
    # `grp_<g>/out`, and `.meta.json` is rewritten each submit, so a rerun under the
    # same dir label with a different --model-noise / --center-mode / chains /
    # --protect-titratable leaves the PREVIOUS run's outputs in place to be stamped with
    # the NEW run's parameters. That is prevented at SUBMIT time, deterministically, by
    # `_check_no_stale_outputs` in run_atomiumph.py -- which is the half that can refuse
    # safely, because it knows exactly which stems this submit will regenerate.
    #
    # Here, all that is available is a FILE MTIME against the sidecar's, and that is not
    # sound enough to fail a parent on:
    #   * a PARTIAL `--force` rerun rewrites the sidecar and regenerates only some
    #     stems; the untouched stems are then "older" though their rows were correct,
    #     and erroring them would destroy good rows in an existing table;
    #   * any copy of a run_dir (rsync, a Volume download, `cp -r` without -p) can
    #     reorder mtimes without a single byte changing. MEASURED while building this:
    #     a plain `shutil.copytree` of a healthy leaf flipped every stem to "stale".
    # So a false positive is cheap to produce and expensive to suffer, while a false
    # NEGATIVE is already covered at submit time. It is recorded as a trust column
    # (`output_predates_meta`) and printed once, and a caller who wants it fatal writes
    # one line in a `-f` module.
    meta_mtime = 0.0
    meta_path = ctx.out_dir / META_FILENAME
    if meta_path.is_file():
        try:
            meta_mtime = meta_path.stat().st_mtime
        except OSError:
            meta_mtime = 0.0

    def stale_files(files: dict[str, list[Path]]) -> list[Path]:
        """Output files older than the sidecar, i.e. left by an earlier submit."""
        if not meta_mtime:
            return []
        out: list[Path] = []
        for paths in files.values():
            for path in paths:
                try:
                    if path.stat().st_mtime < meta_mtime:
                        out.append(path)
                except OSError:
                    continue
        return sorted(out)

    # design_ph.py is always run in --pdb_dir (batch) mode, where a structure that fails
    # does NOT abort the batch: the task exits 0 and that stem simply has no files. One
    # summary line here, rather than a print buried per design, because the count is the
    # signal -- the per-structure reason is in design_ph's trailing failure list in the
    # task's .out.
    absent = [str(n) for n in ctx.ready.index if str(n) not in index]
    if absent:
        shown = ", ".join(absent[:20]) + (" …" if len(absent) > 20 else "")
        print(
            f"atomiumph: {len(absent)} of {len(ctx.ready.index)} parent(s) produced no "
            f"design_ph output and collect no rows: {shown}\n"
            f"  (batch mode does not abort on a per-structure failure -- read the "
            f"failure list at the end of the task's .out)"
        )

    def one(d: DesignCtx) -> Iterable[Collected]:
        files = index.get(d.name)
        if not files:
            # Reported in the summary above. Nothing is minted and the framework leaves
            # the parent uncollected -- the same contract as atomium, and the right one,
            # because collect's ready set can legitimately contain designs a filtered
            # run never submitted.
            return

        base: dict[str, Any] = {
            "join_mode": "",
            "join_key": "",
            "join_checked_fields": "",
            "ext_join_mode": "",
            "ext_checked_fields": "",
            "n_metrics_rows": pd.NA,
            "n_fa_entries": pd.NA,
            "n_ext_entries": pd.NA,
            "n_expected": expected,
            # F10: NA, not 0. On an error row that never parsed the FASTA, zero
            # preamble entries were not OBSERVED -- the file was not read. 0 here would
            # be a measurement, and would make `n_fa_entries + n_preamble_dropped`
            # arithmetic silently wrong.
            "n_preamble_dropped": pd.NA,
            # NA, never False: an error row that never got as far as the centres has
            # not disagreed with anything, it has simply not been checked.
            "centers_seq_check": pd.NA,
            "centers_match_metrics": pd.NA,
            # False = this stem's files are newer than the sidecar, i.e. normal.
            "output_predates_meta": False,
            "native_sequence": "",
            "n_diff_native": pd.NA,
            **run_columns,
            "image_commit": "",
        }

        stale = stale_files(files)
        if stale:
            base["output_predates_meta"] = True
            if d.name not in stale_warned:
                stale_warned.add(d.name)
                print(
                    f"atomiumph: {d.name}: {len(stale)} output file(s) are OLDER than "
                    f"the leaf's {META_FILENAME}. Either this leaf was submitted twice "
                    f"under one dir label and these are the earlier run's results "
                    f"(in which case the run columns on these rows are WRONG), or the "
                    f"run_dir was copied in a way that did not preserve mtimes. "
                    f"Check atomiumph_output_predates_meta before using these rows."
                )

        for kind in ("metrics", "fa", "ext_fa", "centers", "couplings"):
            if len(files.get(kind, [])) > 1:
                where = ", ".join(str(p) for p in files[kind])
                yield _error_row(
                    d.name, f"ambiguous output: {len(files[kind])} {kind} files ({where})",
                    {**base, "join_mode": "failed"},
                )
                return

        metrics_path = files.get("metrics", [None])[0]
        fasta_path = files.get("fa", [None])[0]
        if metrics_path is None or fasta_path is None:
            missing = ".metrics.csv" if metrics_path is None else ".fa"
            # Any file of this stem anchors the walk up to image_commit.txt equally
            # well, and the commit is a trust column worth keeping even on an error
            # row. Both of these can be absent at once, so there may be no anchor.
            anchor = next((p for paths in files.values() for p in paths), None)
            base["image_commit"] = commit_for(anchor, ctx.out_dir) if anchor else ""
            yield _error_row(
                d.name, f"incomplete output: {d.name}{missing} is absent",
                {**base, "join_mode": "failed"},
            )
            return

        base["image_commit"] = commit_for(metrics_path, ctx.out_dir)
        columns, rows, n_overflow = read_metrics(metrics_path)
        # The stem IS the parent name, and design_ph.py writes the native echo as
        # ``><stem>_native`` -- so the stem is what makes the preamble recognisable.
        all_entries = parse_fasta(fasta_path)
        entries, dropped = drop_preamble(all_entries, d.name)
        # F9: KEEP what was dropped. The echo is the sequence every design was built
        # from, and it is the ONLY way to check `n_mutations` against the sequences
        # themselves -- the real data has mutations=4 on rows whose canonical distance
        # to this sequence is 3, and without this column that cannot even be noticed.
        # The native REFERENCE is specifically the ``_native`` echo, never a stage-1
        # seed: under --seed-mode generate the seeds are what the design was built
        # FROM, but the echo is the input structure's own sequence and is the only
        # stable baseline for n_diff_native across seed modes.
        native = next(
            (
                seq
                for head, seq in all_entries[: max(dropped, 0)]
                if first_token(head).casefold()
                in (f"{d.name.casefold()}_native", "native")
            ),
            "",
        )
        # ``join=" "``: .ext.fa is space-separated tokens, so wrapped continuation
        # lines must not be fused (see parse_fasta).
        ext_entries, _ = drop_preamble(
            parse_fasta(files["ext_fa"][0], join=" ") if files.get("ext_fa") else [],
            d.name,
        )
        base.update(
            n_metrics_rows=len(rows),
            n_fa_entries=len(entries),
            n_ext_entries=len(ext_entries),
            n_preamble_dropped=dropped,
            native_sequence=native,
        )

        if n_overflow:
            yield _error_row(
                d.name,
                f"malformed metrics: {n_overflow} row(s) carry more fields than the "
                f"{len(columns)}-column header, so their values are shifted",
                {**base, "join_mode": "failed"},
            )
            return

        join = align_entries(entries, columns, rows)
        if join.mode == "failed":
            # The one failure this tool exists to never paper over: the sequences and
            # the numbers cannot be put in the same row. An error row, not N rows of
            # plausible numbers attached to the wrong sequences.
            yield _error_row(
                d.name, f"cannot map sequences to metrics: {join.note}",
                {**base, "join_mode": "failed",
                 "join_checked_fields": ",".join(join.checked)},
            )
            return

        ext = align_ext(entries, ext_entries)
        base.update(
            join_mode=join.mode,
            join_key=",".join(join.keys),
            join_checked_fields=",".join(join.checked),
            ext_join_mode=ext.mode,
            ext_checked_fields=",".join(ext.checked),
        )

        centers = parse_centers(files.get("centers", [None])[0])
        couplings = files.get("couplings", [None])[0]

        # THE CENTRE TRUST CHECK. `position` is a 0-based index into the binder chain,
        # so the designed sequence must carry the residue the token names at exactly
        # that index. A disagreement means the index base (or the chain) is wrong, and
        # every centre column would then point a reader at the wrong residue -- so the
        # parent FAILS instead of shipping numbers nobody can place.
        seq_checks = [
            check_centers_against_sequence(centers.checks, sequence)
            for (_, sequence), _ in join.pairs
        ]
        bad = next(((v, why) for v, why in seq_checks if v is False), None)
        if bad is not None:
            yield _error_row(
                d.name,
                f"centres do not match the designed sequence: {bad[1]}. The centre "
                f"index base (centers_resnum_key="
                f"{centers.columns['centers_resnum_key']!r}) or the designed chain is "
                f"wrong, so no centre column can be trusted",
                {**base, "join_mode": join.mode,
                 "join_key": ",".join(join.keys),
                 "join_checked_fields": ",".join(join.checked),
                 **centers.columns, "centers_seq_check": False},
            )
            return

        # F2, THE SECOND CENTRE CHECK. centers.json is written ONCE per structure and
        # stamped on every child row; each metrics row carries its OWN `centers` cell.
        # Nothing compared them until now, and with --centers-per-type 2 there is a real
        # set to disagree about. A clear disagreement means the per-parent centre
        # columns describe a different design from the one in the row, so the parent
        # fails exactly as a sequence-check failure does. An unreadable cell (the
        # multi-centre separator is unverified) is NA, never a failure.
        lower_cols = {c.casefold(): c for c in columns}

        # The per-centre energy split and its two self-consistency checks. Done before
        # any row is yielded so a contradiction fails the parent rather than publishing
        # a ddE_mean whose parts do not add up to it.
        per_centre = [
            per_centre_energies(row, lower_cols, centers.columns.get("n_centers"))
            for _, row in join.pairs
        ]
        bad_pc = next((x for x in per_centre if x.note), None)
        if bad_pc is not None:
            yield _error_row(
                d.name,
                f"per-centre energies are inconsistent: {bad_pc.note}",
                {**base, "join_mode": join.mode,
                 "join_key": ",".join(join.keys),
                 "join_checked_fields": ",".join(join.checked),
                 **centers.columns},
            )
            return

        metrics_centers_col = lower_cols.get("centers")
        match_checks = [
            check_centers_against_metrics(
                centers.pairs,
                row.get(metrics_centers_col, "") if metrics_centers_col else None,
            )
            for _, row in join.pairs
        ]
        bad_match = next(((v, why) for v, why in match_checks if v is False), None)
        if bad_match is not None:
            yield _error_row(
                d.name,
                f"centres disagree between centers.json and the metrics rows: "
                f"{bad_match[1]}. The per-parent centre columns would describe a "
                f"different design from the row they sit on",
                {**base, "join_mode": join.mode,
                 "join_key": ",".join(join.keys),
                 "join_checked_fields": ",".join(join.checked),
                 **centers.columns, "centers_seq_check": seq_checks[0][0],
                 "centers_match_metrics": False},
            )
            return

        # The extended FASTA is ONE column's provenance, so it gets its own status. A
        # row's <leaf>_status must gate the metric it describes, not every metric in the
        # row -- an absent .ext.fa must not hide a perfectly collected ddE_mean.
        ext_status = {
            "missing": f"error: {d.name}.ext.fa is absent",
            "failed": f"error: cannot map sequence_ext: {ext.note}",
        }.get(ext.mode, "OK")

        for i, ((header, sequence), row) in enumerate(join.pairs):
            data: dict[str, Any] = dict(base)
            data["entry_index"] = i
            data["sequence"] = sequence
            data["seq_len"] = len(sequence)
            # F9: how far this design moved from the echoed input, computed from the
            # sequences themselves. NA on a length mismatch or an absent echo -- never
            # a truncated count. This is what makes the header's `mutations=` (and so
            # the metrics `n_mutations`) checkable by a reader.
            data["n_diff_native"] = hamming(sequence, native)
            data["centers_match_metrics"] = match_checks[i][0]
            # F3: sequence_ext must be THIS design, token for token. A correct-looking
            # token string belonging to another draw is exactly as wrong as a
            # mis-mapped ddE_mean, and nothing downstream could detect it.
            ext_seq = ext.seqs[i] if ext.seqs else ""
            row_ext_status = ext_status
            if ext_seq:
                ok_ext, why = check_sequence_ext(ext_seq, sequence)
                if not ok_ext:
                    ext_seq = ""
                    row_ext_status = (
                        f"error: sequence_ext does not match the canonical sequence: "
                        f"{why}"
                    )
            data["ext_status"] = row_ext_status
            data["sequence_ext"] = ext_seq
            # Metrics columns pass through VERBATIM under their upstream header names,
            # so a column upstream adds appears here without a code change and nothing
            # loses its provenance to a rename.
            for column in columns:
                key = column
                if key in RESERVED:
                    # VERIFIED collisions: `name` (the design id -- also the join key,
                    # kept away from the table's own index) and `centers` (the metrics
                    # file's own "10:HIS-S" string, kept away from this collector's
                    # parse of centers.json). Both become metrics_<name>; neither is
                    # announced, because the rename is the normal case, not news.
                    key = f"metrics_{column}"
                    if column not in warned and column not in EXPECTED_COLLISIONS:
                        warned.add(column)
                        print(
                            f"atomiumph: metrics column {column!r} collides with a "
                            f"trust column; collected as {key!r}"
                        )
                if column.casefold() in {c.casefold() for c in PER_CENTRE_COLUMNS}:
                    # Keep the column NUMERIC. One centre -> the single value; several
                    # -> NA, with the verbatim string in <name>_centres. Never a
                    # column that is a float on some rows and "a;b" on others.
                    parts = split_per_centre(row.get(column, ""))
                    data[key] = _as_number(parts[0]) if len(parts) == 1 else pd.NA
                else:
                    data[key] = _as_number(row.get(column, ""))
            data.update(per_centre[i].columns)
            data.update(centers.columns)
            data["centers_seq_check"] = seq_checks[i][0]
            data["metrics_path"] = str(metrics_path)
            data["couplings_path"] = str(couplings) if couplings else ""
            data["header"] = header
            yield Collected(
                name=f"{d.name}_p{i}",
                parent=d.name,
                path=fasta_path,
                status="OK",
                data=data,
            )

    return one
