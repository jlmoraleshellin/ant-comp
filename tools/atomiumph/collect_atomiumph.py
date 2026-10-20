"""Collect ``design_ph.py`` outputs into a child table, one row per designed sequence.

The submitter symlinks each input as ``grp_<g>/inputs/<design_name>.pdb`` and
``design_ph.py`` names every output after that stem, so the stem IS the parent row::

    <out_dir>/**/<stem>.metrics.csv, .fa, .ext.fa, .centers.json, .couplings.csv

The search is a recursive glob (the layout under ``--out_folder`` is undocumented); a
stem found in two places is an error, not a first match. ``summary.csv`` is ignored and
``couplings.csv`` stays a file (``couplings_path``), not columns.

THE JOIN IS THE RISK. Sequences come from ``<stem>.fa`` and numbers from
``<stem>.metrics.csv``; if their orders diverge, every design silently inherits
another's ``ddE_mean``. ``join_mode`` says how they were paired:

  1. ``key`` on the design id (header first token == metrics ``name``), or on
     ``(id, lambda)`` when ids repeat; ``lambda`` and ``mutations`` are then
     cross-checked (``join_checked_fields``). The rounded echoes (``dE_complex``,
     ``ddE``) are never compared: the header prints 3 decimals, the CSV 4.
  2. ``key`` on shared identity ``key=value`` header fields (headers without an id).
  3. ``positional``, only when no key exists: counts equal, shared fields compared.
  4. ``failed``: counts differ, key value sets disagree, a cross-checked field
     disagrees, a file is missing or ambiguous, a metrics row overflows its header, or
     a centre check fails. The parent then gets ONE ``error`` row with the evidence,
     never numbers it cannot place.

``dE_complex`` / ``dE_apo`` are ``;``-joined per-centre strings: the numeric column
keeps a single value (NA with several centres) and ``<name>_centres`` the verbatim
string; ``ddE_centres`` and ``ddE_min`` derive from them (design_ph.py:1859-1877).

Under ``--binder-alone`` (sidecar ``atomiumph_binder_alone``) every row has
``binder_alone = True`` and the ``TARGET_DEPENDENT_NUMERIC`` columns are NA. That is an
allowlist of names, so re-read the metrics header on a new pin. Semantics and traps:
``.claude/skills/atomiumph/SKILL.md``.
"""

import csv
import json
import re
from pathlib import Path
from typing import Any, Iterable, NamedTuple

import pandas as pd

from prosapia.core import Collected, CollectCtx, CollectEach, DesignCtx

# Sidecar written by the run (see run_atomiumph.py); the chains are run-wide.
META_FILENAME = ".meta.json"
EXPECTED_KEY = "atomiumph_expected_per_structure"
BINDER_CHAINS_KEY = "atomiumph_binder_chains"
TARGET_CHAINS_KEY = "atomiumph_target_chains"
CENTER_MODE_KEY = "atomiumph_center_mode"
OUTPUT_SCOPE_KEY = "atomiumph_output_scope"
# Model noise decides which centres are picked; older leaves recover it from the tasks.
MODEL_NOISE_KEY = "atomiumph_model_noise"
LAM_SWEEP_KEY = "atomiumph_lam_sweep"
# Changes which positions may move, so two values are not like-for-like.
PROTECT_TITRATABLE_KEY = "atomiumph_protect_titratable"
# Modelled because the upstream client defaults it to 0.5 and design_ph.py to 0.0.
APO_WEIGHT_KEY = "atomiumph_apo_weight"
# The value design_ph was GIVEN (--binder-alone forces 0.0), not the one typed.
CHARGE_PRIOR_KEY = "atomiumph_charge_prior"
# True = no target in the graph: target-dependent numbers are undefined, not measured.
BINDER_ALONE_KEY = "atomiumph_binder_alone"

# Target-dependent columns, blanked under --binder-alone (NA numeric, "" text).
TARGET_DEPENDENT_NUMERIC = (
    "dE_apo",
    "dE_apo_mean",
    "ddE",
    "ddE_mean",
    "ddE_min",
    "interface_fraction",
    "charge_field_min",
    "charge_field_max",
)
TARGET_DEPENDENT_TEXT = ("dE_apo_centres", "ddE_centres")

# Header fields that identify a draw; a shared measurement would contradict a pairing.
IDENTITY_FIELDS = frozenset(
    {"lambda", "lam", "sample", "seed", "index", "design", "id", "name"}
)

# Written per subgroup by atomiumph.sh from the image's baked /opt/atomium_ph/COMMIT.
COMMIT_FILENAME = "image_commit.txt"

# Per-centre ``;``-joined strings (design_ph.py:1859-1860); see per_centre_energies.
PER_CENTRE_COLUMNS = ("dE_complex", "dE_apo")

# Key csv.DictReader files surplus fields under; such a row is shifted, hence an error.
OVERFLOW_KEY = "__overflow__"

# Suffixes design_ph.py appends to the stem; ".ext.fa" must be tried before ".fa".
SUFFIXES = (
    ("metrics", ".metrics.csv"),
    ("couplings", ".couplings.csv"),
    ("centers", ".centers.json"),
    ("ext_fa", ".ext.fa"),
    ("fa", ".fa"),
)

# Directories the run creates that never hold design_ph output.
SKIP_DIRS = ("atomiumph_logs", "atomiumph_tasks", "inputs")

# Preamble (input echo) headers, matched exactly: a substring would drop a real design.
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

# Names this collector owns; a colliding metrics column becomes metrics_<name>.
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
    "binder_alone",
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
    # design_ph echoes apo_weight too (design_ph.py:1869); kept as metrics_apo_weight.
    "apo_weight",
    # Same reasoning as apo_weight; --binder-alone makes it mode-dependent.
    "charge_prior",
    "metrics_path",
    "couplings_path",
    "header",
    # metrics.csv carries `name` (design id) and `centers` ("10:HIS-S"): both renamed.
    "name",
}

# Known collisions, renamed silently; any other collision prints a warning.
EXPECTED_COLLISIONS = {"name", "centers", "apo_weight"}


# --- file discovery ---


def index_outputs(out_dir: Path) -> dict[str, dict[str, list[Path]]]:
    """Map ``stem -> kind -> [paths]`` for every design_ph output under ``out_dir``.

    Recursive (the layout under ``--out_folder`` is undocumented), and lists so an
    ambiguous stem can be reported instead of resolved by luck.
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
    """The image commit recorded beside this output, walking up to ``out_dir``."""
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


# --- parsing ---


def parse_fasta(path: Path, join: str = "") -> list[tuple[str, str]]:
    """``[(header, sequence)]`` in file order; the header excludes the leading '>'.

    ``join`` glues a wrapped sequence's lines: ``""`` for canonical FASTA, ``" "`` for
    ``.ext.fa``, whose tokens would otherwise fuse (``GLU-DHIS-S``).
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
    """``key=value`` fields of a comma-separated FASTA header.

    A bare ``<word>_<int>`` token (``gen_sample_0``) becomes a field so it can serve
    as a join key; any other bare token is ignored.
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

    Coercion happens once, visibly, in ``_as_number``. ``n_overflow`` counts rows with
    MORE fields than the header: their values are shifted, so the caller must error.
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


# --- the join ---


class Join(NamedTuple):
    """The outcome of pairing FASTA entries with metrics rows.

    ``mode`` is "key" | "positional" | "failed"; ``keys`` the fields rows were matched
    on; ``checked`` the shared fields compared value by value (empty with "positional"
    means nothing corroborated the order); ``pairs`` (entry, row) in FASTA order;
    ``note`` why it failed.
    """

    mode: str
    keys: list[str]
    checked: list[str]
    pairs: list[tuple[tuple[str, str], dict[str, str]]]
    note: str = ""


def is_identity_field(name: str, bare: set[str]) -> bool:
    """Is this header field an IDENTIFIER of a draw, rather than a measurement?

    Only identifiers may key or cross-check a join: comparing a rounded ``seq_rec=0.43``
    with the CSV's full precision would fail a correct pairing.
    """
    return name.casefold() in IDENTITY_FIELDS or name in bare


def bare_field_names(header: str) -> set[str]:
    """Fields ``parse_header`` made from bare ``<word>_<int>`` tokens.

    Positional identifiers by construction (``gen_sample_0``).
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

    For a design this is its id (the metrics ``name``); for the echo, ``<stem>_native``.
    """
    return header.split(",")[0].strip()


# Stage-1 seeds under --seed-mode generate (design_ph.py:1772-1779): not designs.
_STAGE1_SEED_RE = r"stage1_seed_s\d+"


def is_preamble_header(header: str, stem: str | None) -> bool:
    """Is this header a non-design entry (the input echo, or a stage-1 seed)?

    Exact match on the whole first token, bare or as ``<stem>_<name>`` (what design_ph
    writes, e.g. ``<stem>_native``, ``<stem>_stage1_seed_s0``). Never a substring:
    ``>native_like_lam0_s0`` is a design.
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


# Metrics columns that may hold the design id (design_ph writes ``name``).
DESIGN_ID_COLUMNS = ("name",)

# Identifiers cross-checked after pairing (never measurements: header 3 dp, CSV 4 dp).
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
    """Pair on the design id: the header's first token and the metrics ``name`` cell.

    ``None`` when no id evidence exists, so the shared-field logic runs. A contradiction
    (different draws, a disagreeing identifier) is a failed join, never a positional
    guess. Only this join tells apart draws identical in sequence and ddE (_s0 vs _s1).
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

    # If ids repeat (multi-lambda sweep?), (id, lambda) is still a key, never position.
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
            # Count a field as checked only where both sides had a value.
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
    """Pair FASTA entries with metrics rows using shared identity fields as evidence."""
    raw = [h for h, _ in entries]
    headers = [parse_header(h) for h in raw]
    bare = set().union(*(bare_field_names(h) for h in raw)) if raw else set()

    # Strongest evidence first: the design id, which the program writes on both sides.
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
                # Unique keys with different value sets: the files do not correspond.
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
        # checked stays []: nothing was compared, so naming fields would overclaim.
        return Join(
            "failed", [], [], [],
            f"{len(entries)} FASTA entries vs {len(rows)} metrics rows, and no shared "
            f"join key in the headers",
        )

    pairs = list(zip(entries, rows))
    # A non-unique shared field (lambda) still corroborates a positional pairing.
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

    As ``Join``, with "missing" for an absent ``.ext.fa``; ``seqs`` are the extended
    sequences in the order of ``entries``.
    """

    mode: str
    checked: list[str]
    seqs: list[str]
    note: str = ""


def align_ext(
    entries: list[tuple[str, str]], ext_entries: list[tuple[str, str]]
) -> ExtJoin:
    """Extended sequences in the order of ``entries``, mirroring ``align_entries``.

    The caller records the outcome in ``ext_status`` / ``ext_checked_fields``.
    """
    if not ext_entries:
        return ExtJoin("missing", [], [])
    raw = [h for h, _ in entries]
    ext_raw = [h for h, _ in ext_entries]
    headers = [parse_header(h) for h in raw]
    ext_headers = [parse_header(h) for h in ext_raw]
    bare = set().union(*(bare_field_names(h) for h in raw + ext_raw))

    # Strongest evidence first, as in _design_id_join (.ext.fa repeats the headers).
    h_ids: list[tuple[str, ...]] = [(first_token(h),) for h in raw]
    e_ids: list[tuple[str, ...]] = [(first_token(h),) for h in ext_raw]
    key_names = ["design_id"]
    usable = (
        all(i[0] for i in h_ids)
        and all(i[0] for i in e_ids)
        and len(set(h_ids)) == len(h_ids)
        and len(set(e_ids)) == len(e_ids)
    )
    # Composite fallback as in _design_id_join; pairing by position would be a guess.
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
                # Both sides must carry a value for the field to count as checked.
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


# --- centers.json ---

_CENTER_LIST_KEYS = ("centers", "selected_centers", "chosen_centers", "pinned_centers")
_CHAIN_KEYS = ("chain", "chain_id", "chain_name", "binder_chain", "asym_id")
# Most explicit first; the matched key is reported in ``centers_resnum_key``.
_RESNUM_KEYS = ("resnum", "res_num", "res_id", "resi", "residue", "seqpos", "position", "pos")
# Keys read as 0-based indices: only these get ``_index0`` and the sequence check.
_INDEX0_KEYS = frozenset({"position", "pos", "seqpos", "index", "idx"})
# Keys that may carry the protonation TOKEN; only a full token (HIS-S, ASP) counts.
_STATE_KEYS = ("token", "state", "target_state", "microstate", "residue_state", "aa")
# One key, not aliases (a bare ``score`` could be anything); collected as the minimum.
_PLACEMENT_KEYS = ("placement_score",)

# Keys that record which positions could move; design_ph writes ``free_positions``.
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

    From the token, never ``--center-mode``: ``HIS-S`` -> ``H``, ``GLU-D`` -> ``E``; a
    bare 3-letter code works. A bare 1-letter value is NOT accepted: a ``state`` of "S"
    is a protonation state, not serine. ``None`` means "cannot check", i.e. NA.
    """
    head = str(token).split("-")[0].strip().upper()
    if len(head) != 3:
        return None
    return THREE_TO_ONE.get(head)


def _center_token(value: Any) -> str | None:
    """``value`` if it is a protonation token (``XXX`` or ``XXX-Y``), else ``None``."""
    if value is None:
        return None
    text = str(value).strip()
    head = text.split("-")[0]
    return text if len(head) == 3 and head.isalpha() else None


def _parse_metrics_centers(cell: Any) -> set[tuple[int, str]] | None:
    """``{(index, token)}`` from a metrics ``centers`` cell (``10:HIS-S``), or ``None``.

    The multi-centre separator is unverified, so ``;`` ``,`` ``|`` and whitespace are
    accepted; ``None`` becomes NA in the caller, never a failure.
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
    """``columns`` for every child row, plus the evidence kept out of them.

    ``checks`` is ``[(index0, letter, token)]`` for sequence-checking each centre;
    ``pairs`` is ``{(index0, token)}`` for cross-checking a metrics row's own cell.
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

    Schema: ``{"centers": [{"position": 10, "token": "HIS-S", "placement_score": 0.7}],
    "free_positions": [...], "center_mode": ...}``. ``position`` is a 0-based index into
    the binder chain; there is no chain, residue number or charge-field range
    (``charge_field_*`` stay NA). ``centers_source`` names the key the centres came
    from; ``unparsed`` leaves only the centre columns NA.
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
                    # Inexpressible centre: the set is incomparable, so NA.
                    pairs_complete = False
                if resnum_key:
                    # Label the numbering convention in the column.
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
        # Minimum, not mean: one weak centre dooms the placement. NA, never 0.
        out["placement_score_min"] = min(placements) if placements else pd.NA
    out["centers_source"] = source

    if isinstance(obj, dict):
        # placement_metric picks the centres and its default varies by commit.
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

    The one check that catches a wrong index base. Returns ``(True, "")``,
    ``(False, why)``, or ``(NA, "")`` when it could not run (no centres or tokens, or a
    multi-chain ``sequence`` under ``--output-scope all``).
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
    """Does this row's ``centers`` cell name the same centres as ``centers.json``?

    Returns ``(True, "")``, ``(False, why)``, or ``(NA, "")`` when either side cannot be
    read as ``(index, token)`` pairs: an unverified separator is never a failure.
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

    One token per residue, and every token that resolves to a residue must BE that
    residue; an unresolvable token is skipped (unrecognised is not disagreement).
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
    """``parts`` as floats, or ``None`` if any is not finite.

    ``nan`` is rejected: design_ph writes it per centre on the apo side with no target.
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
    """Per-centre energy columns for ONE design row.

    A non-empty ``note`` means the row contradicts itself: the caller fails the parent.
    """

    columns: dict[str, Any]
    note: str = ""


def per_centre_energies(
    row: dict[str, str], lower_to_col: dict[str, str], n_centers: Any
) -> PerCentre:
    """Split ``dE_complex`` / ``dE_apo`` into numeric, verbatim and derived columns.

    ``ddE_min`` is the point: ``ddE_mean`` cannot tell two half-switching histidines
    from one that carries all of it. Contradictions fail the parent (``note``): counts
    that disagree with each other, the row's ``centers`` cell or ``centers.json``, or a
    ``mean(per-centre ddE)`` more than 2e-3 from the file's ``ddE_mean``.
    Source: design_ph.py:1859-1868.
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
    # The row's own centres cell; centers.json is compared by centers_match_metrics.
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
        # nan apo (no target) or unparseable: keep the strings, derived columns stay NA.
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

    NA on a length mismatch, never a truncated count (not versions of one design).
    """
    if not a or not b or len(a) != len(b):
        return pd.NA
    return sum(1 for x, y in zip(a, b) if x != y)


def blank_target_columns(data: dict[str, Any]) -> None:
    """Under ``--binder-alone``, blank every target-dependent column in place.

    NA, never 0 or an error: the row is ``OK`` and a zero would read as a measurement.
    Only present columns are touched, and ``placement_metric`` is corrected to
    ``dE_complex``, which is what actually ranked the centres (no apo graph).
    """
    for key in TARGET_DEPENDENT_NUMERIC:
        if key in data:
            data[key] = pd.NA
    for key in TARGET_DEPENDENT_TEXT:
        if key in data:
            data[key] = ""
    if data.get("placement_metric"):
        data["placement_metric"] = "dE_complex"


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


# --- the collector ---


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


# Flags recoverable from the argv tail, by column (--lam has its own TSV field).
_TAIL_FLAGS = {"model_noise": "--checkpoint_model_noise",
               "charge_prior": "--charge_prior",
               "protect_titratable": "--protect_titratable",
               "apo_weight": "--apo_weight"}


def noise_and_lam_from_tasks(out_dir: Path) -> tuple[str, str]:
    """``(model_noise, lam_sweep)`` recovered from the sub-manifests, else ``("", "")``.

    Values that differ across rows are unknown, never a guessed default.
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
    """Run parameters recoverable from the sub-manifests, ``{column: value}``.

    For leaves whose ``.meta.json`` predates a key. A value not identical across every
    row is omitted, never defaulted: an absent key must read as NA.
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

    ``path=""`` so a later ``-i atomiumph_path`` cannot hand a .metrics.csv to a tool.
    """
    return Collected(name=f"{parent}_p0", parent=parent, path="",
                     status=f"error: {message}", data=data)


def collect_atomiumph(ctx: CollectCtx) -> CollectEach:
    """Per-parent collector: one child row ``<parent>_p<i>`` per designed sequence."""
    index = index_outputs(ctx.out_dir)
    meta = read_meta(ctx.out_dir)
    expected = meta.get(EXPECTED_KEY)
    expected = expected if isinstance(expected, int) else pd.NA
    # --binder-alone: no target in the graph; absent (older leaves) means False.
    binder_alone = bool(meta.get(BINDER_ALONE_KEY, False))
    run_columns = {
        "binder_alone": binder_alone,
        "binder_chains": str(meta.get(BINDER_CHAINS_KEY, "")),
        "target_chains": str(meta.get(TARGET_CHAINS_KEY, "")),
        # Run-wide, and both decide how a row must be read.
        "center_mode": str(meta.get(CENTER_MODE_KEY, "")),
        "output_scope": str(meta.get(OUTPUT_SCOPE_KEY, "")),
    }
    # Model noise decides which centres are picked; older leaves fall back to the tasks.
    sidecar = {
        "model_noise": str(meta.get(MODEL_NOISE_KEY, "") or ""),
        "lam_sweep": str(meta.get(LAM_SWEEP_KEY, "") or ""),
        "protect_titratable": str(meta.get(PROTECT_TITRATABLE_KEY, "") or ""),
        "apo_weight": str(meta.get(APO_WEIGHT_KEY, "") or ""),
        "charge_prior": str(meta.get(CHARGE_PRIOR_KEY, "") or ""),
    }
    if not all(sidecar.values()):
        from_tasks = run_params_from_tasks(ctx.out_dir)
        for key, value in sidecar.items():
            sidecar[key] = value or from_tasks.get(key, "")
    run_columns.update(sidecar)
    warned: set[str] = set()
    stale_warned: set[str] = set()

    # Advisory only: mtimes are unsound (partial reruns, copies); submit is the guard.
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

    # Batch mode does not abort on a per-structure failure: report the count once.
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
            # Reported above; the framework leaves the parent uncollected.
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
            # NA, not 0: an error row that never parsed the FASTA observed no preamble.
            "n_preamble_dropped": pd.NA,
            # NA, never False: an error row has not disagreed, it was just not checked.
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
            # Any file of the stem anchors the image_commit.txt walk, if there is one.
            anchor = next((p for paths in files.values() for p in paths), None)
            base["image_commit"] = commit_for(anchor, ctx.out_dir) if anchor else ""
            yield _error_row(
                d.name, f"incomplete output: {d.name}{missing} is absent",
                {**base, "join_mode": "failed"},
            )
            return

        base["image_commit"] = commit_for(metrics_path, ctx.out_dir)
        columns, rows, n_overflow = read_metrics(metrics_path)
        # The stem is the parent name, and design_ph writes the echo as <stem>_native.
        all_entries = parse_fasta(fasta_path)
        entries, dropped = drop_preamble(all_entries, d.name)
        # Keep the `_native` echo (never a stage-1 seed): the n_diff_native baseline.
        native = next(
            (
                seq
                for head, seq in all_entries[: max(dropped, 0)]
                if first_token(head).casefold()
                in (f"{d.name.casefold()}_native", "native")
            ),
            "",
        )
        # join=" ": .ext.fa is space-separated tokens, so wrapped lines must not fuse.
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
            # Never paper over this: an error row, not numbers on the wrong sequences.
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

        # Centre trust check: the sequence must carry the token's residue at its index.
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

        # Second check: centers.json vs each row's own centers cell; unreadable is NA.
        lower_cols = {c.casefold(): c for c in columns}

        # Per-centre energy split and consistency checks, before any row is yielded.
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

        # .ext.fa has its own status: its absence must not hide a good ddE_mean.
        ext_status = {
            "missing": f"error: {d.name}.ext.fa is absent",
            "failed": f"error: cannot map sequence_ext: {ext.note}",
        }.get(ext.mode, "OK")

        for i, ((header, sequence), row) in enumerate(join.pairs):
            data: dict[str, Any] = dict(base)
            data["entry_index"] = i
            data["sequence"] = sequence
            data["seq_len"] = len(sequence)
            # Distance from the echoed input (NA on a length mismatch or absent echo).
            data["n_diff_native"] = hamming(sequence, native)
            data["centers_match_metrics"] = match_checks[i][0]
            # sequence_ext must be THIS design, token for token.
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
            # Metrics columns pass through verbatim; new upstream columns need no code.
            for column in columns:
                key = column
                if key in RESERVED:
                    # Known collisions (name, centers) become metrics_<name> quietly.
                    key = f"metrics_{column}"
                    if column not in warned and column not in EXPECTED_COLLISIONS:
                        warned.add(column)
                        print(
                            f"atomiumph: metrics column {column!r} collides with a "
                            f"trust column; collected as {key!r}"
                        )
                if column.casefold() in {c.casefold() for c in PER_CENTRE_COLUMNS}:
                    # Numeric only: several centres -> NA (see <name>_centres).
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
            # Last: blank what is undefined without a target.
            if binder_alone:
                blank_target_columns(data)
            yield Collected(
                name=f"{d.name}_p{i}",
                parent=d.name,
                path=fasta_path,
                status="OK",
                data=data,
            )

    return one
