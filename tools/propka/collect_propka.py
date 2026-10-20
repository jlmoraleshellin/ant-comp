#!/usr/bin/env python3
"""
Collect propka results into the table.

Reads the per-design <name>.tsv files propka_worker.py wrote under
<run_dir>/<table>/propka[_<label>]/ and merges them back as <prefix>_status,
<prefix>_path (the structure PROPKA actually read: a converted PDB when the input was a
CIF) and one <prefix>_<field> column per metric. Every column and how to read it:
``.claude/skills/propka/SKILL.md``.

Spelling is BRITISH throughout (centre, centres): propka_center1_pka_complex is a
KeyError.

  trust    propka_version, binder_chains and binder_chains_n_res / _n_het_res /
           _sequence, expect_binder_len, chain_ids_found, centres_source, free_only, centres_n /
           centres_reported / centres_his_ok / centres_pka_ok, the stripped-hydrogen and
           warning counts (gate on n_propka_warnings_FREE; the complex count is the
           target's), the binder's terminal groups per run and term_groups_match, the pH
           echo (ph_values, reference_ph, release_factors, f_model) and n_binder_his_no_pka.
  science  per centre k (first MAX_CENTRE_COLUMNS): centre<k>_chain / _resnum /
           _pka_complex / _pka_free / _dpka / _buried_* / _n_coupled_*; per_centre_pka_detail,
           min_dpka, n_centres_coupled, n_binder_his, binder_his_pka; the pH layer
           (f_, f_ceiling_, switch_ <pH tag> <site set>, ph_half_, ph_release_,
           f_dips_below_) and the free-state fold stability (stab_*, with stab_error).

Under ``free_only`` everything that needs the complex run is NA by construction.
Otherwise NA (empty) rather than 0 marks a value not measured: every 'error: ...' row
and every 'missing' row has all science columns NA and keeps only the echo and the trust
facts established before the failure. 0 is a measurement.

Statuses: 'OK'; 'error: ...' (unresolvable structure or centres, a missing or
unreadable file, an absent binder chain, no partner chain without --free-only, a
residue number that is not in the chain, a centre that is not a histidine or whose
number repeats, a binder length disagreeing with --expect-binder-len-col, a PROPKA
failure, a centre PROPKA reported no pKa for, a crashed worker); 'missing' when the
design has no TSV.
"""

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# The pH-dependence column names depend on the run's --ph-values / --release-factors, so
# they are rebuilt per row from that row's own echo, with the SAME functions the worker
# named them with. Nothing here is a hand-kept copy; the worker must keep importing only
# the standard library (it also runs in the tool image) because a module-level failure
# there breaks tool DISCOVERY for every sapia command.
from .propka_worker import (
    CENTRE_COMPLEX_FIELDS,
    CENTRE_FIELDS,
    CENTRE_INT_FIELDS,
    CENTRE_STR_FIELDS,
    DEFAULT_PH_VALUES,
    DEFAULT_REFERENCE_PH,
    DEFAULT_RELEASE_FACTORS,
    MAX_CENTRE_COLUMNS,
    derived_columns,
    parse_ph_values,
    parse_reference_ph,
    parse_release_factors,
)

# The pH columns for the DEFAULT flags; a run with other flags writes other names.
DEFAULT_DERIVED = derived_columns(
    list(DEFAULT_PH_VALUES), list(DEFAULT_RELEASE_FACTORS), DEFAULT_REFERENCE_PH
)
_DEFAULT_DERIVED_SET = set(DEFAULT_DERIVED)
# Among the pH columns ph_release_saturated is the only STRING, the f_dips_below_* flags
# the only INTEGERS, and the rest floats. Matched by name: a non-default
# --release-factors writes names no static list can hold.
DERIVED_STR = {"ph_release_saturated"}
DERIVED_INT_PREFIX = "f_dips_below_"

_CENTRE_STR = [
    f"centre{k}_{f}"
    for k in range(1, MAX_CENTRE_COLUMNS + 1)
    for f in CENTRE_FIELDS
    if f in CENTRE_STR_FIELDS
]
_CENTRE_INT = [
    f"centre{k}_{f}"
    for k in range(1, MAX_CENTRE_COLUMNS + 1)
    for f in CENTRE_FIELDS
    if f in CENTRE_INT_FIELDS
]
_CENTRE_FLOAT = [
    f"centre{k}_{f}"
    for k in range(1, MAX_CENTRE_COLUMNS + 1)
    for f in CENTRE_FIELDS
    if f not in CENTRE_STR_FIELDS | CENTRE_INT_FIELDS
]

STR_COLUMNS = [
    "propka_version",
    "binder_chains",
    "binder_chains_sequence",
    "chain_ids_found",
    "centres_source",
    "free_state_model",
    "ph_values",
    "reference_ph",
    "release_factors",
    "f_model",
    "ph_release_saturated",
    "binder_term_groups_complex",
    "binder_term_groups_free",
    "per_centre_pka_detail",
    "binder_his_pka",
    # Empty when the fold-stability pass succeeded; else why every other stab_* is NA.
    "stab_error",
    *_CENTRE_STR,
]
INT_COLUMNS = [
    "binder_chains_n_res",
    "binder_chains_n_het_res",
    "expect_binder_len",
    "free_only",
    "centres_n",
    "centres_reported",
    "centres_his_ok",
    "centres_pka_ok",
    "n_h_stripped",
    "n_h_stripped_guessed",
    "n_his_renamed",
    "n_propka_warnings_complex",
    "n_propka_warnings_free",
    "n_binder_his",
    "n_binder_his_no_pka",
    "n_centres_coupled",
    "term_groups_match",
    # stab_n_his_groups MUST equal centres_n: it catches a centre that matched no group.
    "stab_n_groups_titrated",
    "stab_n_his_groups",
    *[c for c in DEFAULT_DERIVED if c.startswith(DERIVED_INT_PREFIX)],
    *_CENTRE_INT,
]
FLOAT_COLUMNS = [
    "min_dpka",
    "stab_profile_check_maxdiff",
    "stab_share_check",
    "stab_opt_ph",
    "stab_opt_dg",
    *_CENTRE_FLOAT,
    *[
        c
        for c in DEFAULT_DERIVED
        if c not in DERIVED_STR and not c.startswith(DERIVED_INT_PREFIX)
    ],
]
RESULT_COLUMNS = STR_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS

# The MEASUREMENT, as opposed to the echo and the trust facts: every one is emptied on
# an error row.
SCIENCE_COLUMNS = [
    *[f"centre{k}_{f}" for k in range(1, MAX_CENTRE_COLUMNS + 1) for f in CENTRE_FIELDS],
    # A per-centre result of the centre loop, so an error inside it leaves it half-counted.
    "n_centres_coupled",
    "per_centre_pka_detail",
    "min_dpka",
    "n_binder_his",
    "binder_his_pka",
    "stab_n_groups_titrated",
    "stab_n_his_groups",
    "stab_profile_check_maxdiff",
    "stab_share_check",
    "stab_opt_ph",
    "stab_opt_dg",
    "stab_error",
    *DEFAULT_DERIVED,
]

# Columns that need the complex run: NA by construction on a free_only row.
_COMPLEX_ONLY = {
    "n_propka_warnings_complex",
    "binder_term_groups_complex",
    "term_groups_match",
    "n_centres_coupled",
    "per_centre_pka_detail",
    "min_dpka",
    "ph_release_saturated",
}

# What an OK row MAY leave NA. The per-centre columns are checked against this row's own
# centres_reported instead (a row with fewer centres leaves the rest NA, never 0); the
# stability block is optional on purpose, so a propka API change cannot turn sound pKa
# values into an error row.
OPTIONAL_ON_OK = {
    "expect_binder_len",
    "binder_term_groups_complex",
    "binder_term_groups_free",
    "stab_n_groups_titrated",
    "stab_n_his_groups",
    "stab_profile_check_maxdiff",
    "stab_share_check",
    "stab_opt_ph",
    "stab_opt_dg",
    "stab_error",
    *_CENTRE_STR,
    *_CENTRE_INT,
    *_CENTRE_FLOAT,
}


def _needs_complex(col: str) -> bool:
    """Is this column NA by construction when only the free state was run?"""
    if col in _COMPLEX_ONLY:
        return True
    if col.startswith("centre"):
        return col.split("_", 1)[1] in CENTRE_COMPLEX_FIELDS
    return col.startswith(
        ("f_ph", "switch_ph", "ph_release_", "f_dips_below_", "ph_half_complex_")
    )


def _derived_optional(col: str) -> bool:
    """A pH column that may be NA on an OK row (name-based: names follow the flags).

    ph_release_*, f_dips_below_*, ph_half_* and stab_* have their own legitimate NA
    cases, and every ``*_all`` column is NA as a block when a binder His has no pKa
    (checked conditionally in one()). Only f_*_designed columns are always required.
    """
    if col == "ph_release_saturated" or col.startswith(
        ("ph_release_", "ph_half_", "f_dips_below_", "stab_")
    ):
        return True
    return col.startswith(("f_ph", "f_ceiling_ph", "switch_ph")) and col.endswith("_all")


def _empty() -> dict[str, Any]:
    """Every DEFAULT-named column as NA, for rows with no TSV to read."""
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _report_stab_errors(out_dir: Path) -> None:
    """Print how many OK rows of this leaf carry a stab_error, ALWAYS (including 0).

    The stability block is non-fatal, so a systematic failure would otherwise read as
    a quietly empty stab_* column, i.e. as 'no design loses stability'.
    """
    rows = 0
    with_error = 0
    for tsv_path in sorted(out_dir.glob("*.tsv")):
        try:
            frame = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
            continue
        if frame.empty or "status" not in frame.columns:
            continue
        row = frame.iloc[0]
        if str(row["status"]) != "OK":
            continue
        rows += 1
        if str(row.get("stab_error", "")).strip():
            with_error += 1
    if rows:
        print(f"propka: stab_error: {with_error} of {rows} OK rows carry a stab_error")
        if with_error:
            print(
                "  (their stab_* columns are all NA; the pKa columns are unaffected. "
                "An NA stab_loss_* is NOT 'no stability loss')"
            )


def collect_propka(ctx: CollectCtx) -> CollectEach:
    """Per-design collector; variants on one table are told apart by ``--dir-label``."""
    _report_stab_errors(ctx.out_dir)

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row = result_df.iloc[0]
        status = str(row["status"])

        # The pH column names come from this row's own echo. A row with no readable
        # echo (the worker-crash fallback TSV) can only be an error row, which is
        # blanked under the default names.
        try:
            row_ph_values = parse_ph_values(str(_present(row, "ph_values") or ""))
            derived = derived_columns(
                row_ph_values,
                parse_release_factors(str(_present(row, "release_factors") or "")),
                parse_reference_ph(str(_present(row, "reference_ph") or ""), row_ph_values),
            )
        except ValueError:
            if status == "OK":
                yield Collected(
                    status="error: OK row without a readable ph_values / reference_ph "
                    "/ release_factors echo",
                    path="",
                    data=_empty(),
                )
                return
            derived = list(DEFAULT_DERIVED)
        # The DEFAULT-named pH columns are dropped before this row's own are added, as
        # the worker's metric_columns_for does: keeping them would fail a legitimate
        # --ph-values 6.5 run for a missing f_ph7p4_designed it never wrote.
        columns = [c for c in RESULT_COLUMNS if c not in _DEFAULT_DERIVED_SET] + derived
        science = [c for c in SCIENCE_COLUMNS if c not in _DEFAULT_DERIVED_SET] + derived
        data: dict[str, Any] = {c: None for c in columns}

        # The echo and the trust facts survive an ERROR row (they diagnose a wrong chain
        # letter or a wrong binder length); no science column is filled on one.
        for col in columns:
            if (value := _present(row, col)) is not None:
                try:
                    if col in INT_COLUMNS or col.startswith(DERIVED_INT_PREFIX):
                        data[col] = int(value)
                    elif col in STR_COLUMNS or col in DERIVED_STR:
                        data[col] = str(value)
                    else:
                        data[col] = float(value)
                except ValueError:
                    yield Collected(
                        status=f"error: unreadable {col} cell {value!r}",
                        path="",
                        data={c: None for c in columns},
                    )
                    return

        if status != "OK":
            # Every science column stays NA: a 0 here would read as a measurement.
            for col in science:
                data[col] = None
            yield Collected(status=status, path="", data=data)
            return

        def breach(message: str) -> Collected:
            return Collected(
                status=f"error: {message}", path="", data={c: None for c in columns}
            )

        free_only = data.get("free_only") == 1
        missing = [
            c
            for c in columns
            if c not in OPTIONAL_ON_OK
            and not _derived_optional(c)
            and not (free_only and _needs_complex(c))
            and data.get(c) is None
        ]
        if missing:
            # A hole in an OK row is a contract breach, not a value.
            yield breach(f"OK row missing {','.join(missing)}")
            return
        for name in ("centres_his_ok", "centres_pka_ok"):
            if data[name] != data["centres_n"]:
                # Cannot come from this worker (it raises): a foreign TSV, never silent.
                yield breach(
                    f"OK row claims {data[name]} of {data['centres_n']} centres in {name}"
                )
                return
        n_reported = int(data["centres_reported"])
        if n_reported != min(int(data["centres_n"]), MAX_CENTRE_COLUMNS):
            yield breach(
                f"OK row reports {n_reported} centres of {data['centres_n']}, expected "
                f"{min(int(data['centres_n']), MAX_CENTRE_COLUMNS)}"
            )
            return
        # A hole inside the reported centres is a breach; past them it is the normal shape.
        absent = [
            f"centre{k}_{field}"
            for k in range(1, n_reported + 1)
            for field in CENTRE_FIELDS
            if field not in ("buried_complex", "buried_free")
            and not (free_only and field in CENTRE_COMPLEX_FIELDS)
            and data.get(f"centre{k}_{field}") is None
        ]
        # The 'all' site set is NA as a WHOLE when a binder His lacked a pKa, and that is
        # the only legitimate reason: otherwise a worker that stopped writing those
        # columns would collect as a clean OK row. f_ceiling_* is all a free_only row has.
        if data.get("n_binder_his_no_pka") == 0:
            prefixes = ("f_ceiling_ph",) if free_only else (
                "f_ph", "f_ceiling_ph", "switch_ph"
            )
            absent += [
                c
                for c in derived
                if c.startswith(prefixes) and c.endswith("_all") and data.get(c) is None
            ]
        # A dips flag is NA exactly when its release pH is NA; where the release pH is
        # a number the flag says whether "X-fold weaker at and below it" is true.
        absent += [
            f"f_dips_below_{c[len('ph_release_'):]}"
            for c in derived
            if c.startswith("ph_release_")
            and c != "ph_release_saturated"
            and data.get(c) is not None
            and data.get(f"f_dips_below_{c[len('ph_release_'):]}") is None
        ]
        if absent:
            yield breach(f"OK row missing {','.join(absent)}")
            return

        path = _present(row, "path")
        yield Collected(status=status, path="" if path is None else str(path), data=data)

    return one
