#!/usr/bin/env python3
"""
Collect af3conf results into the table.

Reads the per-design <name>.tsv files af3conf_worker.py wrote under
<run_dir>/<table>/af3conf[_<label>]/.submissions/<NN>/ (newest submission that lists a
design in its af3conf_tasks/ answers for it; a leaf from before submission dirs is read
as the oldest one) and merges them back as <prefix>_status,
<prefix>_path (the predicted structure that was scored) and one <prefix>_<field> column
per metric. Every name shared with boltzconf is IDENTICAL to boltzconf's, so
af3conf_ipsae_binder_max sits beside boltzconf_ipsae_binder_max. Names, units, NA rules
and how to read them: ``.claude/skills/af3conf/SKILL.md``.

  trust      n_tokens, n_res_struct, chains ('binder|targets', an echo), binder_n_res,
             n_target_res, plddt_bfac_maxdiff (per-atom pLDDT vs B-factor; above 0.05 is
             an error), chain_pae_min_maxdiff / chain_pae_min_margin (does the summary's
             POSITIONAL chain index match the letters; above 0.1 is an error, a margin
             <= 0.1 means the check could not have caught a swap), stem_matches_name,
             binder_seq_match (NA = UNCHECKED; the only numeric kept on an error row),
             ranking_score, input_source ('own table' or 'lineage').
  binder     binder_plddt, binder_plddt_min, chains_ptm_binder.
  interface  iptm_binder_mean / _max / _top2 and iptm_<x> (the symmetrised binder-target
             chain_pair_iptm), and boltzconf's ipSAE port on the RAW ASYMMETRIC PAE: per
             target chain x ``ipsae_*_<x>``, binder-level
             ``ipsae_binder_mean/_max/_top2`` and the key ``ipsae_pae_cutoff``.
  ipSAE min  DERIVED HERE, not written by the worker, and OUR aggregate (the reference's
             ipsae.py has only "asym" and "max" rows, no min); higher is better, 0-1:
             ``ipsae_min_<x>`` = min(ipsae_binder_to_<x>, ipsae_<x>_to_binder), the
             weaker direction of that one pair, NA when EITHER direction is NA; and
             ``ipsae_binder_min`` = the MAX over target chains of ipsae_min_<x>, i.e.
             the best chain's weaker direction (NOT the weakest chain), NA when no chain
             has both directions. Both are a pure function of two columns every af3conf
             TSV carries, so a TSV written before them collects them too
             (``ipsae_min_columns``, source-identical in collect_boltzconf.py).

NA (empty) rather than 0 marks a value not measured: every 'error: ...' row
(binder_seq_match aside) and every 'missing' row has all numeric columns NA, and a
per-chain ipSAE score is NA for a pair with no residue pair under the cutoff. Measured
zeros are ``ipsae_nres_binder_<x>`` / ``ipsae_nres_<x>``. ``chains`` and
``input_source`` stay on error rows: they are what was ASKED.

Statuses: 'OK'; 'error: ...' (unresolvable path, missing/unreadable structure, a file
not named <prefix>_model.cif, an absent chain, a nucleic acid among the named chains, a
missing sibling JSON, a token/residue/atom mismatch, a chain count or chain_pair_pae_min
the summary disagrees with, a non-finite value, a crashed worker); 'missing' when the
newest submission listing the design has no TSV for it, or no submission has one.

Per-chain columns are dynamic (``iptm_<x>`` and nine ``ipsae_*_<x>`` per --target-chains
entry: eight from the TSV, ``ipsae_min_<x>`` derived here): names come from the
per-design TSV header, and every per-chain column the table already carries under this
leaf is blanked on every row yielded. """

import math
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx
from prosapia.core.naming import submission_dirs

TASKS_DIRNAME = "af3conf_tasks"

# Bare column names, leaf-prefixed by the driver; kept in step with the worker by name.
FLOAT_COLUMNS = [
    "plddt_bfac_maxdiff",
    "chain_pae_min_maxdiff",
    "chain_pae_min_margin",
    "ranking_score",
    "binder_plddt",
    "binder_plddt_min",
    "chains_ptm_binder",
    "iptm_binder_mean",
    "iptm_binder_max",
    "iptm_binder_top2",
    "ipsae_binder_mean",
    "ipsae_binder_max",
    "ipsae_binder_top2",
    "ipsae_pae_cutoff",
]
INT_COLUMNS = [
    "n_tokens",
    "n_res_struct",
    "binder_n_res",
    "n_target_res",
    "stem_matches_name",
    "binder_seq_match",
]
STR_COLUMNS = ["chains", "input_source"]
# Fixed names computed HERE from the TSV (never read from it, so never required on OK).
DERIVED_COLUMNS = ["ipsae_binder_min"]
RESULT_COLUMNS = STR_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS + DERIVED_COLUMNS

# A fixed ipsae_*/iptm_* name missing here would be parsed as a per-chain column.
KNOWN = set(RESULT_COLUMNS) | {"name", "status", "path"}

# The forward-direction per-chain prefix; its suffix is the chain key.
IPSAE_FORWARD_PREFIX = "ipsae_binder_to_"

# Legitimately NA on an OK row: no pair under the PAE cutoff, or no --sequence-column.
OPTIONAL_ON_OK = {
    "ipsae_binder_mean", "ipsae_binder_max", "ipsae_binder_top2", "binder_seq_match",
}

# Numeric values kept on an ERROR row when measured before the error.
KEPT_ON_ERROR = ["binder_seq_match"]


def ipsae_min_columns(
    directions: dict[str, tuple[float | None, float | None]],
) -> dict[str, float | None]:
    """ipsae_min_<k> per target chain and ipsae_binder_min: OUR aggregate, not ipsae.py's.

    ``directions`` maps a chain key k to (ipsae_binder_to_<k>, ipsae_<k>_to_binder).
    A non-finite direction (NaN/inf) counts as None, so neither min nor max depends on
    argument order. ipsae_min_<k> is the WEAKER direction of that one binder-vs-k
    pair, None when either direction is None (a direction with no sub-cutoff pair has
    no ipSAE, and a 0.0 would read as measured). ipsae_binder_min is the MAX over
    chains of ipsae_min_<k>: the best chain on which BOTH directions are DEFINED (each
    has >= 1 sub-cutoff PAE pair; the weaker may still be low), NOT the weakest chain
    (a min over chains would be None or ~0 for any one-protomer binder); None when no
    chain has both. Higher is better for both, 0-1. The reference writes only "asym"
    and "max" (= the larger direction) rows; it has no min. Source-identical in
    collect_boltzconf.py and collect_af3conf.py.
    """

    def finite(value: float | None) -> float | None:
        return value if value is not None and math.isfinite(value) else None

    out: dict[str, float | None] = {}
    for key, (forward, reverse) in directions.items():
        forward, reverse = finite(forward), finite(reverse)
        out[f"ipsae_min_{key}"] = (
            None if forward is None or reverse is None else min(forward, reverse)
        )
    defined = [v for v in out.values() if v is not None]
    out["ipsae_binder_min"] = max(defined) if defined else None
    return out


def _ipsae_chain_keys(header: Iterable[str]) -> list[str]:
    """Chain keys of the TSV's ipSAE block, read off ``ipsae_binder_to_<k>``."""
    return [c[len(IPSAE_FORWARD_PREFIX) :] for c in header if c.startswith(IPSAE_FORWARD_PREFIX)]


def _is_per_chain(column: str) -> bool:
    """A per-target-chain column: an iptm_*/ipsae_* name that is not in KNOWN."""
    if column in KNOWN:
        return False
    return column.startswith("iptm_") or column.startswith("ipsae_")


def _dynamic_columns(header: Iterable[str]) -> list[str]:
    """The per-target-chain columns present in this design's TSV."""
    return [c for c in header if _is_per_chain(c)]


def _stale_per_chain(ctx: CollectCtx) -> list[str]:
    """Per-chain columns the TABLE already carries under this leaf.

    Without them a 'missing' or 'error: empty tsv' row (no header to read), or a rerun
    with different --target-chains, would keep an earlier run's values; every yielded
    row blanks them. """
    leaf = ctx.leaf
    out: list[str] = []
    for column in ctx.df.columns:
        text = str(column)
        if text.startswith(f"{leaf}_"):
            bare = text[len(leaf) + 1 :]
            if _is_per_chain(bare):
                out.append(bare)
    return out


def _empty(dynamic: Iterable[str]) -> dict[str, Any]:
    return {c: None for c in [*RESULT_COLUMNS, *dict.fromkeys(dynamic)]}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _submitted(sub: Path) -> set[str]:
    """The design names a submission's task files list (first field of each line)."""
    tasks_dir = sub / TASKS_DIRNAME
    if not tasks_dir.is_dir():
        if any(sub.glob("*.tsv")):
            raise RuntimeError(
                f"submission {sub} holds per-design TSVs but no {TASKS_DIRNAME}/: the "
                f"designs it ran are unknown, so a killed one would collect an older "
                f"submission's result. Restore {TASKS_DIRNAME}/ or move the submission aside."
            )
        return set()
    names: set[str] = set()
    for task_file in tasks_dir.glob("task_*.tsv"):
        for line in task_file.read_text().splitlines():
            if line.strip():
                names.add(line.split("\t", 1)[0])
    return names


def collect_af3conf(ctx: CollectCtx) -> CollectEach:
    """Per-design collector; variants on one table are told apart by ``--dir-label``."""

    stale = _stale_per_chain(ctx)
    oldest = ctx.out_dir == submission_dirs(ctx.leaf_dir)[-1]
    flat = ctx.out_dir == ctx.leaf_dir
    submitted = set() if flat else _submitted(ctx.out_dir)
    # The driver stops reading flat TSVs once .submissions exists; they predate every submission.
    legacy = ctx.leaf_dir if oldest and not flat else None

    def _tsv(name: str) -> Path | None:
        for base in (ctx.out_dir, None if name in submitted else legacy):
            if base is not None and (path := base / f"{name}.tsv").is_file():
                return path
        return None

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = _tsv(d.name)
        if tsv_path is None:
            # Silence hands the design to an older submission; a listed one answers here.
            if oldest or d.name in submitted:
                yield Collected(status="missing", path="", data=_empty(stale))
            return

        result_df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty(stale))
            return

        ipsae_keys = _ipsae_chain_keys(result_df.columns)
        dynamic = list(
            dict.fromkeys(
                [
                    *_dynamic_columns(result_df.columns),
                    *(f"ipsae_min_{k}" for k in ipsae_keys),  # derived below
                    *stale,
                ]
            )
        )
        row = result_df.iloc[0]
        status = str(row["status"])
        data = _empty(dynamic)
        for col in STR_COLUMNS:
            if (value := _present(row, col)) is not None:
                data[col] = str(value)
        if status != "OK":
            # Numeric columns stay NA, except binder_seq_match when the worker wrote it.
            for col in KEPT_ON_ERROR:
                if (value := _present(row, col)) is not None:
                    data[col] = int(float(value))
            yield Collected(status=status, path="", data=data)
            return

        # An OK row must carry every always-measurable column (iptm_<x>, ipsae_nres_*).
        in_header = set(result_df.columns)
        optional = {
            *OPTIONAL_ON_OK,
            *(
                c
                for c in dynamic
                if (c.startswith("ipsae_") and not c.startswith("ipsae_nres_"))
                or c not in in_header
            ),
        }
        missing = [
            c
            for c in [*INT_COLUMNS, *FLOAT_COLUMNS, *dynamic]
            if c not in optional and _present(row, c) is None
        ]
        # ipsae_min_<k> needs both directions IN THE HEADER (each may be NA as a value).
        missing += [
            f"ipsae_{k}_to_binder"
            for k in ipsae_keys
            if f"ipsae_{k}_to_binder" not in in_header
        ]
        if missing:
            # A hole in an OK row is a contract breach, not a value.
            yield Collected(
                status=f"error: OK row missing {','.join(missing)}",
                path="",
                data=data,
            )
            return
        for col in INT_COLUMNS:
            if (value := _present(row, col)) is not None:
                data[col] = int(float(value))
        for col in FLOAT_COLUMNS + dynamic:
            if (value := _present(row, col)) is not None:
                data[col] = float(value)
        # From the two parsed direction columns, so an old TSV collects them too.
        data.update(
            ipsae_min_columns(
                {
                    k: (data[f"ipsae_binder_to_{k}"], data[f"ipsae_{k}_to_binder"])
                    for k in ipsae_keys
                }
            )
        )
        path = _present(row, "path")
        yield Collected(status=status, path="" if path is None else str(path), data=data)

    return one
