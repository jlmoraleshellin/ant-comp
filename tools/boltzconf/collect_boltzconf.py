#!/usr/bin/env python3
"""
Collect boltzconf results into the table.

Reads the per-design <name>.tsv files boltzconf_worker.py wrote under
<run_dir>/<table>/boltzconf[_<label>]/ and merges them back as <prefix>_status,
<prefix>_path (the predicted structure that was scored) and one <prefix>_<field>
column per metric. Names, units, NA rules and how to read them:
``.claude/skills/boltzconf/SKILL.md``.

  trust      n_tokens and n_res_struct (equal on an OK row by construction), chains
             ('binder|targets', an echo), binder_n_res, n_target_res, plddt_scale,
             plddt_bfac_maxdiff (token -> residue check: ~1e-5 when right, above 0.05
             an error, NA without B-factors; it does NOT check the JSON's chain-index
             -> letter map, only the chain COUNT), stem_matches_name, input_source
             ('own table' or 'lineage', the latter only with --allow-lineage).
  binder     binder_plddt, binder_plddt_min, chains_ptm_binder, pae_binder_self.
  interface  iptm_binder_mean / _max / _top2 and iptm_<x> (the symmetrised binder-target
             pair_chains_iptm); ipae_mean / ipae_median / ipae_<x>_mean (symmetrised
             PAE over interface residue PAIRS within contact_cutoff, NA when there are
             none); n_if_pairs (a real 0 means the binder touches nothing),
             n_if_res_binder, contact_cutoff (an ipae only compares at the same cutoff).
  ipSAE      Dunbrack's interface-restricted pTM, ported from the reference (URL, md5
             and fetch date in boltzconf_worker.py's docstring) and computed on the RAW
             ASYMMETRIC PAE, unlike the symmetrised columns above; no distance cutoff.
             ``ipsae_binder_to_<x>`` and ``ipsae_<x>_to_binder`` (both directions kept),
             ``ipsae_max_<x>``, the ``n0res`` support behind each score and
             ``d0res_max``, ``ipsae_nres_binder_<x>`` / ``ipsae_nres_<x>`` (breadth of
             the binder->x direction only; a measured 0 is exactly when
             ``ipsae_binder_to_<x>`` is NA), ``ipsae_binder_mean/_max/_top2`` (a no-pair
             chain contributes 0.0, as in the reference; NA only if none has a pair) and
             the comparability key ``ipsae_pae_cutoff``.
  ipSAE min  DERIVED HERE, not written by the worker, and OUR aggregate (the reference's
             ipsae.py has only "asym" and "max" rows, no min); higher is better, 0-1:
             ``ipsae_min_<x>`` = min(ipsae_binder_to_<x>, ipsae_<x>_to_binder), the
             weaker direction of that one pair, NA when EITHER direction is NA; and
             ``ipsae_binder_min`` = the MAX over target chains of ipsae_min_<x>, i.e.
             the best chain's weaker direction (NOT the weakest chain), NA when no chain
             has both directions. Both are a pure function of two columns every TSV
             since the ipSAE edit carries, so a TSV written before them collects them
             too (``ipsae_min_columns``, source-identical in collect_af3conf.py).

NA (empty) rather than 0 marks a value not measured: every 'error: ...' and 'missing'
row has all numeric columns NA, and ipae_* are NA when there are no pairs to average;
``n_if_pairs`` tells the two apart (0 with status OK is a measured zero). ``chains`` is
kept on error rows because it is what was ASKED.

Statuses: 'OK'; 'error: ...' (unresolvable path, missing/unreadable structure, an absent
chain, a missing sibling JSON/npz, a token count or JSON chain count that disagrees with
the structure, pLDDT that disagrees with the B-factors, a crashed worker); 'missing'
when the design has no TSV.

Per-chain columns are dynamic (``iptm_<x>``, ``ipae_<x>_mean`` and nine ``ipsae_*_<x>``
per --target-chains entry: eight from the TSV, ``ipsae_min_<x>`` derived here): names
come from the per-design TSV header, and a 'missing' design reads back NA.
"""

import math
from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names, leaf-prefixed by the driver; kept in step with the worker by name.
FLOAT_COLUMNS = [
    "plddt_bfac_maxdiff",
    "binder_plddt",
    "binder_plddt_min",
    "chains_ptm_binder",
    "iptm_binder_mean",
    "iptm_binder_max",
    "iptm_binder_top2",
    "ipae_mean",
    "ipae_median",
    "pae_binder_self",
    "contact_cutoff",
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
    "n_if_pairs",
    "n_if_res_binder",
    "stem_matches_name",
]
STR_COLUMNS = ["chains", "plddt_scale", "input_source"]
# Fixed names computed HERE from the TSV (never read from it, so never required on OK).
DERIVED_COLUMNS = ["ipsae_binder_min"]
RESULT_COLUMNS = STR_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS + DERIVED_COLUMNS

# Per-chain columns come from the TSV header and from the table (see _stale_per_chain).
# A fixed ipsae_* name missing here would be parsed as a per-chain column.
KNOWN = set(RESULT_COLUMNS) | {"name", "status", "path"}

# The forward-direction per-chain prefix; its suffix is the chain key.
IPSAE_FORWARD_PREFIX = "ipsae_binder_to_"


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
    """A per-target-chain column: anything in the three dynamic families not in KNOWN.

    Every fixed ipsae_* name is in KNOWN, so any other is a per-chain template.
    """
    if column in KNOWN:
        return False
    if column.startswith("iptm_"):
        return True
    if column.startswith("ipae_") and column.endswith("_mean"):
        return True
    return column.startswith("ipsae_")


def _dynamic_columns(header: Iterable[str]) -> list[str]:
    """The per-target-chain columns present in this design's TSV."""
    return [c for c in header if _is_per_chain(c)]


def _stale_per_chain(ctx: CollectCtx) -> list[str]:
    """Per-chain columns the TABLE already carries under this leaf.

    Without them a 'missing' or 'error: empty tsv' row (no header to read), or a rerun
    with different --target-chains, would keep an earlier run's values; every yielded
    row blanks them. """
    leaf = ctx.out_dir.name
    out: list[str] = []
    for column in getattr(ctx, "df", pd.DataFrame()).columns:
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


def collect_boltzconf(ctx: CollectCtx) -> CollectEach:
    """Per-design collector; variants on one table are told apart by ``--dir-label``."""

    # Scanned once per collect: the per-chain columns already on the table.
    stale = _stale_per_chain(ctx)

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty(stale))
            return

        result_df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty(stale))
            return

        ipsae_keys = _ipsae_chain_keys(result_df.columns)
        dynamic = [
            *_dynamic_columns(result_df.columns),
            *(f"ipsae_min_{k}" for k in ipsae_keys),  # derived below, NA on non-OK rows
            *stale,
        ]
        row = result_df.iloc[0]
        status = str(row["status"])
        data = _empty(dynamic)
        for col in STR_COLUMNS:
            if (value := _present(row, col)) is not None:
                data[col] = str(value)
        if status != "OK":
            # Numeric columns stay NA: a 0.0 would read as a measurement.
            yield Collected(status=status, path="", data=data)
            return

        # An OK row must carry every always-measurable column (iptm_<x>, ipsae_nres_*).
        in_header = set(result_df.columns)
        optional = {
            "ipae_mean",
            "ipae_median",
            "plddt_bfac_maxdiff",
            "pae_binder_self",
            # NA when no target chain has a residue pair under the PAE cutoff.
            "ipsae_binder_mean",
            "ipsae_binder_max",
            "ipsae_binder_top2",
            *(
                c
                for c in dynamic
                if c.startswith("ipae_")
                # NA when no sub-cutoff pair; ipsae_nres_* never optional (0 is real).
                or (c.startswith("ipsae_") and not c.startswith("ipsae_nres_"))
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
