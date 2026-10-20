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

NA (empty) rather than 0 marks a value not measured: every 'error: ...' and 'missing'
row has all numeric columns NA, and ipae_* are NA when there are no pairs to average;
``n_if_pairs`` tells the two apart (0 with status OK is a measured zero). ``chains`` is
kept on error rows because it is what was ASKED.

Statuses: 'OK'; 'error: ...' (unresolvable path, missing/unreadable structure, an absent
chain, a missing sibling JSON/npz, a token count or JSON chain count that disagrees with
the structure, pLDDT that disagrees with the B-factors, a crashed worker); 'missing'
when the design has no TSV.

Per-chain columns are dynamic (``iptm_<x>``, ``ipae_<x>_mean`` and eight ``ipsae_*_<x>``
per --target-chains entry): names come from the per-design TSV header, and a 'missing'
design reads back NA.
"""

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
RESULT_COLUMNS = STR_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS

# Per-chain columns come from the TSV header and from the table (see _stale_per_chain).
KNOWN = set(RESULT_COLUMNS) | {"name", "status", "path"}


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

        dynamic = [*_dynamic_columns(result_df.columns), *stale]
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
        path = _present(row, "path")
        yield Collected(status=status, path="" if path is None else str(path), data=data)

    return one
