#!/usr/bin/env python3
"""
Collect af3conf results into the table.

Reads the per-design <name>.tsv files af3conf_worker.py wrote under
<run_dir>/<table>/af3conf[_<label>]/ and merges them back as <prefix>_status,
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

NA (empty) rather than 0 marks a value not measured: every 'error: ...' row
(binder_seq_match aside) and every 'missing' row has all numeric columns NA, and a
per-chain ipSAE score is NA for a pair with no residue pair under the cutoff. Measured
zeros are ``ipsae_nres_binder_<x>`` / ``ipsae_nres_<x>``. ``chains`` and
``input_source`` stay on error rows: they are what was ASKED.

Statuses: 'OK'; 'error: ...' (unresolvable path, missing/unreadable structure, a file
not named <prefix>_model.cif, an absent chain, a nucleic acid among the named chains, a
missing sibling JSON, a token/residue/atom mismatch, a chain count or chain_pair_pae_min
the summary disagrees with, a non-finite value, a crashed worker); 'missing' when the
design has no TSV.

Per-chain columns are dynamic (``iptm_<x>`` and eight ``ipsae_*_<x>`` per
--target-chains entry): names come from the per-design TSV header, and every per-chain
column the table already carries under this leaf is blanked on every row yielded. """

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

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
RESULT_COLUMNS = STR_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS

KNOWN = set(RESULT_COLUMNS) | {"name", "status", "path"}

# Legitimately NA on an OK row: no pair under the PAE cutoff, or no --sequence-column.
OPTIONAL_ON_OK = {
    "ipsae_binder_mean", "ipsae_binder_max", "ipsae_binder_top2", "binder_seq_match",
}

# Numeric values kept on an ERROR row when measured before the error.
KEPT_ON_ERROR = ["binder_seq_match"]


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


def collect_af3conf(ctx: CollectCtx) -> CollectEach:
    """Per-design collector; variants on one table are told apart by ``--dir-label``."""

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

        dynamic = list(dict.fromkeys([*_dynamic_columns(result_df.columns), *stale]))
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
