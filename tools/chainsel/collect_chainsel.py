#!/usr/bin/env python3
"""
Collect chainsel results into the table.

Scans <run_dir>/<table>/chainsel/ for the per-design TSV files written by
chainsel.sh, and merges the extraction back into the table as <prefix>_status /
<prefix>_path (the extracted structure -- the whole point of the tool, and what a
downstream comparison consumes) plus one <prefix>_<field> column per field:

    n_chains    number of OUTPUT chains actually written. With --merge-groups
                that is the number of GROUPS, so it is smaller than the number of
                source chains ('A+C:A,B+D:B' on a 4-chain file writes 2).
    n_res       residues actually written, summed over the output chains. A merge
                does NOT change it: a merged chain contributes all of its
                segments' residues, so 41 + 49 reads as 90 on one chain.
    chains      the chain IDs actually written, comma-joined and in output order
                (AFTER --rename-to or a group's ':<out_id>', so this is what a
                viewer, usalign or ringfit --target-chains will see)
    n_atoms     heavy atoms actually written, summed the same way
    chain_n_res     SPLIT MODE ONLY (--split-groups), else NA: residues written
                    per output chain, in output-ID order -- 'A=152,...,D=60'
    chain_mean_b    SPLIT MODE ONLY, else NA: mean B-factor over the written atoms
                    of each output chain, 2 decimals -- 'A=1.00,...,D=0.00'. With
                    chain_n_res, the trust columns that say which segment is the
                    binder (RFdiffusion: B 0) and which the target (B 1)

Every one of those is "actually written", not "requested": a design whose requested
chain was absent gets an 'error: chain E not in ...' status, an empty path and NA
columns. That is the contract -- a shorter file would be indistinguishable from a
correct one downstream.

NA (empty) rather than 0 marks a field that does not apply: the extraction failed, or
no TSV was written at all. Do not read a blank as "zero chains" -- and do not read a
blank path as "nothing to extract".

Usage:
    sapia collect chainsel outputs/20260928_100140_7ojg_binder --table table1
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (chainsel_<name>).
INT_COLUMNS = ["n_chains", "n_res", "n_atoms"]
STR_COLUMNS = ["chains", "chain_n_res", "chain_mean_b"]
RESULT_COLUMNS = INT_COLUMNS + STR_COLUMNS


def _empty() -> dict[str, Any]:
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _fields(row: pd.Series) -> dict[str, Any]:
    """One result row -> collected columns, NA-safe (a missing or blank cell is a
    field that did not apply, not a zero)."""
    data = _empty()
    for col in INT_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = int(float(value))
    for col in STR_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = str(value)
    return data


def collect_chainsel(ctx: CollectCtx) -> CollectEach:
    """Per-design chainsel collector. The framework iterates ready designs and stamps
    status/path (keyed by the tool leaf); this reads one design's one-row <name>.tsv
    and adds the extraction columns. Variants -- a different chain set off the same
    table -- are distinguished via --dir-label, matching the output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"

        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        # Text columns as str. A non-split TSV has no chain_n_res/chain_mean_b
        # column at all; pandas ignores dtype keys for absent columns, and
        # _present() turns the absence into NA.
        result_df = pd.read_csv(
            tsv_path, sep="\t", dtype={c: str for c in STR_COLUMNS}
        )
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row_data = result_df.iloc[0]
        status = str(row_data["status"])

        if status != "OK":
            # An absent chain, an empty chain, an unreadable file: no path, no
            # counts. Leaving the counts at 0 would read as a successful
            # extraction of nothing.
            yield Collected(status=status, path="", data=_empty())
            return

        out_path = row_data.get("path")
        yield Collected(
            status=status,
            path="" if out_path is None or pd.isna(out_path) else str(out_path),
            data=_fields(row_data),
        )

    return one
