#!/usr/bin/env python3
"""
Collect mkcomplex results into the table.

Reads the per-design TSVs mkcomplex.sh wrote under <run_dir>/<table>/mkcomplex/ and adds
<prefix>_status plus:

    sequence            the assembled multi-chain string, e.g. 'SEQA/SEQB/DESIGN'
                        (feed it to a predictor with ``-i mkcomplex_sequence``)
    n_chains            number of separator-delimited chains
    chain_lens          per-chain residue counts, comma-joined, in chain order
    design_chain_index  0-based index of the design's chain (the FIRST copy under
                        ``--design-copies N``)
    total_len           residues summed over all chains
    design_copies       copies of the design's sequence emitted (``--design-copies``)

<prefix>_path is always empty (no structure is written). A design whose sequence is
empty, already multi-chain or carries a non-standard residue (in a fixed chain or its
own) gets an 'error: ...' status and empty columns; the task still exited 0, so the
rest of the array collects normally.

``design_copies`` is read from the TSV, never inferred: a run made before the flag
existed re-collects as NA, not 1, since an invented 1 would read as a measured one.

Usage:
    sapia collect mkcomplex outputs/20260928_100140_7ojg_binder --table table1
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (mkcomplex_<name>).
INT_COLUMNS = ["n_chains", "design_chain_index", "total_len", "design_copies"]
STR_COLUMNS = ["sequence", "chain_lens"]
RESULT_COLUMNS = STR_COLUMNS + INT_COLUMNS


def _empty() -> dict[str, Any]:
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _fields(row: pd.Series) -> dict[str, Any]:
    """One result row -> columns; a blank cell is a missing field, not a zero."""
    data = _empty()
    for col in STR_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = str(value)
    for col in INT_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = int(float(value))
    return data


def collect_mkcomplex(ctx: CollectCtx) -> CollectEach:
    """Per-design collector: adds the assembled-sequence columns from ``<name>.tsv``."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"

        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(tsv_path, sep="\t", dtype={"sequence": str})
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row_data = result_df.iloc[0]
        # A TSV without a status column is truncated/foreign: error row, don't raise.
        status = str(row_data.get("status", "error: unknown"))

        if status != "OK":
            yield Collected(status=status, path="", data=_empty())
            return

        # No structure is produced: the result is the sequence column, not a path.
        yield Collected(status=status, path="", data=_fields(row_data))

    return one
