#!/usr/bin/env python3
"""
Collect posrmsd results into the table.

Reads ``<run_dir>/<table>/posrmsd[_<label>]/<prefix>/<name>.tsv`` (written by the
task, or at submit time for a design that could not be submitted) and writes
``posrmsd[_<label>]_<prefix>_<column>``. Like usalign, the status is
PER-COMPARISON (``..._<prefix>_status``), carried in ``data`` with ``status=None``
so the framework's bare ``<leaf>_status`` is not stamped.

STATUS VALUES
    OK                 every requested position resolved in both structures
    warn: ...          numbers computed, but the chains differ in length, so the
                       same-ordering premise is not guaranteed; excluded by == "OK"
    error: ...         nothing computed (names the structure, chain and index)
    missing_input: ... not submitted: --col-b or --positions-column resolved to
                       nothing on this row and all its ancestors
    missing            no result file at all (not selected by --col-a, not yet
                       finished, or killed before it wrote)

TRUST COLUMNS -- read before any RMSD:
    a_chain, b_chain, a_len, b_len, n_positions, seq_match_frac,
    a_motif_seq, b_motif_seq, a_resnum_offset, b_resnum_offset,
    a_gapless, b_gapless
SCIENCE:
    rmsd, max_dev, n_segments, worst_seg_rmsd, seg_ranges, seg_rmsds,
    seg0_rmsd, seg1_rmsd, ... (only for indices some design has)
    sup_path            structure A moved into B's frame (mmCIF)

Every non-OK/non-warn row carries NA (never 0) in every metric column, including
per-segment columns left over from an earlier collect of the same comparison.

Usage:
    sapia collect posrmsd outputs/RUN -t table1 \\
        --col-a boltz_path --col-b chainsel_binder_path
"""

import re
from argparse import ArgumentParser
from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectArgs, CollectCtx, CollectEach, Collected, DesignCtx

FLOAT_COLUMNS = ["seq_match_frac", "rmsd", "max_dev", "worst_seg_rmsd"]
INT_COLUMNS = [
    "a_len",
    "b_len",
    "n_positions",
    "a_resnum_offset",
    "b_resnum_offset",
    "a_gapless",
    "b_gapless",
    "n_segments",
]
STR_COLUMNS = [
    "a_chain",
    "b_chain",
    "a_motif_seq",
    "b_motif_seq",
    "seg_ranges",
    "seg_rmsds",
]
METRIC_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS
SEG_COLUMN = re.compile(r"^seg(\d+)_rmsd$")
PREFIX_RE = re.compile(r"^[A-Za-z0-9_]+$")


class PosRmsdCollectArgs(CollectArgs):
    col_a: str | None
    col_b: str | None
    output_prefix: str | None


def add_collect_posrmsd_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--col-a",
        type=str,
        default=None,
        help="Same --col-a as the run (used to derive the prefix when "
        "--output-prefix is not set).",
    )
    parser.add_argument(
        "--col-b",
        type=str,
        default=None,
        help="Same --col-b as the run (used to derive the prefix).",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="Same --output-prefix as the run. Defaults to "
        "'<col_a minus _path>_vs_<col_b minus _path>'.",
    )


def resolve_prefix(col_a: str | None, col_b: str | None, output_prefix: str | None) -> str:
    """Identical to run_posrmsd.resolve_prefix (kept in step by a test)."""
    if output_prefix:
        prefix = output_prefix.strip()
    elif col_a and col_b:
        prefix = f"{col_a.strip().removesuffix('_path')}_vs_{col_b.strip().removesuffix('_path')}"
    else:
        raise ValueError("Provide --output-prefix, or both --col-a and --col-b.")
    if not PREFIX_RE.match(prefix):
        raise ValueError(
            f"comparison prefix {prefix!r} must be letters, digits and '_' only "
            f"(it names a directory and every column); pass --output-prefix"
        )
    return prefix


def _present(row: pd.Series, col: str) -> str | None:
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return str(value).strip()


def collect_posrmsd(ctx: CollectCtx) -> CollectEach:
    """Per-design posrmsd collector for one named comparison (prefix)."""
    prefix = resolve_prefix(
        getattr(ctx.args, "col_a", None),
        getattr(ctx.args, "col_b", None),
        getattr(ctx.args, "output_prefix", None),
    )
    results_dir = ctx.out_dir / prefix
    if not results_dir.is_dir():
        raise FileNotFoundError(
            f"Results dir not found: {results_dir}. collect needs the same "
            f"--col-a/--col-b or --output-prefix (and -l) as the run."
        )

    def key(col: str) -> str:
        return f"{prefix}_{col}"

    # Per-segment columns this comparison already has in the table: NA'd on every
    # row unless re-filled, so a re-collect never leaves a stale segment value
    # beside a new status.
    leaf = ctx.out_dir.name
    stale_seg = re.compile(rf"^{re.escape(leaf)}_{re.escape(prefix)}_(seg\d+_rmsd)$")
    existing_seg = [
        m.group(1) for c in map(str, ctx.df.columns) if (m := stale_seg.match(c))
    ]

    def empty() -> dict[str, Any]:
        data: dict[str, Any] = {key(c): None for c in METRIC_COLUMNS}
        data.update({key(c): None for c in existing_seg})
        data[key("sup_path")] = None
        return data

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = results_dir / f"{d.name}.tsv"
        data = empty()
        if not tsv_path.is_file():
            data[key("status")] = "missing"
            yield Collected(status=None, data=data)
            return

        # Everything as text: no NA-sniffing of a chain or sequence that happens
        # to spell 'NA', and explicit conversion below.
        result_df = pd.read_csv(
            tsv_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""]
        )
        if result_df.empty:
            data[key("status")] = "error: empty tsv"
            yield Collected(status=None, data=data)
            return

        row = result_df.iloc[0]
        status = str(row["status"])
        data[key("status")] = status
        if status != "OK" and not status.startswith("warn:"):
            yield Collected(status=None, data=data)
            return

        for col in FLOAT_COLUMNS:
            if (v := _present(row, col)) is not None:
                data[key(col)] = float(v)
        for col in INT_COLUMNS:
            if (v := _present(row, col)) is not None:
                data[key(col)] = int(float(v))
        for col in STR_COLUMNS:
            if (v := _present(row, col)) is not None:
                data[key(col)] = v
        for col in row.index:
            if SEG_COLUMN.match(str(col)) and (v := _present(row, str(col))) is not None:
                data[key(str(col))] = float(v)
        data[key("sup_path")] = _present(row, "sup_path")
        yield Collected(status=None, data=data)

    return one
