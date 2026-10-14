#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) measuring, per design, the CA RMSD of structure A
onto structure B over an EXPLICIT list of chain positions.

PREMISE. A and B carry the same residue ordering on the compared chains, so the
correspondence is "position N in A = position N in B" -- nothing is aligned or
searched for (see ``posrmsd_worker.py`` for the full premise and for when the
numbers mean nothing). B is the REFERENCE (fixed frame); A is MOVED onto it.

ROWS AND LINEAGE -- mirrors the built-in ``usalign``:

* ``--col-a`` selects the rows: every row of ``-t`` whose ``--col-a`` cell is
  present (minus rows already ``OK`` under this comparison, unless ``--force``).
  ``-i/--input-column`` is not used.
* ``--col-b`` and ``--positions-column`` are resolved per design with
  ``ctx.lookup`` -- the same lineage walk usalign uses for ``--col-b``: the row's
  own value if set, else the nearest ancestor table's (``parent_name`` /
  ``parent_table``). So a Boltz prediction on a child table is compared with the
  backbone column, and the motif positions, that live on its parent.
* A row on which either resolves to nothing is NOT submitted, and is NOT dropped
  silently either: a per-design TSV with status
  ``missing_input: <which flag> <column> has no value on this row or any
  ancestor`` is written straight into the results dir here, so it collects as a
  named situation, distinct from ``missing`` (no output at all) and ``error:``.
  If a column resolves on NONE of the selected rows the submit is refused -- that
  is a typo or a wrong table, not per-design data.

POSITIONS. ``--positions`` (one literal spec for every row) and
``--positions-column`` (a per-design spec from a lineage-resolved column) are
mutually exclusive; exactly one is required. The grammar is prosapia's position
mini-language, parsed HERE with prosapia's own ``parse_chain_positions`` /
``parse_positions`` (``{expr}`` islands resolved up the lineage): 1-indexed,
inclusive ``start:end`` ranges and single positions, comma-separated, ONE chain
group -- a ``/`` is refused. Each comma-separated fragment is one SEGMENT. The
result is handed to the task as JSON, so the task never parses the grammar.
Refused as well: open-ended ranges (no chain length is known at submit), a
fragment that expands to nothing (``18:11``), a position < 1, a position listed
twice, and fewer than 3 positions in total. A literal brace-free ``--positions``
that fails these checks refuses the whole submit; a per-design spec that fails
them writes that design an ``error:`` row and the rest are submitted.

PREFIX CONTRACT (same as usalign). Every column is
``posrmsd[_<dir-label>]_<prefix>_<metric>`` with ``<prefix>`` from
``--output-prefix`` or ``<col_a minus _path>_vs_<col_b minus _path>``; the
per-design files live in ``<out_dir>/<prefix>/``. ``sapia collect posrmsd`` must
be given the same ``--col-a``/``--col-b`` or ``--output-prefix``.

Usage:
    sapia run posrmsd outputs/RUN -t table1 \\
        --col-a boltz_path --col-b chainsel_binder_path \\
        --positions-column motifrmsd_design_fixed_positions \\
        --chain-a B
"""

import csv
import json
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, cast

import pandas as pd
from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.utils import parse_chain_positions, parse_positions, resolve_template

MIN_POSITIONS = 3
# spec.py's default_input_column: matches no column (as for usalign), so the
# framework's filter_ready leaves the frame unchanged on both run and collect.
NO_INPUT_COLUMN = "not applicable"
# A chain id as gemmi reports it. Anything else (a comma list, a range) cannot
# name one chain, so it is refused at submit.
CHAIN_RE = re.compile(r"^[A-Za-z0-9]+$")
# The prefix names a directory and is embedded in every column name.
PREFIX_RE = re.compile(r"^[A-Za-z0-9_]+$")


class PosRmsdArgs(CommonArgs):
    col_a: str
    col_b: str
    positions: str
    positions_column: str
    chain_a: str
    chain_b: str
    output_prefix: str | None


def add_run_posrmsd_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--col-a",
        type=str,
        required=True,
        help="Column holding structure A -- the design/prediction, the structure "
        "that is MOVED onto B. Selects which rows run (a row with no value here is "
        "not submitted). Read from the -t table itself, not up the lineage.",
    )
    parser.add_argument(
        "--col-b",
        type=str,
        required=True,
        help="Column holding structure B -- the REFERENCE/backbone, the fixed "
        "frame. Resolved per design up the lineage (row's own value, else the "
        "nearest ancestor's), like usalign's --col-b.",
    )
    parser.add_argument(
        "--positions",
        type=str,
        default="",
        help="One position spec for every row: 1-indexed, inclusive, "
        "comma-separated ranges in FILE ORDER within one chain, e.g. "
        "'11:18,45:54,87:92' (prosapia's position language; {expr} islands "
        "allowed). Each fragment is one segment. Mutually exclusive with "
        "--positions-column; exactly one is required.",
    )
    parser.add_argument(
        "--positions-column",
        type=str,
        default="",
        help="Column holding a per-design position spec (same grammar as "
        "--positions), e.g. motifrmsd_design_fixed_positions. Resolved up the "
        "lineage. Mutually exclusive with --positions.",
    )
    parser.add_argument(
        "--chain-a",
        type=str,
        default="",
        help="Chain of A to index into. Empty (default): A must hold exactly one "
        "polymer chain, otherwise that design is an error naming the chains -- a "
        "chain is never guessed.",
    )
    parser.add_argument(
        "--chain-b",
        type=str,
        default="",
        help="Chain of B to index into; same rule as --chain-a.",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="Prefix for the results subdir and every column. Defaults to "
        "'<col_a minus _path>_vs_<col_b minus _path>'. collect must be given the "
        "same value (or the same --col-a/--col-b).",
    )


def resolve_prefix(col_a: str | None, col_b: str | None, output_prefix: str | None) -> str:
    """The comparison prefix, shared verbatim with the collector."""
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


def _check_chain(value: str, flag: str) -> str:
    value = (value or "").strip()
    if value and not CHAIN_RE.match(value):
        raise ValueError(
            f"{flag} {value!r} is not a single chain id (letters/digits only, e.g. "
            f"'A'); each structure is indexed on exactly one chain"
        )
    return value


def parse_segments(spec: str, lookup: Any, name: str) -> tuple[list[str], list[list[int]]]:
    """Split a one-chain position spec into its comma-separated SEGMENTS.

    Each fragment is expanded by prosapia's own ``parse_chain_positions`` -- the
    function ``parse_positions`` (and so ProteinMPNN's ``--fixed-positions``)
    delegates to -- after the same ``resolve_template`` / ``[...]`` handling
    ``parse_positions`` applies. The concatenation is then checked against
    ``parse_positions`` on the original spec, so the semantics cannot drift from
    prosapia's: a disagreement raises.
    """
    resolved = resolve_template(spec, lookup, name).strip().strip("[]").strip()
    if not resolved:
        raise ValueError(f"design {name!r}: empty position spec {spec!r}")
    if "/" in resolved:
        raise ValueError(
            f"design {name!r}: position spec {spec!r} has a '/' chain break; "
            f"posrmsd addresses ONE chain per structure (--chain-a/--chain-b), so "
            f"give one chain's positions only"
        )
    labels: list[str] = []
    segments: list[list[int]] = []
    for token in resolved.split(","):
        token = token.strip()
        if not token:  # parse_chain_positions skips empty tokens too
            continue
        positions = parse_chain_positions(token, name=name)
        if not positions:
            raise ValueError(
                f"design {name!r}: fragment {token!r} of {spec!r} expands to no "
                f"positions (start > end?)"
            )
        if min(positions) < 1:
            raise ValueError(
                f"design {name!r}: fragment {token!r} of {spec!r} has a position "
                f"< 1; positions are 1-indexed"
            )
        labels.append(
            str(positions[0]) if len(positions) == 1 else f"{positions[0]}:{positions[-1]}"
        )
        segments.append(positions)

    flat = [p for seg in segments for p in seg]
    reference = parse_positions(spec, lookup, name)
    if reference != [flat]:
        raise RuntimeError(
            f"design {name!r}: posrmsd's segment split of {spec!r} ({flat}) "
            f"disagrees with prosapia.parse_positions ({reference}); refusing"
        )
    dupes = sorted({p for p in flat if flat.count(p) > 1})
    if dupes:
        raise ValueError(
            f"design {name!r}: position spec {spec!r} lists position(s) "
            f"{', '.join(map(str, dupes[:6]))} more than once (overlapping "
            f"fragments); one CA would enter the fit twice"
        )
    if len(flat) < MIN_POSITIONS:
        raise ValueError(
            f"design {name!r}: position spec {spec!r} gives {len(flat)} "
            f"position(s); a superposition needs at least {MIN_POSITIONS}"
        )
    return labels, segments


def _no_lookup(name: str, column: str) -> Any:
    raise ValueError(f"no table lookup available for column {column!r}")


def _lookup_str(lookup: Any, name: str, column: str) -> str | None:
    """A lineage-resolved cell as a stripped string, or None when absent/blank."""
    value = lookup(name, column)
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return None
    text = str(value).strip()
    return text or None


def _write_status_tsv(path: Path, name: str, status: str) -> None:
    """A per-design result written at submit time (no task ran). Same minimal
    shape as the task script's crash fallback; the collector NAs every metric."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["name", "status", "sup_path"])
        writer.writerow([name, status.replace("\t", " ").replace("\n", " "), ""])


def build_posrmsd_manifest(ctx: ManifestCtx[PosRmsdArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool

    args = ctx.args
    col_a = (args.col_a or "").strip()
    col_b = (args.col_b or "").strip()
    pos_literal = (args.positions or "").strip()
    pos_column = (args.positions_column or "").strip()

    if pos_literal and pos_column:
        raise ValueError(
            "--positions (one spec for every row) and --positions-column (a "
            "per-design spec from a table column) are mutually exclusive; pass "
            "exactly one."
        )
    if not pos_literal and not pos_column:
        raise ValueError("Pass exactly one of --positions or --positions-column.")
    if not col_a or not col_b:
        raise ValueError("--col-a and --col-b are both required.")
    if (getattr(args, "input_column", NO_INPUT_COLUMN) or NO_INPUT_COLUMN) != NO_INPUT_COLUMN:
        raise ValueError(
            f"posrmsd selects rows by --col-a; -i/--input-column "
            f"({args.input_column!r}) is not used and would make collect iterate a "
            f"different row set. Drop -i."
        )
    chain_a = _check_chain(args.chain_a, "--chain-a")
    chain_b = _check_chain(args.chain_b, "--chain-b")
    prefix = resolve_prefix(col_a, col_b, args.output_prefix)

    df = ctx.df
    if col_a not in df.columns:
        available = ", ".join(sorted(map(str, df.columns))) or "<none>"
        raise ValueError(
            f"--col-a {col_a!r} is not a column of the -t table (it selects the "
            f"rows, so it is not resolved up the lineage). Available: {available}"
        )

    # A brace-free literal is the same for every row: validate it once, here,
    # so a typo refuses the submit instead of writing N identical error rows.
    literal_segments = None
    if pos_literal and "{" not in pos_literal:
        literal_segments = parse_segments(pos_literal, _no_lookup, "<--positions>")

    cells = df[col_a]
    ready = df[cells.notna() & (cells.astype(str).str.strip() != "")]
    # The driver leaf-prefixes collected data, so this comparison's status lives
    # at <leaf>_<prefix>_status (as for usalign).
    status_col = f"{ctx.out_dir.name}_{prefix}_status"
    if args.force:
        print("Re-running posrmsd for all selected designs (including OK).")
    elif status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    results_dir = ctx.out_dir / prefix

    # Resolve everything first, write nothing until the whole batch is known to
    # be sane: a column that resolves on no row at all refuses the submit.
    plan: list[tuple[str, str, str | None, str | None]] = []
    for name in ready.index:
        name = cast(str, name)
        src_a = str(ready.at[name, col_a]).strip()
        src_b = _lookup_str(ctx.lookup, name, col_b)
        spec = pos_literal if pos_literal else _lookup_str(ctx.lookup, name, pos_column)
        plan.append((name, src_a, src_b, spec))

    if plan:
        if all(b is None for _, _, b, _ in plan):
            raise ValueError(
                f"--col-b {col_b!r} has no value on any of the {len(plan)} selected "
                f"row(s) or their ancestors -- a typo, or the column lives on a "
                f"table that is not an ancestor of -t. Nothing submitted."
            )
        if pos_column and all(s is None for _, _, _, s in plan):
            raise ValueError(
                f"--positions-column {pos_column!r} has no value on any of the "
                f"{len(plan)} selected row(s) or their ancestors. Nothing submitted."
            )

    manifest_rows: list[tuple[str, ...]] = []
    n_missing = n_error = 0
    for name, src_a, src_b, spec in plan:
        result_tsv = results_dir / f"{name}.tsv"
        missing = []
        if src_b is None:
            missing.append(f"--col-b {col_b!r} has no value on this row or any ancestor")
        if spec is None:
            missing.append(
                f"--positions-column {pos_column!r} has no value on this row or any ancestor"
            )
        if missing:
            n_missing += 1
            print(f"{name}: missing_input -- {'; '.join(missing)}")
            _write_status_tsv(result_tsv, name, "missing_input: " + "; ".join(missing))
            continue

        assert spec is not None and src_b is not None
        try:
            labels, segments = (
                literal_segments
                if literal_segments is not None
                else parse_segments(spec, ctx.lookup, name)
            )
        except ValueError as e:
            n_error += 1
            print(f"{name}: position spec refused -- {e}")
            _write_status_tsv(result_tsv, name, f"error: position spec: {e}")
            continue

        # A stale result from an earlier run must not survive a task that is
        # killed before it writes: remove it, so that design collects 'missing'.
        result_tsv.unlink(missing_ok=True)
        (results_dir / f"{name}_sup.cif").unlink(missing_ok=True)

        payload = json.dumps({"labels": labels, "segments": segments}, separators=(",", ":"))
        manifest_rows.append(
            (
                name,
                str(volume_path(src_a)),
                str(volume_path(src_b)),
                chain_a,  # may be empty: a mid-line empty field keeps its place
                chain_b,
                payload,
                prefix,  # kept last: always non-empty
            )
        )

    if n_missing or n_error:
        print(
            f"WARNING: {n_missing} design(s) not submitted for missing input and "
            f"{n_error} for an unusable position spec; each has a status row in "
            f"{results_dir} and will collect as such."
        )
    return manifest_rows
