#!/usr/bin/env python3
"""
Submit ONE task that clusters a whole set of design sequences with MMseqs2
``easy-cluster``, and annotate every design with the cluster it landed in.

PREMISE AND SCOPE. The question is "how much sequence redundancy is in this set of
designs?". The input is ONE amino-acid sequence per row (one chain, one-letter code),
in a column of THIS table. The answer is a property of the SET, not of any one design:
a design's cluster depends on every other sequence that was clustered with it. So:

  * The numbers mean nothing outside the set they were computed on. Add or remove one
    design and representatives, cluster IDs and sizes can all change.
    ``seqcluster_cluster_id`` in particular is a label local to one run (1-based, by
    first appearance of the representative in MMseqs2's output), NOT a portable or
    hashed identifier. Never compare IDs across two runs or two dir labels.
  * The numbers are meaningless if the column is not the sequence you mean to compare:
    a multi-chain complex string (``/`` or ``:`` separated -- ProteinMPNN joins the
    chains of a multi-chain design with ``/``) is REFUSED here, because a target chain
    shared by every design would make every design look alike. Pass a single-chain
    (binder) sequence column.
  * No structure, no fold, no function: two designs in different clusters can share a
    fold, and two in one cluster can differ at every interface position that matters.
  * It records, it does not gate. The MMseqs2 parameters default to MMseqs2's own
    defaults for ``easy-cluster`` (min-seq-id 0.0, -c 0.8, cov-mode 0, cluster-mode 0,
    as printed by ``mmseqs easy-cluster -h`` for release 18-8cc5c); there is no
    "correct" threshold baked in.

THE ONE-TASK ARCHITECTURE. Every other tool in this repo shards its ready set into
independent tasks. This one deliberately does not: the builder returns EXACTLY ONE
manifest row, whose first field points at a members file listing every ready design's
name and sequence. There is no ``--designs-per-task`` flag, because there is no valid
value other than "all of them". (A manifest row is a task -- see
``prosapia.core.base_run.run_from_args``: ``n = len(manifest_rows)`` tasks are
submitted -- so returning one row IS the framework's way to say "one task". The driver
therefore prints ``Submitting 1 designs`` however many designs are clustered.)

THE READY SET, AND WHY THERE IS NO PARTIAL RESUME. The ready set is every row of the
(``-f``-filtered) table whose input column is non-null and non-blank. It is NOT
``ctx.ready``: that property drops rows this tool already marked ``OK``, and
re-clustering only the not-yet-OK rows would cluster a SUBSET on its own and stamp it
with cluster IDs that collide with the earlier run's. So every run re-clusters the
whole ready set, ``--force`` or not, and the collector re-annotates every row (rows
not in this run's set are reset to ``missing`` with NA, so an earlier run's cluster
labels never survive next to this run's). To keep two clusterings side by side, use
two ``-l/--dir-label`` values.

One case the collector CANNOT reset: a row that carries an earlier run's labels but
has no value in THIS run's input column (the collector only ever visits rows with a
present input column). That happens when ``-i`` changes under the same dir label. The
builder refuses it up front, naming the rows, rather than leave two clusterings mixed
under one set of columns.

``action: update`` -- a cluster assignment is a property of designs that already
exist; this annotates the table in place.

``default_input_column`` is the ``"not applicable"`` sentinel (as in chainsel,
ssprofile, epitope): there is no generic sequence column, so the builder raises unless
``-i/--input-column`` names one the table actually has. The column is read from THIS
table only (no lineage walk): the designs being clustered are this table's rows, and a
sequence inherited from a parent would be shared by all its children.

Usage:
    sapia run seqcluster outputs/20261010_binders \
        --table table1 \
        --input-column proteinmpnn_redesign_full_sequence
"""

from __future__ import annotations

import json
import uuid
from argparse import ArgumentParser
from typing import cast

import pandas as pd

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.core.naming import status_column

# The sentinel default_input_column (see the module docstring).
NO_DEFAULT_COLUMN = "not applicable"

# On-disk layout. Kept in step BY HAND with seqcluster_worker.py and
# collect_seqcluster.py (the worker runs standalone in the task image and the
# neighbours' run/collect modules never import siblings); test_seqcluster.py asserts
# the three agree.
TASK_DIRNAME = "seqcluster_task"
MEMBERS_FILENAME = "members.tsv"
RUN_FILENAME = "run.json"
RESULT_DIRNAME = "seqcluster_result"
ASSIGNMENTS = "assignments.tsv"

# MMseqs2 easy-cluster defaults, release 18-8cc5c (read from `mmseqs easy-cluster -h`).
# The tool ALWAYS passes these explicitly, so the recorded params never depend on
# whatever defaults a different binary version might have.
MMSEQS_DEFAULTS = {
    "min_seq_id": 0.0,
    "coverage": 0.8,
    "cov_mode": 0,
    "cluster_mode": 0,
}

# One-letter amino-acid code MMseqs2 accepts (20 standard + X unknown + the ambiguity /
# rare codes B Z J U O). Anything else -- a chain separator, a gap, a stop -- is refused.
AA_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYXBZJUO")
# Chain separators: ProteinMPNN joins multi-chain designs with '/'; ':' is the
# ColabFold/Boltz-style complex separator.
CHAIN_SEPARATORS = frozenset("/:")


class SeqclusterArgs(CommonArgs):
    min_seq_id: float
    coverage: float
    cov_mode: int
    cluster_mode: int


def add_run_seqcluster_args(parser: ArgumentParser) -> None:
    d = MMSEQS_DEFAULTS
    parser.add_argument(
        "--min-seq-id",
        type=float,
        default=d["min_seq_id"],
        help="MMseqs2 --min-seq-id: minimum sequence identity of an alignment for "
        f"two sequences to cluster (0.0-1.0). Default {d['min_seq_id']} (MMseqs2's own "
        "default, i.e. no identity floor: clustering is then decided by -c and the "
        "E-value alone).",
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=d["coverage"],
        help="MMseqs2 -c: minimum fraction of aligned (covered) residues, interpreted "
        f"per --cov-mode (0.0-1.0). Default {d['coverage']} (MMseqs2's own default).",
    )
    parser.add_argument(
        "--cov-mode",
        type=int,
        choices=range(0, 6),
        default=d["cov_mode"],
        help="MMseqs2 --cov-mode: 0 coverage of query and target; 1 of target; 2 of "
        "query; 3 target length >= c x query length; 4 query length >= c x target "
        "length; 5 shorter length >= c x longer length. Default "
        f"{d['cov_mode']} (MMseqs2's own default).",
    )
    parser.add_argument(
        "--cluster-mode",
        type=int,
        choices=range(0, 4),
        default=d["cluster_mode"],
        help="MMseqs2 --cluster-mode: 0 Set-Cover (greedy); 1 connected component "
        "(BLASTclust-like, transitive); 2,3 greedy by sequence length (CD-HIT-like). "
        f"Default {d['cluster_mode']} (MMseqs2's own default).",
    )


def _validate_params(args: SeqclusterArgs) -> None:
    if not 0.0 <= args.min_seq_id <= 1.0:
        raise ValueError(f"--min-seq-id must be within 0.0-1.0, got {args.min_seq_id}.")
    if not 0.0 <= args.coverage <= 1.0:
        raise ValueError(f"--coverage must be within 0.0-1.0, got {args.coverage}.")


def manifest_params(args: SeqclusterArgs) -> tuple[str, str, str, str]:
    """The four MMseqs2 values exactly as they go into the manifest (and the worker)."""
    return (
        repr(float(args.min_seq_id)),
        repr(float(args.coverage)),
        str(int(args.cov_mode)),
        str(int(args.cluster_mode)),
    )


def params_string(min_seq_id: str, coverage: str, cov_mode: str, cluster_mode: str) -> str:
    """The provenance string recorded in ``<leaf>_params``. seqcluster_worker.py builds
    the same string from the same four manifest fields; the collector requires the two
    to agree, so a result produced under different parameters can never be collected
    as this run's."""
    return f"min_seq_id={min_seq_id};c={coverage};cov_mode={cov_mode};cluster_mode={cluster_mode}"


def _sequence_columns(df: pd.DataFrame) -> list[str]:
    return [str(c) for c in df.columns if "seq" in str(c).lower()]


def clean_sequence(raw: object) -> str:
    """Strip surrounding whitespace and upper-case. Does NOT drop inner characters."""
    return str(raw).strip().upper()


def bad_characters(seq: str) -> set[str]:
    return set(seq) - AA_ALPHABET


def _preview(items: list[str]) -> str:
    return ", ".join(items[:5]) + (" ..." if len(items) > 5 else "")


def build_seqcluster_manifest(ctx: ManifestCtx[SeqclusterArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0
    _validate_params(ctx.args)

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN or column not in ctx.df.columns:
        candidates = ", ".join(_sequence_columns(ctx.df)) or "(none)"
        what = (
            "seqcluster has no default input column"
            if column == NO_DEFAULT_COLUMN
            else f"column {column!r} is not in table '{ctx.args.table}'"
        )
        raise ValueError(
            f"{what}: pass -i/--input-column with a single-chain amino-acid sequence "
            f"column of THIS table (no lineage walk). Sequence-like columns here: "
            f"{candidates}."
        )

    # The ready set: every (filtered) row with a present sequence. NOT ctx.ready, which
    # drops already-OK rows -- a set-level clustering cannot be resumed on a subset.
    col = ctx.df[column]
    present = col.notna() & (col.astype(str).str.strip() != "")
    ready = ctx.df[present]

    status_col = status_column(ctx.out_dir.name)
    if status_col in ctx.df.columns:
        prior = ctx.df[status_col]
        # Rows labelled by an earlier run (anything but 'missing', which the collector
        # writes with every data column NA) that this run's collector will never visit,
        # because they have no value in THIS input column. Their old labels would sit
        # next to this run's, under the same columns, with colliding cluster IDs.
        stale = ctx.df.index[prior.notna() & (prior != "missing") & ~present]
        if len(stale):
            raise ValueError(
                f"{len(stale)} row(s) carry results in '{status_col}' from an earlier "
                f"seqcluster run but have no {column!r}: {_preview([str(s) for s in stale])}. "
                f"The collector only visits rows with a present input column, so those "
                f"labels would survive next to this run's and their cluster IDs would "
                f"collide. This happens when -i changes under the same dir label: use a "
                f"different -l/--dir-label. Nothing was submitted."
            )
        if (prior == "OK").any() and not ctx.args.force:
            print(
                "seqcluster: rows already OK from an earlier run will be RE-CLUSTERED with "
                "the rest -- a set-level clustering is never resumed on a subset (--force "
                "is implied). Use a different -l/--dir-label to keep the earlier result."
            )

    n_absent = int((~present).sum())
    if n_absent:
        print(
            f"seqcluster: {n_absent} row(s) have no {column}; they are not clustered "
            f"and will collect as 'missing'."
        )

    members: list[tuple[str, str]] = []
    multichain: list[str] = []
    bad: list[str] = []
    for name in ready.index:
        name = cast(str, name)
        if "\t" in name or "\n" in name:
            raise ValueError(f"design name {name!r} contains a tab or newline.")
        seq = clean_sequence(ready.at[name, column])
        if CHAIN_SEPARATORS & set(seq):
            multichain.append(name)
            continue
        if chars := bad_characters(seq):
            bad.append(f"{name} ({''.join(sorted(chars))!r})")
            continue
        members.append((name, seq))

    if multichain:
        raise ValueError(
            f"{len(multichain)} design(s) in {column!r} hold a MULTI-CHAIN sequence (a '/' "
            f"or ':' chain separator): {_preview(multichain)}. ProteinMPNN joins the chains "
            f"of a multi-chain design with '/', so a complex column carries the target "
            f"chain(s) too, and a target shared by every design would make every design "
            f"cluster together. seqcluster clusters ONE chain per design and does not "
            f"split complexes: pass a single-chain (binder) sequence column. Nothing was "
            f"submitted."
        )
    if bad:
        raise ValueError(
            f"{len(bad)} design(s) in {column!r} contain characters outside the "
            f"one-letter amino-acid code (ACDEFGHIKLMNPQRSTVWY plus X B Z J U O): "
            f"{_preview(bad)}. Nothing was submitted."
        )

    task_dir = ctx.out_dir / TASK_DIRNAME
    if not members:
        # Leave no earlier run's members/run.json behind for collect to pair with an
        # earlier result: with no run set, every row collects as 'missing'.
        for f in (MEMBERS_FILENAME, RUN_FILENAME):
            (task_dir / f).unlink(missing_ok=True)
        return []

    # Exactly one task: a members file listing the whole set, and one manifest row.
    task_dir.mkdir(parents=True, exist_ok=True)
    members_file = task_dir / MEMBERS_FILENAME
    with open(members_file, "w") as f:
        for name, seq in members:
            f.write(f"{name}\t{seq}\n")
    token = uuid.uuid4().hex
    mparams = manifest_params(ctx.args)
    (task_dir / RUN_FILENAME).write_text(
        json.dumps(
            {
                "token": token,
                "params": params_string(*mparams),
                "input_column": column,
                "n_designs": len(members),
            },
            indent=2,
        )
        + "\n"
    )
    # Drop the previous run's result now, so a collect issued before this task
    # finishes cannot pair the new members with the old assignments. (The token,
    # echoed by the worker and checked by the collector, is the real guard: under
    # Modal the task container sees the Volume through its own mount.)
    (ctx.out_dir / RESULT_DIRNAME / ASSIGNMENTS).unlink(missing_ok=True)

    print(f"seqcluster: clustering {len(members)} designs in ONE task")
    return [
        (
            str(volume_path(members_file)),
            token,
            *mparams,
            # Always non-empty. New fields go on the END: seqcluster.sh cuts by index.
            column,
        )
    ]
