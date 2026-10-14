#!/usr/bin/env python3
"""
Collect a seqcluster run into the table.

seqcluster is SET-LEVEL: one task clustered every design of the run together, so this
collector validates ONE result for the whole run, once, and then annotates every row
from it. Files read, under <run_dir>/<table>/seqcluster[_<label>]/:

  seqcluster_task/members.tsv          the set the run submitted (name<TAB>sequence)
  seqcluster_task/run.json             the submit's token, params and input column
  seqcluster_result/assignments.tsv    the worker's one-row-per-member result

Columns (bare names; the driver leaf-prefixes them to seqcluster[_<label>]_<name>):

  TRUST COLUMNS -- read these first:
    n_designs_clustered  designs fed into the ONE combined task (same on every member
                         row). Compare with the number you meant to cluster: a filter
                         or a sparse sequence column makes it smaller than the table.
    seq_len              length of THIS design's sequence as clustered. A complex-
                         sized length means the wrong column. Kept on error rows.
    input_column         the table column the sequences came from  } provenance, the
    params               the MMseqs2 clustering parameters used     } same on every
    mmseqs_version       the MMseqs2 build (commit hash)            } row of a run

  SCIENCE COLUMNS:
    cluster_rep          table `name` of this design's cluster representative
    cluster_id           1-based integer per representative, by first appearance in
                         MMseqs2's cluster TSV. LOCAL TO THIS ONE RUN: a rerun on a
                         different set can renumber every cluster. Not a hash.
    cluster_size         designs (within this run's set) sharing that representative
    is_rep               True on the representative's own row (one per cluster)
    n_clusters           distinct clusters in this run's set (same on every row)

Statuses:
    'OK'          this design is a member of the validated clustering.
    'error: ...'  the combined task failed, has not finished, or its result does not
                  belong to this run / is internally inconsistent. EVERY member of the
                  run gets the SAME status (the task is all-or-nothing); the science
                  columns are NA, while seq_len / input_column / params are kept (they
                  come from this run's own members list and run.json, not from the
                  result being refused).
    'missing'     the row was not in this run's set (filtered out with -f, a blank
                  sequence, or no run has been submitted under this dir label) --
                  expected, not a failure. EVERY data column is written as NA, so an
                  earlier run's labels never survive on a row this run did not cluster.

RERUN SEMANTICS. A cluster label is meaningful only together with the labels of the
rest of its set, so a partial re-collect would mix two clusterings under one set of
columns. The factory therefore forces ``ctx.args.force = True`` before the driver
reads ``ctx.ready`` (``by_design`` calls the factory first, then iterates
``ctx.ready``; see prosapia/core/base_collect.py), and every collect re-annotates every
row with a present input column. ``--force`` on the command line changes nothing.

NA is never 0: a 0 or False here is a measurement; an empty cell is the absence of one.

Usage:
    sapia collect seqcluster outputs/20261010_binders --table table1
"""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Kept in step BY HAND with run_seqcluster.py / seqcluster_worker.py (test-asserted).
TASK_DIRNAME = "seqcluster_task"
MEMBERS_FILENAME = "members.tsv"
RUN_FILENAME = "run.json"
RESULT_DIRNAME = "seqcluster_result"
ASSIGNMENTS = "assignments.tsv"

# Every data column, bare. A row this run did not cluster gets every one of them NA.
TRUST_COLUMNS = ["n_designs_clustered", "seq_len", "input_column", "params", "mmseqs_version"]
SCIENCE_COLUMNS = ["cluster_rep", "cluster_id", "cluster_size", "is_rep", "n_clusters"]
DATA_COLUMNS = TRUST_COLUMNS + SCIENCE_COLUMNS
# The worker's assignments.tsv header, in order (seqcluster_worker.FIELDS), test-asserted.
ASSIGNMENT_FIELDS = [
    "name", "status", "token", "seq_len", "cluster_rep", "cluster_id", "cluster_size",
    "is_rep", "n_designs_clustered", "n_clusters", "mmseqs_version", "params", "input_column",
]


class _Refused(Exception):
    """The run's result cannot be collected: the message becomes every member's status."""


def _empty() -> dict[str, Any]:
    return {c: None for c in DATA_COLUMNS}


def _read_members(path: Path) -> dict[str, str]:
    """name -> sequence, in file order. Raises _Refused on a malformed file."""
    members: dict[str, str] = {}
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 2 or not parts[0] or not parts[1]:
            raise _Refused(f"error: malformed line {lineno} in {TASK_DIRNAME}/{MEMBERS_FILENAME}")
        if parts[0] in members:
            raise _Refused(f"error: duplicate name {parts[0]} in {TASK_DIRNAME}/{MEMBERS_FILENAME}")
        members[parts[0]] = parts[1]
    return members


def _int(value: str, what: str) -> int:
    try:
        f = float(value)
    except (TypeError, ValueError):
        raise _Refused(f"error: assignments.tsv {what}={value!r} is not a number") from None
    if not f.is_integer():
        raise _Refused(f"error: assignments.tsv {what}={value!r} is not an integer")
    return int(f)


def _bool(value: str) -> bool:
    if value in ("True", "true", "1"):
        return True
    if value in ("False", "false", "0"):
        return False
    raise _Refused(f"error: assignments.tsv is_rep={value!r} is not a boolean")


def _validated(
    frame: pd.DataFrame, members: dict[str, str], run: dict[str, Any]
) -> dict[str, dict[str, Any]] | tuple[str, str | None]:
    """Check the whole assignments table against this run, all-or-nothing.

    Returns name -> data for an all-OK result, or ``(status, mmseqs_version)`` when the
    worker itself reported a failure for this run. Raises _Refused on anything that
    makes the result not this run's, or not a valid partition of this run's set.
    """
    missing_cols = [c for c in ASSIGNMENT_FIELDS if c not in frame.columns]
    if missing_cols:
        raise _Refused(f"error: assignments.tsv lacks column(s) {', '.join(missing_cols)}")

    tokens = set(frame["token"])
    if tokens != {run["token"]}:
        shown = ", ".join(sorted(t or "(blank)" for t in tokens)) or "(no rows)"
        raise _Refused(
            f"error: stale result -- assignments.tsv token {shown} is not this run's "
            f"{run['token']}; it was written by another submit. The task for this run "
            f"has not written its result yet, or failed before it could"
        )

    names = list(frame["name"])
    dup = [n for n, k in Counter(names).items() if k > 1]
    if dup:
        raise _Refused(f"error: assignments.tsv lists {dup[0]} more than once")
    if set(names) != set(members):
        extra = sorted(set(names) - set(members))
        absent = sorted(set(members) - set(names))
        raise _Refused(
            f"error: assignments.tsv covers {len(names)} design(s), the run submitted "
            f"{len(members)} (absent: {absent[:3]}, unknown: {extra[:3]})"
        )
    for col, want in (("params", run["params"]), ("input_column", run["input_column"])):
        got = set(frame[col])
        if got != {want}:
            raise _Refused(f"error: assignments.tsv {col} {sorted(got)} is not this run's {want!r}")

    statuses = set(frame["status"])
    if statuses != {"OK"}:
        if len(statuses) == 1 and (s := statuses.pop()).startswith("error"):
            versions = sorted(set(frame["mmseqs_version"]) - {""})
            return s, (versions[0] if len(versions) == 1 else None)  # worker's own failure
        raise _Refused(f"error: assignments.tsv carries unexpected status(es) {sorted(statuses)}")

    rows: dict[str, dict[str, Any]] = {}
    n = len(members)
    for rec in frame.to_dict("records"):
        name = rec["name"]
        data = {
            "n_designs_clustered": _int(rec["n_designs_clustered"], "n_designs_clustered"),
            "seq_len": _int(rec["seq_len"], "seq_len"),
            "input_column": rec["input_column"],
            "params": rec["params"],
            "mmseqs_version": rec["mmseqs_version"] or None,
            "cluster_rep": rec["cluster_rep"],
            "cluster_id": _int(rec["cluster_id"], "cluster_id"),
            "cluster_size": _int(rec["cluster_size"], "cluster_size"),
            "is_rep": _bool(rec["is_rep"]),
            "n_clusters": _int(rec["n_clusters"], "n_clusters"),
        }
        if data["seq_len"] != len(members[name]):
            raise _Refused(
                f"error: assignments.tsv seq_len {data['seq_len']} for {name} is not the "
                f"submitted sequence's length {len(members[name])}"
            )
        if data["n_designs_clustered"] != n:
            raise _Refused(
                f"error: assignments.tsv n_designs_clustered {data['n_designs_clustered']} "
                f"is not the {n} designs submitted"
            )
        if data["cluster_rep"] not in members:
            raise _Refused(f"error: cluster_rep {data['cluster_rep']!r} of {name} is not a run member")
        if data["is_rep"] != (data["cluster_rep"] == name):
            raise _Refused(f"error: is_rep of {name} disagrees with its cluster_rep")
        rows[name] = data

    # Set-level consistency: the rows must form a partition of the run's set.
    sizes = Counter(d["cluster_rep"] for d in rows.values())
    reps = {name for name, d in rows.items() if d["is_rep"]}
    if reps != set(sizes):
        raise _Refused("error: the representatives flagged is_rep are not the cluster_rep values")
    id_of: dict[str, int] = {}
    for d in rows.values():
        if id_of.setdefault(d["cluster_rep"], d["cluster_id"]) != d["cluster_id"]:
            raise _Refused(f"error: cluster {d['cluster_rep']} carries two cluster_id values")
    if sorted(id_of.values()) != list(range(1, len(id_of) + 1)):
        raise _Refused("error: cluster_id values are not exactly 1..n_clusters")
    for d in rows.values():
        if d["cluster_size"] != sizes[d["cluster_rep"]]:
            raise _Refused(
                f"error: cluster_size {d['cluster_size']} of cluster {d['cluster_rep']} is "
                f"not its {sizes[d['cluster_rep']]} member rows"
            )
        if d["n_clusters"] != len(sizes):
            raise _Refused(f"error: n_clusters {d['n_clusters']} is not the {len(sizes)} clusters found")
    return rows


def collect_seqcluster(ctx: CollectCtx) -> CollectEach:
    """Set-level seqcluster collector: validate the run's one result once, then hand
    each design its row. Variants are distinguished via --dir-label, matching the
    output dir."""
    # Always re-collect every row (see RERUN SEMANTICS above). This must happen here,
    # BEFORE by_design reads ctx.ready, which filters out OK rows unless args.force.
    ctx.args.force = True

    # Backstop for the builder's stale-label guard (which sees only the -f-filtered
    # frame): rows holding an earlier run's labels that this collect will not visit.
    ready_names = set(map(str, ctx.ready.index))
    if ctx.status_col in ctx.df.columns:
        prior = ctx.df[ctx.status_col]
        stale = [
            str(n)
            for n in ctx.df.index[prior.notna() & (prior != "missing")]
            if str(n) not in ready_names
        ]
        if stale:
            print(
                f"WARNING: seqcluster: {len(stale)} row(s) keep labels from an EARLIER run "
                f"because they have no value in this run's input column, so this collect "
                f"cannot reset them: {', '.join(stale[:5])}{' ...' if len(stale) > 5 else ''}. "
                f"Their cluster_id values belong to a different clustering. Re-run under a "
                f"different -l/--dir-label."
            )

    task_dir = ctx.out_dir / TASK_DIRNAME
    members_path = task_dir / MEMBERS_FILENAME
    run_path = task_dir / RUN_FILENAME
    assignments_path = ctx.out_dir / RESULT_DIRNAME / ASSIGNMENTS

    members: dict[str, str] = {}
    run: dict[str, Any] = {}
    per_row: dict[str, Collected] = {}

    if members_path.is_file() or run_path.is_file():
        # A run set exists. From here on, every member gets OK or the SAME error --
        # never 'missing', which is reserved for rows outside the run's set.
        try:
            if members_path.is_file():
                members = _read_members(members_path)
            if not (members_path.is_file() and run_path.is_file()):
                raise _Refused(
                    f"error: {TASK_DIRNAME}/ holds only one of {MEMBERS_FILENAME} and "
                    f"{RUN_FILENAME}; the submit was interrupted"
                )
            run = json.loads(run_path.read_text())
            if not isinstance(run, dict):
                raise _Refused(f"error: {TASK_DIRNAME}/{RUN_FILENAME} is not a JSON object")
            for key in ("token", "params", "input_column"):
                if not run.get(key):
                    raise _Refused(f"error: {TASK_DIRNAME}/{RUN_FILENAME} has no {key!r}")
            if run.get("n_designs") != len(members):
                raise _Refused(
                    f"error: {RUN_FILENAME} records {run.get('n_designs')} designs but "
                    f"{MEMBERS_FILENAME} lists {len(members)}"
                )
            if not assignments_path.is_file():
                raise _Refused(
                    f"error: no {RESULT_DIRNAME}/{ASSIGNMENTS} -- the clustering task has "
                    f"not finished, or died before writing a result (see the task's .err)"
                )
            frame = pd.read_csv(
                assignments_path, sep="\t", dtype=str, keep_default_na=False
            )
            result = _validated(frame, members, run)
        except _Refused as e:
            result = (str(e), None)
        except (OSError, ValueError) as e:  # unreadable / unparsable file
            result = (f"error: cannot read the run's files: {type(e).__name__}: {e}", None)

        if not members:
            # The run set exists but cannot be read, so its members are unknown. Marking
            # rows 'missing' would claim they were never submitted; raise instead.
            msg = result[0] if isinstance(result, tuple) else "empty members list"
            raise RuntimeError(
                f"seqcluster: cannot collect -- this run's member list is unreadable "
                f"({msg}). Re-run seqcluster under this dir label."
            )
        path = str(assignments_path) if assignments_path.is_file() else ""
        if isinstance(result, tuple):
            status = " ".join(result[0].split())  # one line: it lands in a TSV cell
            for name, seq in members.items():
                data = _empty()
                # Trust/provenance from THIS run's own files, never from a refused result.
                data["seq_len"] = len(seq)
                data["input_column"] = run.get("input_column") or None
                data["params"] = run.get("params") or None
                data["mmseqs_version"] = result[1]
                per_row[name] = Collected(status=status, path=path, data=data)
            print(f"seqcluster: {len(members)} run member(s) -> {status}")
        else:
            for name, data in result.items():
                per_row[name] = Collected(status="OK", path=path, data=data)
            print(f"seqcluster: {len(members)} run member(s) -> OK")

    def one(d: DesignCtx) -> Iterable[Collected]:
        hit = per_row.get(d.name)
        if hit is None:
            # Not in this run's set: explicit NA in EVERY data column (an empty emission
            # would let the driver write only status='missing' and leave old labels).
            yield Collected(status="missing", path="", data=_empty())
            return
        yield hit

    return one
