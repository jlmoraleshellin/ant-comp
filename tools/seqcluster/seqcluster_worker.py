#!/usr/bin/env python3
"""
seqcluster worker: cluster ONE whole set of designs with MMseqs2 ``easy-cluster`` and
write one assignment row per design.

This is the only per-task step of seqcluster, and there is exactly one task per run:
clustering is a SET-LEVEL operation (a sequence's cluster depends on every other
sequence in the set), so the manifest builder hands this worker every ready design at
once. Nothing here is per-design-independent; the task is all-or-nothing.

Input: a members TSV written by the builder, ``name<TAB>sequence`` per line, no header.

Steps:
  1. write a combined FASTA. Each record's header is ``>s<i> <name>``: the first word
     (``s1``, ``s2``, ...) is the ONLY thing MMseqs2 uses as the identifier in its
     cluster TSV, so the mapping back to table names never depends on how MMseqs2
     tokenises a design name (whitespace, ``|``-delimited UniProt-style headers). The
     name is kept as the description so the FASTA stays human-readable.
  2. run ``mmseqs easy-cluster <fasta> <prefix> <tmp> --min-seq-id .. -c .. --cov-mode
     .. --cluster-mode .. --threads ..`` in a node-local temp dir (never on the runs
     Volume / group storage).
  3. parse ``<prefix>_cluster.tsv`` (representative ID, member ID; one line per member,
     the representative listed as its own member) and VALIDATE it: every ID known,
     every design a member exactly once, every representative a member of its own
     cluster. Any violation fails the whole batch -- a partial or mis-mapped clustering
     is never written as OK.
  4. write ``assignments.tsv`` with one row per design.

Output, under ``<out_dir>/seqcluster_result/``:
  assignments.tsv        what collect_seqcluster.py reads (one row per member)
  input.fasta            the exact FASTA clustered
  clusters_cluster.tsv   MMseqs2's raw cluster TSV (s<i> IDs, not names)
  clusters_rep_seq.fasta MMseqs2's representative sequences
  mmseqs.log             MMseqs2 stdout+stderr

On failure every member gets the SAME ``error: ...`` status and the worker exits 1, so
the task's exit code says the run failed AND the table says why. If the worker dies
before writing anything (unreadable members file, a Python crash), there is no
assignments.tsv at all, and the collector marks every member ``error:`` from its own
copy of the members list.

Every row carries the run token passed on the command line; the collector compares it
with the token the builder recorded, so a result left by an EARLIER run is never
collected as this one's.

Stdlib only: the Modal image is debian_slim + the static mmseqs binary.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

RESULT_DIRNAME = "seqcluster_result"
ASSIGNMENTS = "assignments.tsv"

# Columns of assignments.tsv, in order. collect_seqcluster.py reads them by name.
FIELDS = [
    "name",
    "status",
    "token",
    "seq_len",
    "cluster_rep",
    "cluster_id",
    "cluster_size",
    "is_rep",
    "n_designs_clustered",
    "n_clusters",
    "mmseqs_version",
    "params",
    "input_column",
]


class ClusterError(Exception):
    """The batch cannot be assigned: the message becomes every member's status."""


def read_members(path: Path) -> list[tuple[str, str]]:
    """``name<TAB>sequence`` lines -> ordered list. Raises on any malformed line,
    duplicate name or empty sequence (the builder should never write one)."""
    members: list[tuple[str, str]] = []
    seen: set[str] = set()
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ClusterError(f"malformed members line {lineno}")
            name, seq = parts
            if name in seen:
                raise ClusterError(f"duplicate design name in members file: {name}")
            seen.add(name)
            members.append((name, seq))
    if not members:
        raise ClusterError("members file is empty")
    return members


def member_ids(n: int) -> list[str]:
    """The FASTA identifiers MMseqs2 will echo back: s1..sN, in members order."""
    return [f"s{i}" for i in range(1, n + 1)]


def write_fasta(members: list[tuple[str, str]], path: Path) -> list[str]:
    ids = member_ids(len(members))
    with open(path, "w") as f:
        for id_, (name, seq) in zip(ids, members):
            f.write(f">{id_} {name}\n{seq}\n")
    return ids


def parse_cluster_tsv(text: str, ids: list[str]) -> dict[str, str]:
    """MMseqs2 ``_cluster.tsv`` -> {member_id: rep_id}. Raises ClusterError on any
    count or mapping violation (see ``rep_order`` for the representatives' order)."""
    known = set(ids)
    rep_of: dict[str, str] = {}
    for lineno, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            raise ClusterError(f"cluster TSV line {lineno} has {len(parts)} fields, expected 2")
        rep, mem = parts[0].strip(), parts[1].strip()
        for x in (rep, mem):
            if x not in known:
                raise ClusterError(f"cluster TSV names an unknown sequence ID {x!r}")
        if mem in rep_of:
            raise ClusterError(f"sequence {mem} is a member of more than one cluster")
        rep_of[mem] = rep
    missing = [i for i in ids if i not in rep_of]
    if missing:
        raise ClusterError(
            f"{len(missing)} of {len(ids)} designs absent from the cluster TSV "
            f"(first: {missing[0]})"
        )
    for rep in set(rep_of.values()):
        if rep_of.get(rep) != rep:
            raise ClusterError(f"representative {rep} is not a member of its own cluster")
    return rep_of


def rep_order(text: str) -> list[str]:
    """Representative IDs in the order they first appear in the cluster TSV."""
    order: list[str] = []
    seen: set[str] = set()
    for line in text.splitlines():
        if not line.strip():
            continue
        rep = line.split("\t")[0].strip()
        if rep not in seen:
            seen.add(rep)
            order.append(rep)
    return order


def assign(
    members: list[tuple[str, str]], ids: list[str], text: str
) -> list[dict[str, object]]:
    """Validated per-design assignments, keyed back to table names.

    cluster_id is 1-based, by first appearance of the representative in the cluster
    TSV. It is stable only within this one run's output.
    """
    rep_of = parse_cluster_tsv(text, ids)
    name_of = {id_: name for id_, (name, _seq) in zip(ids, members)}
    cid = {rep: i for i, rep in enumerate(rep_order(text), 1)}
    size: dict[str, int] = {}
    for rep in rep_of.values():
        size[rep] = size.get(rep, 0) + 1
    n = len(members)
    n_clusters = len(cid)
    rows: list[dict[str, object]] = []
    for id_, (name, seq) in zip(ids, members):
        rep = rep_of[id_]
        rows.append(
            {
                "name": name,
                "status": "OK",
                "seq_len": len(seq),
                "cluster_rep": name_of[rep],
                "cluster_id": cid[rep],
                "cluster_size": size[rep],
                "is_rep": rep == id_,
                "n_designs_clustered": n,
                "n_clusters": n_clusters,
            }
        )
    if sum(size.values()) != n:  # belt and braces: parse_cluster_tsv already enforces it
        raise ClusterError("cluster sizes do not sum to the number of designs")
    return rows


def write_assignments(path: Path, rows: list[dict[str, object]], common: dict[str, object]) -> None:
    tmp = path.with_suffix(".tsv.partial")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            full = {k: "" for k in FIELDS}
            full.update(common)
            full.update({k: ("" if v is None else v) for k, v in row.items()})
            w.writerow(full)
    tmp.replace(path)  # atomic: a reader never sees half a file


def error_rows(members: list[tuple[str, str]], reason: str) -> list[dict[str, object]]:
    status = reason if reason.startswith("error:") else f"error: {reason}"
    status = " ".join(status.split())  # one line, no tabs: it lands in a TSV cell
    return [{"name": n, "status": status, "seq_len": len(s)} for n, s in members]


def mmseqs_version(mmseqs: str) -> str:
    try:
        p = subprocess.run([mmseqs, "version"], capture_output=True, text=True, timeout=60)
        return p.stdout.strip() or "unknown"
    except Exception:  # noqa: BLE001 -- provenance only; the cluster call is what must succeed
        return "unknown"


def n_threads() -> int:
    for var in ("SEQCLUSTER_THREADS", "SLURM_CPUS_PER_TASK"):
        v = os.environ.get(var, "").strip()
        if v.isdigit() and int(v) > 0:
            return int(v)
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def cluster_args(min_seq_id: str, coverage: str, cov_mode: str, cluster_mode: str) -> list[str]:
    return [
        "--min-seq-id", min_seq_id,
        "-c", coverage,
        "--cov-mode", cov_mode,
        "--cluster-mode", cluster_mode,
    ]


def run(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    result_dir = out_dir / RESULT_DIRNAME
    result_dir.mkdir(parents=True, exist_ok=True)
    assignments = result_dir / ASSIGNMENTS

    # MUST match run_seqcluster.params_string() for the same four manifest fields: the
    # collector refuses a result whose params differ from the run's (test-asserted).
    params = (
        f"min_seq_id={args.min_seq_id};c={args.coverage};"
        f"cov_mode={args.cov_mode};cluster_mode={args.cluster_mode}"
    )
    common: dict[str, object] = {
        "token": args.token,
        "params": params,
        "input_column": args.input_column,
    }

    try:
        members = read_members(Path(args.members))
    except (OSError, ClusterError) as e:
        # Without a readable members list there is nobody to write an error row for.
        # No assignments.tsv is written; the collector turns its absence into an
        # 'error:' on every member of the run (it reads the members list itself).
        print(f"seqcluster: cannot read members file: {e}", file=sys.stderr)
        return 1

    mmseqs = args.mmseqs
    common["mmseqs_version"] = mmseqs_version(mmseqs)
    print(f"seqcluster: {len(members)} designs, {params}, mmseqs {common['mmseqs_version']}")

    work = Path(tempfile.mkdtemp(prefix="seqcluster_"))
    try:
        fasta = work / "input.fasta"
        ids = write_fasta(members, fasta)
        prefix = work / "clusters"
        cmd = [
            mmseqs, "easy-cluster", str(fasta), str(prefix), str(work / "tmp"),
            *cluster_args(args.min_seq_id, args.coverage, args.cov_mode, args.cluster_mode),
            "--threads", str(n_threads()),
        ]
        print("seqcluster: " + " ".join(cmd))
        log = result_dir / "mmseqs.log"
        try:
            with open(log, "w") as lf:
                proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
            rc = proc.returncode
        except OSError as e:  # binary missing / not executable
            raise ClusterError(f"could not execute mmseqs ({mmseqs}): {e}") from e
        shutil.copy2(fasta, result_dir / "input.fasta")
        if rc != 0:
            tail = log.read_text(errors="replace").strip().splitlines()[-5:]
            print("\n".join(tail), file=sys.stderr)
            raise ClusterError(f"mmseqs easy-cluster exited {rc} (see seqcluster_result/mmseqs.log)")
        tsv = Path(f"{prefix}_cluster.tsv")
        if not tsv.is_file():
            raise ClusterError("mmseqs easy-cluster wrote no _cluster.tsv")
        shutil.copy2(tsv, result_dir / "clusters_cluster.tsv")
        rep_fa = Path(f"{prefix}_rep_seq.fasta")
        if rep_fa.is_file():
            shutil.copy2(rep_fa, result_dir / "clusters_rep_seq.fasta")
        rows = assign(members, ids, tsv.read_text())
    except ClusterError as e:
        print(f"seqcluster: {e}", file=sys.stderr)
        write_assignments(assignments, error_rows(members, str(e)), common)
        return 1
    finally:
        shutil.rmtree(work, ignore_errors=True)

    write_assignments(assignments, rows, common)
    n_clusters = rows[0]["n_clusters"] if rows else 0
    print(f"seqcluster: {len(rows)} designs -> {n_clusters} clusters")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--members", required=True)
    p.add_argument("--token", required=True)
    p.add_argument("--min-seq-id", required=True)
    p.add_argument("--coverage", required=True)
    p.add_argument("--cov-mode", required=True)
    p.add_argument("--cluster-mode", required=True)
    p.add_argument("--input-column", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--mmseqs", default=os.environ.get("MMSEQS", "mmseqs"))
    return run(p.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
