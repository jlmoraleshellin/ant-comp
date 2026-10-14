#!/usr/bin/env python3
"""
Self-test for seqcluster. MMseqs2 is replaced by a STUB executable (a few lines of
Python that clusters by a sequence prefix) for everything except the optional
real-binary section at the end.

WHAT THIS PROVES, AND WHAT IT CANNOT.

It proves the plumbing and the contracts: that the builder emits exactly ONE task for
the whole set and ignores the already-OK resume filter, refuses a missing column, a
multi-chain ('/' or ':') sequence, foreign characters and an input-column change that
would strand an earlier run's labels; that the worker maps MMseqs2's s<i> identifiers
back to table names (names with spaces and '|' included), numbers clusters by first
appearance in the cluster TSV and not by size, and turns every malformed or failed
MMseqs2 result into the SAME error on every member; that the collector forces a full
re-collect even without --force, writes 'missing' with EVERY data column NA for rows
outside the run set, refuses a result whose token / params / membership is not this
run's as an 'error:' (never 'missing'), and re-checks that the result is a partition;
that the .sh cuts its manifest fields and clears a stale result; and -- through the
REAL prosapia run/collect drivers with only the executor faked -- that the columns land
in the TSV with dtypes the skill's filter expressions work on.

It CANNOT prove the premise or the real-data mapping. The stub agrees with this file
by construction, and the optional real-binary check uses synthetic sequences. UNKNOWN
until a worker runs it on real data:

  * whether the intended input column (e.g. proteinmpnn_redesign_full_sequence) holds
    ONE chain per row in the real table, or a '/'-joined complex (which is refused);
  * whether the real column is upper-case one-letter code with no gap/stop characters;
  * whether the pinned static sse4.1 binary runs on the Modal host the task lands on
    (the image build smoke-tests it on the BUILD host only);
  * whether MMseqs2 with its default parameters groups real ProteinMPNN redesigns the
    way the caller expects (min-seq-id 0.0 means coverage and E-value decide alone);
  * whether the Modal Volume shows the builder's unlink of the old assignments.tsv to
    the task container (the run token makes this a non-issue for correctness).

Run:  uv run python tools/seqcluster/test_seqcluster.py
      SEQCLUSTER_TEST_MMSEQS=/path/to/mmseqs uv run python tools/seqcluster/test_seqcluster.py
      (the second form adds the real-binary checks; without it they are SKIPPED)
"""

import contextlib
import importlib.util
import io
import json
import os
import random
import subprocess
import sys
import tempfile
from argparse import ArgumentParser, Namespace
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


worker = _load("seqcluster_worker")
run_mod = _load("run_seqcluster")
collect_mod = _load("collect_seqcluster")

from prosapia.core import CollectCtx, ManifestCtx  # noqa: E402
from prosapia.core import base_run  # noqa: E402
from prosapia.core.base_collect import (  # noqa: E402
    build_collect_parser,
    by_design,
    collect_from_args,
)
from prosapia.core.base_run import build_run_parser, run_from_args  # noqa: E402
from prosapia.core.data_manager import DataManager, Table  # noqa: E402
from prosapia.core.naming import RUN_META_FILENAME  # noqa: E402
from prosapia.core.tool_registry import _load as load_tool  # noqa: E402

FAILURES: list[str] = []
CHECKS = 0
SKIPPED: list[str] = []
LEAF = "seqcluster"
SCIENCE = collect_mod.SCIENCE_COLUMNS


def check(cond: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not cond:
        FAILURES.append(label)


def raises(fn, *needles: str) -> bool:
    try:
        fn()
    except Exception as e:  # noqa: BLE001
        msg = str(e)
        missing = [n for n in needles if n not in msg]
        if missing:
            FAILURES.append(f"error raised but lacks {missing}: {msg[:300]}")
        return True
    return False


def quiet(fn, *a, **kw):
    """Call fn, swallowing stdout; return (result, stdout)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **kw)
    return out, buf.getvalue()


# ---------------------------------------------------------------------------
# The MMseqs2 stub
# ---------------------------------------------------------------------------

STUB = r'''#!{python}
"""Stand-in for `mmseqs`. Clusters by the first STUB_K residues (default 8). Within a
family the representative is the LAST member; families are written to the cluster
TSV in REVERSE order of first appearance -- so neither the FASTA order nor cluster
size can be mistaken for the numbering. STUB_MODE corrupts the output on purpose."""
import os, sys
args = sys.argv[1:]
if args[:1] == ["version"]:
    print("stub-0000"); sys.exit(0)
assert args[0] == "easy-cluster", args
fasta, prefix, tmp = args[1:4]
if os.environ.get("STUB_ARGLOG"):
    open(os.environ["STUB_ARGLOG"], "w").write("\t".join(args))
mode = os.environ.get("STUB_MODE", "ok")
if mode == "fail":
    print("stub: boom", file=sys.stderr); sys.exit(7)
if mode == "notsv":
    sys.exit(0)
ids, seqs = [], {}
cur = None
for line in open(fasta):
    line = line.rstrip("\n")
    if line.startswith(">"):
        cur = line[1:].split()[0]; ids.append(cur); seqs[cur] = ""
    else:
        seqs[cur] += line
k = int(os.environ.get("STUB_K", "8"))
fam = {}
for i in ids:
    fam.setdefault(seqs[i][:k], []).append(i)
lines = []
for members in reversed(list(fam.values())):
    rep = members[-1]
    lines += [f"{rep}\t{rep}"] + [f"{rep}\t{m}" for m in members if m != rep]
if mode == "drop":
    lines = lines[:-1]
elif mode == "unknown":
    lines.append("s999\ts999")
elif mode == "twice":
    first = lines[0].split("\t")[0]
    lines.append(f"{first}\t{lines[-1].split(chr(9))[1]}")
elif mode == "threefields":
    lines[0] += "\textra"
open(prefix + "_cluster.tsv", "w").write("\n".join(lines) + "\n")
open(prefix + "_rep_seq.fasta", "w").write("".join(f">{l.split()[0]}\n" for l in lines))
'''


def make_stub(tmp: Path) -> Path:
    p = tmp / "mmseqs_stub"
    p.write_text(STUB.replace("{python}", sys.executable))
    p.chmod(0o755)
    return p


def family(prefix: str, n: int) -> list[str]:
    """n sequences sharing an 8-residue prefix (one stub family), distinct tails."""
    return [prefix + "ACDEFGHIK"[: (j % 9) + 1] + "W" * j for j in range(n)]


# ---------------------------------------------------------------------------
# Helpers around the builder, worker and collector
# ---------------------------------------------------------------------------


def manifest_args(**kw) -> Namespace:
    return Namespace(
        run_dir=kw.pop("run_dir", Path("/nonexistent")), table="table1",
        input_column=kw.pop("input_column", "seq"), force=kw.pop("force", False),
        gpus_per_task=1, min_seq_id=kw.pop("min_seq_id", 0.0),
        coverage=kw.pop("coverage", 0.8), cov_mode=kw.pop("cov_mode", 0),
        cluster_mode=kw.pop("cluster_mode", 0),
    )


def build(df: pd.DataFrame, out_dir: Path, **kw):
    out_dir.mkdir(parents=True, exist_ok=True)
    args = manifest_args(**kw)
    ctx = ManifestCtx(df=df, args=args, out_dir=out_dir,
                      lookup=lambda name, column: None, write_meta=lambda **_: None)
    rows, printed = quiet(run_mod.build_seqcluster_manifest, ctx)
    return rows, args, printed


def run_worker(row: tuple, out_dir: Path, mmseqs: Path | str, mode: str = "ok",
               arglog: Path | None = None) -> int:
    env_backup = dict(os.environ)
    os.environ["STUB_MODE"] = mode
    if arglog:
        os.environ["STUB_ARGLOG"] = str(arglog)
    try:
        rc, _ = quiet(worker.main, [
            "--members", row[0], "--token", row[1], "--min-seq-id", row[2],
            "--coverage", row[3], "--cov-mode", row[4], "--cluster-mode", row[5],
            "--input-column", row[6], "--out-dir", str(out_dir), "--mmseqs", str(mmseqs),
        ])
    finally:
        os.environ.clear()
        os.environ.update(env_backup)
    return rc


def assignments(out_dir: Path) -> pd.DataFrame:
    return pd.read_csv(out_dir / "seqcluster_result" / "assignments.tsv", sep="\t",
                       dtype=str, keep_default_na=False)


def collect(df: pd.DataFrame, out_dir: Path, input_column: str = "seq", force: bool = False):
    (out_dir / RUN_META_FILENAME).write_text(json.dumps({"input_column": input_column}))
    args = Namespace(run_dir=out_dir.parent, table="table1", dir_label="", force=force)
    ctx = CollectCtx(
        df=df, args=args, table_name="table1", out_dir=out_dir,
        status_col=f"{LEAF}_status", path_col=f"{LEAF}_path", parent_table=None,
        parent_df=pd.DataFrame(), lookup=lambda n, c: None, creates_table=False,
        default_input_column="not applicable",
    )
    updates, printed = quiet(by_design(collect_mod.collect_seqcluster), ctx)
    return updates, printed, args


def table(seqs: dict[str, object], **extra) -> pd.DataFrame:
    df = pd.DataFrame({"seq": pd.Series(seqs, dtype=object)})
    for col, values in extra.items():
        df[col] = pd.Series(values, dtype=object)
    df.index.name = "name"
    return df


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_layout_and_defaults() -> None:
    for const in ("TASK_DIRNAME", "MEMBERS_FILENAME", "RUN_FILENAME", "RESULT_DIRNAME", "ASSIGNMENTS"):
        check(getattr(run_mod, const) == getattr(collect_mod, const), f"run/collect agree on {const}")
    for const in ("RESULT_DIRNAME", "ASSIGNMENTS"):
        check(getattr(run_mod, const) == getattr(worker, const), f"run/worker agree on {const}")
    check(worker.FIELDS == collect_mod.ASSIGNMENT_FIELDS,
          f"worker FIELDS == collector ASSIGNMENT_FIELDS ({worker.FIELDS})")
    sh = (HERE / "seqcluster.sh").read_text()
    check(f"{run_mod.RESULT_DIRNAME}/{run_mod.ASSIGNMENTS}" in sh, ".sh clears the same assignments path")
    check(set(collect_mod.DATA_COLUMNS) == set(worker.FIELDS) - {"name", "status", "token"},
          "collector DATA_COLUMNS are the worker's data fields")

    p = ArgumentParser()
    run_mod.add_run_seqcluster_args(p)
    a = p.parse_args([])
    check((a.min_seq_id, a.coverage, a.cov_mode, a.cluster_mode) == (0.0, 0.8, 0, 0),
          f"flag defaults are MMseqs2 easy-cluster's own ({a})")
    check(run_mod.params_string(*run_mod.manifest_params(a))
          == "min_seq_id=0.0;c=0.8;cov_mode=0;cluster_mode=0", "default params string")

    tool = load_tool(HERE / "spec.py")
    check(tool.name == "seqcluster" and tool.action == "update"
          and tool.default_input_column == run_mod.NO_DEFAULT_COLUMN,
          f"spec: name/action/default column ({tool.name}, {tool.action})")


def test_builder(tmp: Path) -> None:
    seqs = {"des 1": "mkvlaa", "a|b": "MKVLAG", "c": "WWWWWW", "nan": None, "blank": "  "}
    rows, args, printed = build(table(seqs), tmp / "b1" / LEAF)
    check(len(rows) == 1, f"ONE manifest row for the whole set ({len(rows)})")
    row = rows[0]
    check(len(row) == 7 and all(isinstance(f, str) and f for f in row), f"7 non-empty fields ({row})")
    check(args.gpus_per_task == 0, "builder forces gpus_per_task = 0")
    members = [l.split("\t") for l in Path(row[0]).read_text().splitlines()]
    check(members == [["des 1", "MKVLAA"], ["a|b", "MKVLAG"], ["c", "WWWWWW"]],
          f"members: present rows only, in table order, upper-cased ({members})")
    check(row[2:6] == ("0.0", "0.8", "0", "0") and row[6] == "seq", f"param/column fields ({row[2:]})")
    run = json.loads((tmp / "b1" / LEAF / "seqcluster_task" / "run.json").read_text())
    check(run == {"token": row[1], "params": "min_seq_id=0.0;c=0.8;cov_mode=0;cluster_mode=0",
                  "input_column": "seq", "n_designs": 3}, f"run.json ({run})")
    check("2 row(s) have no seq" in printed, "absent sequences reported")
    rows2, _, _ = build(table(seqs), tmp / "b1" / LEAF)
    check(rows2[0][1] != row[1], "every submit gets a fresh token")

    # ctx.ready would drop these OK rows; the builder must not.
    df = table({"x": "MKV", "y": "MKW", "z": "MKY"}, seqcluster_status={"x": "OK", "y": "OK", "z": None})
    rows, _, printed = build(df, tmp / "b2" / LEAF)
    n = len(Path(rows[0][0]).read_text().splitlines())
    check(n == 3, f"already-OK rows are RE-CLUSTERED, not skipped ({n} of 3)")
    check("RE-CLUSTERED" in printed, "re-cluster notice printed when earlier OK rows exist")

    # The builder drops the previous result so it cannot be paired with new members.
    res = tmp / "b2" / LEAF / "seqcluster_result"
    res.mkdir(parents=True, exist_ok=True)
    (res / "assignments.tsv").write_text("stale\n")
    build(df, tmp / "b2" / LEAF)
    check(not (res / "assignments.tsv").exists(), "previous assignments.tsv removed at submit")

    good = table({"x": "MKV"}, proteinmpnn_sequence={"x": "MKV"})
    check(raises(lambda: build(good, tmp / "b3" / LEAF, input_column="not applicable"),
                 "no default input column", "proteinmpnn_sequence"),
          "no -i: raises, listing sequence-like columns")
    check(raises(lambda: build(good, tmp / "b3" / LEAF, input_column="boltz_seq"), "'boltz_seq' is not in table"),
          "absent column raises (no lineage walk)")

    mc = table({"m1": "MKV/GGG", "m2": "MKV:GGG", "ok": "MKV"})
    check(raises(lambda: build(mc, tmp / "b4" / LEAF), "MULTI-CHAIN", "m1", "m2", "ProteinMPNN", "'/'"),
          "multi-chain sequences refused, naming the designs and ProteinMPNN's '/' join")
    check(not (tmp / "b4" / LEAF / "seqcluster_task" / "members.tsv").exists(), "refusal writes no members file")
    check(raises(lambda: build(table({"s": "MKV*"}), tmp / "b5" / LEAF), "s ('*')"), "stop character refused")
    check(raises(lambda: build(table({"s": "MKV"}), tmp / "b5" / LEAF, min_seq_id=1.5), "--min-seq-id"),
          "out-of-range --min-seq-id refused")
    check(raises(lambda: build(table({"s": "MKV"}), tmp / "b5" / LEAF, coverage=-0.1), "--coverage"),
          "out-of-range --coverage refused")

    # Stale-label guard: an earlier run labelled 'old', which has no value in THIS column.
    stale = table({"new": "MKV", "old": None}, seqcluster_status={"new": "OK", "old": "OK"})
    check(raises(lambda: build(stale, tmp / "b6" / LEAF), "old", "-l/--dir-label"),
          "row with earlier labels but no input -> refused (labels would survive)")
    fine = table({"new": "MKV", "old": None}, seqcluster_status={"new": "OK", "old": "missing"})
    check(not raises(lambda: build(fine, tmp / "b6" / LEAF)), "a 'missing' row with no input is no obstacle")

    # Nothing to submit: earlier task files are cleared, so collect cannot pair them.
    rows, _, _ = build(table({"s": "MKV"}), tmp / "b7" / LEAF)
    rows, _, _ = build(table({"s": None}), tmp / "b7" / LEAF)
    check(rows == [], "no present sequence -> no task")
    check(not (tmp / "b7" / LEAF / "seqcluster_task" / "members.tsv").exists()
          and not (tmp / "b7" / LEAF / "seqcluster_task" / "run.json").exists(),
          "empty submit clears the earlier members.tsv / run.json")


def test_worker(tmp: Path) -> None:
    stub = make_stub(tmp)
    # Families of sizes 4, 1, 2, in that table order. The stub writes them to the TSV
    # in REVERSE order (2, 1, 4), so by first appearance the pair is cluster 1, the solo
    # 2 and the quad 3 -- different from table order (quad first) AND from size order
    # (quad first), so neither can pass for the real numbering.
    names = ["q1", "q 2", "q|3", "q4", "solo", "pair a", "pair|b"]
    seqs = family("DDDDDDDD", 4) + family("AAAAAAAA", 1) + family("CCCCCCCC", 2)
    df = table(dict(zip(names, seqs)))
    rows, _, _ = build(df, tmp / "w1" / LEAF)
    log = tmp / "arglog"
    rc = run_worker(rows[0], tmp / "w1" / LEAF, stub, arglog=log)
    check(rc == 0, f"worker exits 0 on a good clustering ({rc})")
    a = assignments(tmp / "w1" / LEAF).set_index("name")
    check(list(a.index) == names, f"one row per member, table order ({list(a.index)})")
    check(set(a["status"]) == {"OK"} and set(a["token"]) == {rows[0][1]}, "all OK, token echoed")
    check(a.loc["q1", "cluster_rep"] == "q4" and a.loc["pair a", "cluster_rep"] == "pair|b"
          and a.loc["solo", "cluster_rep"] == "solo",
          f"s<i> IDs mapped back to names incl. spaces and '|' ({a['cluster_rep'].to_dict()})")
    check(a.loc["q1", "cluster_id"] == "3" and a.loc["pair a", "cluster_id"] == "1"
          and a.loc["solo", "cluster_id"] == "2",
          f"cluster_id by first appearance in the TSV, not size/table order ({a['cluster_id'].to_dict()})")
    check(a.loc["q 2", "cluster_size"] == "4" and a.loc["pair|b", "cluster_size"] == "2"
          and a.loc["solo", "cluster_size"] == "1", "cluster sizes")
    check(sorted(a.index[a["is_rep"] == "True"]) == ["pair|b", "q4", "solo"], "one is_rep per cluster")
    check(set(a["n_clusters"]) == {"3"} and set(a["n_designs_clustered"]) == {"7"}, "set-level counts")
    check(a.loc["q4", "seq_len"] == str(len(seqs[3])), "seq_len is the clustered length")
    check(set(a["params"]) == {run_mod.params_string(*rows[0][2:6])},
          "worker params string == builder params_string (collector requires equality)")
    check(set(a["mmseqs_version"]) == {"stub-0000"} and set(a["input_column"]) == {"seq"}, "provenance")
    fasta = (tmp / "w1" / LEAF / "seqcluster_result" / "input.fasta").read_text().splitlines()
    check(fasta[0] == ">s1 q1" and fasta[2] == ">s2 q 2",
          f"FASTA identifiers are s<i>, the name only a description ({fasta[:3]})")
    argv = log.read_text().split("\t")
    check(argv[4:12] == ["--min-seq-id", "0.0", "-c", "0.8", "--cov-mode", "0", "--cluster-mode", "0"]
          and argv[12] == "--threads", f"easy-cluster called with explicit params ({argv[4:]})")
    check((tmp / "w1" / LEAF / "seqcluster_result" / "clusters_cluster.tsv").is_file(), "raw cluster TSV kept")

    rows, _, _ = build(table({"only": "MKVLAAGG"}), tmp / "w2" / LEAF)
    check(run_worker(rows[0], tmp / "w2" / LEAF, stub) == 0, "single sequence clusters")
    a = assignments(tmp / "w2" / LEAF).iloc[0]
    check(a["cluster_rep"] == "only" and a["is_rep"] == "True" and a["cluster_size"] == "1"
          and a["n_clusters"] == "1", "a lone sequence is its own cluster of 1")

    rows, _, _ = build(df, tmp / "w3" / LEAF)
    for mode, needle in [("fail", "exited 7"), ("notsv", "no _cluster.tsv"), ("drop", "absent"),
                         ("unknown", "unknown sequence ID"), ("twice", "more than one cluster"),
                         ("threefields", "3 fields")]:
        rc = run_worker(rows[0], tmp / "w3" / LEAF, stub, mode=mode)
        a = assignments(tmp / "w3" / LEAF)
        st = set(a["status"])
        check(rc == 1 and len(st) == 1 and next(iter(st)).startswith("error:") and needle in next(iter(st)),
              f"{mode}: exit 1, ONE shared error naming the cause ({rc}, {st})")
        check(all((a[c] == "").all() for c in SCIENCE) and (a["seq_len"] != "").all()
              and set(a["token"]) == {rows[0][1]},
              f"{mode}: science blank, seq_len and token kept")
    rc = run_worker(rows[0], tmp / "w3" / LEAF, tmp / "no_such_mmseqs")
    check(rc == 1 and "could not execute mmseqs" in assignments(tmp / "w3" / LEAF)["status"][0],
          "missing binary -> error rows")
    (tmp / "w3" / LEAF / "seqcluster_result" / "assignments.tsv").unlink()
    bad = list(rows[0])
    bad[0] = str(tmp / "nope.tsv")
    check(run_worker(tuple(bad), tmp / "w3" / LEAF, stub) == 1
          and not (tmp / "w3" / LEAF / "seqcluster_result" / "assignments.tsv").exists(),
          "unreadable members file -> exit 1, no result (the collector reports it)")


def _ready_run(tmp: Path, sub: str, df: pd.DataFrame, mode: str = "ok", build_df=None):
    out = tmp / sub / LEAF
    rows, _, _ = build(df if build_df is None else build_df, out)
    run_worker(rows[0], out, make_stub(tmp), mode=mode)
    return out, rows[0]


def test_collector(tmp: Path) -> None:
    names = ["f1", "f2", "g1", "g2", "g3"]
    seqs = family("EEEEEEEE", 2) + family("FFFFFFFF", 3)
    df = table(dict(zip(names, seqs)))
    out, row = _ready_run(tmp, "c1", df)
    path = str(out / "seqcluster_result" / "assignments.tsv")

    up, _, args = collect(df, out)
    check(set(up) == set(names), f"every member collected ({sorted(up)})")
    check(all(up[n][f"{LEAF}_status"] == "OK" and up[n][f"{LEAF}_path"] == path for n in names), "OK + path")
    g = up["g2"]
    check(g[f"{LEAF}_cluster_rep"] == "g3" and g[f"{LEAF}_cluster_id"] == 1 and g[f"{LEAF}_cluster_size"] == 3
          and g[f"{LEAF}_is_rep"] is False and up["g3"][f"{LEAF}_is_rep"] is True
          and g[f"{LEAF}_n_clusters"] == 2 and g[f"{LEAF}_n_designs_clustered"] == 5,
          f"science columns typed and correct ({g})")
    check(set(g) == {f"{LEAF}_{c}" for c in collect_mod.DATA_COLUMNS} | {f"{LEAF}_status", f"{LEAF}_path"},
          "exact collected column set")
    check(args.force is True, "collector forced args.force = True")

    # The ctx.ready override: all rows already OK and no --force -> still re-collected.
    ok_df = df.copy()
    ok_df[f"{LEAF}_status"] = "OK"
    up, _, _ = collect(ok_df, out, force=False)
    check(set(up) == set(names), f"already-OK rows re-collected without --force ({len(up)} of 5)")

    # A row filtered out of the run (-f) but present in the table -> missing, ALL NA.
    out2, _ = _ready_run(tmp, "c2", df, build_df=df.drop(index=["g1"]))
    up, _, _ = collect(df, out2)
    m = up["g1"]
    check(m[f"{LEAF}_status"] == "missing" and all(m[f"{LEAF}_{c}"] is None for c in collect_mod.DATA_COLUMNS),
          f"row outside the run set -> 'missing' with every data column explicitly NA ({m})")
    check(up["g2"][f"{LEAF}_cluster_size"] == 2, "the run's own set (4 designs) is what was clustered")

    # A row with no sequence is never visited; if it holds earlier labels, say so loudly.
    nseq = df.copy()
    nseq.loc["g1", "seq"] = None
    nseq[f"{LEAF}_status"] = pd.Series({"g1": "OK"}, dtype=object)
    up, printed, _ = collect(nseq, out2)
    check("g1" not in up and "WARNING" in printed and "g1" in printed,
          "stale labels on a row the collect cannot visit -> loud warning")

    def all_error(up, needle):
        st = {up[n][f"{LEAF}_status"] for n in names}
        s = next(iter(st))
        return (len(st) == 1 and s.startswith("error:") and needle in s
                and all(up[n][f"{LEAF}_{c}"] is None for n in names for c in SCIENCE)
                and all(up[n][f"{LEAF}_seq_len"] == len(seqs[i]) for i, n in enumerate(names))
                and up["f1"][f"{LEAF}_params"] == "min_seq_id=0.0;c=0.8;cov_mode=0;cluster_mode=0")

    # Token mismatch: the result was written for ANOTHER submit -> error, never missing.
    run_json = out / "seqcluster_task" / "run.json"
    good_run = run_json.read_text()
    run_json.write_text(good_run.replace(row[1], "f" * 32))
    up, _, _ = collect(df, out)
    check(all_error(up, "stale result"), f"token mismatch -> every member error ({up['f1']})")
    run_json.write_text(good_run)

    # No result yet (collect before the task finished, or the task died silently).
    a_path = out / "seqcluster_result" / "assignments.tsv"
    good_a = a_path.read_text()
    a_path.unlink()
    up, _, _ = collect(df, out)
    check(all_error(up, "no seqcluster_result/assignments.tsv"), "absent result -> every member error")
    check(all(up[n][f"{LEAF}_path"] == "" for n in names), "no path when there is no result file")
    a_path.write_text(good_a)

    # The worker's own failure passes through, provenance kept.
    out3, _ = _ready_run(tmp, "c3", df, mode="fail")
    up, _, _ = collect(df, out3)
    check(all_error(up, "exited 7") and up["g1"][f"{LEAF}_mmseqs_version"] == "stub-0000",
          "worker error status propagated to every member, mmseqs_version kept")

    # Tampered results: each must be refused as a whole.
    frame = pd.read_csv(a_path, sep="\t", dtype=str, keep_default_na=False)

    def tampered(edit, needle, label):
        f = frame.copy()
        edit(f)
        f.to_csv(a_path, sep="\t", index=False)
        up, _, _ = collect(df, out)
        check(all_error(up, needle), f"{label} -> every member error ({up['f1'][f'{LEAF}_status']})")

    tampered(lambda f: f.__setitem__("cluster_size", f["cluster_size"].where(f.name != "g1", "2")),
             "cluster_size", "cluster_size disagreeing with member rows")
    tampered(lambda f: f.__setitem__("is_rep", f["is_rep"].where(f.name != "g1", "True")),
             "is_rep", "is_rep disagreeing with cluster_rep")
    tampered(lambda f: f.__setitem__("params", "min_seq_id=0.5;c=0.8;cov_mode=0;cluster_mode=0"),
             "params", "params not this run's")
    tampered(lambda f: f.__setitem__("seq_len", f["seq_len"].where(f.name != "f1", "999")),
             "seq_len", "seq_len not the submitted length")
    tampered(lambda f: f.__setitem__("name", f["name"].where(f.name != "f1", "intruder")),
             "covers", "membership not the run's")
    tampered(lambda f: f.__setitem__("status", f["status"].where(f.name != "f1", "error: x")),
             "unexpected status", "mixed OK/error statuses")
    tampered(lambda f: f.__setitem__("cluster_id", f["cluster_id"].replace({"2": "5"})),
             "1..n_clusters", "cluster_id not 1..n")
    tampered(lambda f: f.drop(f.index, inplace=True), "stale result", "header-only result")
    a_path.write_text(good_a)
    up, _, _ = collect(df, out)
    check(all(up[n][f"{LEAF}_status"] == "OK" for n in names), "restored result collects OK again")

    # No run under this label at all -> every ready row missing, all NA.
    empty_out = tmp / "c4" / LEAF
    empty_out.mkdir(parents=True)
    up, _, _ = collect(df, empty_out)
    check(all(up[n][f"{LEAF}_status"] == "missing" for n in names), "no run set -> all missing")

    # Interrupted submits.
    (out / "seqcluster_task" / "run.json").unlink()
    up, _, _ = collect(df, out)
    check(all_error_no_params(up, names, "interrupted"), "members without run.json -> every member error")
    (out / "seqcluster_task" / "run.json").write_text(good_run)
    (out / "seqcluster_task" / "members.tsv").unlink()
    check(raises(lambda: collect(df, out), "member list is unreadable"),
          "run.json without members -> collect RAISES (members unknown, never 'missing')")


def all_error_no_params(up, names, needle) -> bool:
    st = {up[n][f"{LEAF}_status"] for n in names}
    return len(st) == 1 and needle in next(iter(st)) and next(iter(st)).startswith("error:")


PRELUDE = 'OUT_DIR="$SH_OUT"; SAPIA_TASK_ID=1; SAPIA_LINE="$SH_LINE"\nsapia_activate() { :; }\n'


def sh_env(tmp: Path, out: Path, line: str, **kw) -> dict:
    prelude = tmp / "prelude.sh"
    prelude.write_text(PRELUDE)
    return dict(os.environ, SAPIA_PRELUDE=str(prelude), SAPIA_TOOL_DIR=str(HERE),
                PIPELINE_PYTHON=kw.get("py", sys.executable), MMSEQS=str(kw.get("mmseqs", "")),
                SH_OUT=str(out), SH_LINE=line, STUB_MODE="ok")


def test_task_script(tmp: Path) -> None:
    stub = make_stub(tmp)
    df = table(dict(zip(["a", "b", "c"], family("GGGGGGGG", 2) + family("HHHHHHHH", 1))))
    out = tmp / "sh" / LEAF
    rows, _, _ = build(df, out)
    res = out / "seqcluster_result"
    res.mkdir(parents=True, exist_ok=True)
    (res / "assignments.tsv").write_text("name\tstatus\ttoken\nb\tOK\tSTALE\n")
    env = sh_env(tmp, out, "\t".join(rows[0]), mmseqs=stub)
    p = subprocess.run(["bash", str(HERE / "seqcluster.sh")], env=env, capture_output=True, text=True)
    check(p.returncode == 0, f"task script exits 0 ({p.returncode}: {p.stderr[-300:]})")
    a = assignments(out)
    check(set(a["token"]) == {rows[0][1]} and len(a) == 3 and set(a["status"]) == {"OK"},
          "task script: fields cut correctly, fresh result over a stale one")

    (res / "assignments.tsv").write_text("stale\n")
    fake = tmp / "fakepy"
    fake.write_text("#!/bin/bash\nexit 3\n")
    fake.chmod(0o755)
    env = sh_env(tmp, out, "\t".join(rows[0]), mmseqs=stub, py=str(fake))
    p = subprocess.run(["bash", str(HERE / "seqcluster.sh")], env=env, capture_output=True, text=True)
    check(p.returncode == 3, f"crashed worker -> task exits with the worker's code ({p.returncode})")
    check(not (res / "assignments.tsv").exists(), "crashed worker leaves NO result, not a stale one")
    up, _, _ = collect(df, out)
    check(all(up[n][f"{LEAF}_status"].startswith("error:") for n in "abc"),
          "...which the collector turns into an error on every member")


def test_drivers_end_to_end(tmp: Path) -> None:
    """The REAL run/collect drivers; only the executor is faked (it runs the .sh)."""
    tool = load_tool(HERE / "spec.py")
    stub = make_stub(tmp)
    run_dir = tmp / "run"
    run_dir.mkdir()
    names = ["d1", "d2", "d3", "e1", "e2", "x"]
    seqs = family("KKKKKKKK", 3) + family("LLLLLLLL", 2) + ["MMMMMMMMQQ"]
    df = table(dict(zip(names, seqs)))
    df.loc["nos"] = [None]
    with DataManager(run_dir) as (dm, (_read, save), registry):
        quiet(registry.register_table, Table("table0", gen=0))
        save("table0", df)

    submits: list = []

    def fake_submit(sctx):
        for r in sctx.rows:
            env = sh_env(tmp, sctx.out_dir, "\t".join(r), mmseqs=stub)
            submits.append(subprocess.run(["bash", str(sctx.script)], env=env,
                                          capture_output=True, text=True))

    real_get = base_run.get_executor
    base_run.get_executor = lambda name: fake_submit
    try:
        rp = ArgumentParser(parents=[build_run_parser(tool.metadata, tool.default_script,
                                                      tool.default_input_column, tool.add_run_args_fn)])
        cp = ArgumentParser(parents=[build_collect_parser(tool.add_collect_args_fn)])

        def sapia_run(*extra):
            a = rp.parse_args([str(run_dir), "-t", "table0", *extra])
            return quiet(run_from_args, tool.metadata, tool.build_manifest_fn, a)[1]

        def sapia_collect():
            a = cp.parse_args([str(run_dir), "-t", "table0"])
            quiet(collect_from_args, tool.metadata, tool.collect_fn, a)
            return pd.read_csv(run_dir / "table0.tsv", sep="\t", index_col="name")

        printed = sapia_run("-i", "seq")
        check("Submitting 1 designs" in printed, "driver submits exactly ONE task")
        check(len(submits) == 1 and submits[0].returncode == 0,
              f"the one task succeeded ({[s.returncode for s in submits]}, {submits[-1].stderr[-300:]})")
        t = sapia_collect()
        check((t.loc[names, f"{LEAF}_status"] == "OK").all() and pd.isna(t.loc["nos", f"{LEAF}_status"]),
              "TSV: every sequenced row OK; the sequence-less row untouched")

        # The skill's filter expressions, against the TSV as prosapia reads it back.
        reps = t[(t.seqcluster_status == "OK") & (t.seqcluster_is_rep == True)]  # noqa: E712
        check(sorted(reps.index) == ["d3", "e2", "x"], f"'one design per cluster' filter ({list(reps.index)})")
        sizes = t[t.seqcluster_status == "OK"].groupby("seqcluster_cluster_id").size().to_dict()
        check(sorted(sizes.values()) == [1, 2, 3], f"groupby cluster_id diversity ({sizes})")
        big = t[(t.seqcluster_status == "OK") & (t.seqcluster_cluster_size >= 2)]
        check(sorted(big.index) == ["d1", "d2", "d3", "e1", "e2"], "cluster_size filter")
        # The skill's "best design per cluster" filter, verbatim, with a stand-in score.
        df_s = t.copy()
        df_s["boltz_iptm"] = pd.Series({"d1": 0.9, "d2": 0.5, "d3": 0.7, "e1": 0.2, "e2": 0.8, "x": 0.1})
        ok = df_s[df_s["seqcluster_status"] == "OK"]
        best = ok.sort_values("boltz_iptm", ascending=False).groupby("seqcluster_cluster_id").head(1)
        check(sorted(best.index) == ["d1", "e2", "x"], f"'best per cluster' filter ({list(best.index)})")
        # The skill's diversity reads, verbatim.
        ok = t[t["seqcluster_status"] == "OK"]
        desc = ok.groupby("seqcluster_cluster_id").size().describe()
        ratio = ok["seqcluster_n_clusters"].iloc[0] / ok["seqcluster_n_designs_clustered"].iloc[0]
        check(desc["count"] == 3 and abs(ratio - 3 / 6) < 1e-12, f"diversity reads ({ratio})")

        # Rerun WITHOUT --force: everything is re-clustered under a new token.
        tok1 = json.loads((run_dir / "table0" / LEAF / "seqcluster_task" / "run.json").read_text())["token"]
        printed = sapia_run("-i", "seq")
        tok2 = json.loads((run_dir / "table0" / LEAF / "seqcluster_task" / "run.json").read_text())["token"]
        check("Submitting 1 designs" in printed and "RE-CLUSTERED" in printed and tok1 != tok2,
              "rerun without --force re-submits the whole set")
        t = sapia_collect()
        check((t.loc[names, f"{LEAF}_status"] == "OK").all(), "rerun collects every row OK again")

        # A filtered run: 'x' leaves the set -> missing with its old labels wiped in the TSV.
        flt = tmp / "drop_x.py"
        flt.write_text("def apply_filter(df):\n    return df.drop(index=['x'])\n")
        sapia_run("-i", "seq", "-f", str(flt))
        t = sapia_collect()
        check(t.loc["x", f"{LEAF}_status"] == "missing"
              and all(pd.isna(t.loc["x", f"{LEAF}_{c}"]) for c in collect_mod.DATA_COLUMNS),
              f"filtered-out row: 'missing', every column NA in the TSV ({t.loc['x'].to_dict()})")
        check(int(t.loc["d1", f"{LEAF}_n_designs_clustered"]) == 5, "n_designs_clustered = 5 after the filter")

        # Changing -i under the same label is refused (the sequence-less row is fine).
        t2 = pd.read_csv(run_dir / "table0.tsv", sep="\t", index_col="name")
        t2["seq2"] = t2["seq"]
        t2.loc["d1", "seq2"] = None
        t2.to_csv(run_dir / "table0.tsv", sep="\t")
        n_before = len(submits)
        check(raises(lambda: sapia_run("-i", "seq2"), "d1", "-l/--dir-label") and len(submits) == n_before,
              "-i change stranding labelled rows -> refused before submit")
        printed = sapia_run("-i", "seq2", "-l", "v2")
        check("Submitting 1 designs" in printed, "...and accepted under a new dir label")
    finally:
        base_run.get_executor = real_get


def test_real_mmseqs(tmp: Path) -> None:
    binary = os.environ.get("SEQCLUSTER_TEST_MMSEQS")
    if not binary:
        SKIPPED.append("real-binary checks (set SEQCLUSTER_TEST_MMSEQS=/path/to/mmseqs)")
        return
    rng = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    names, seqs, fam_of = [], [], {}
    for f in range(4):
        base = [rng.choice(aa) for _ in range(80)]
        for c in range(3):
            s = list(base)
            for pos in rng.sample(range(80), 3):
                s[pos] = rng.choice([x for x in aa if x != s[pos]])
            name = f"fam{f} copy|{c}"
            names.append(name)
            seqs.append("".join(s))
            fam_of[name] = f
    rows, _, _ = build(table(dict(zip(names, seqs))), tmp / "real" / LEAF)
    rc = run_worker(rows[0], tmp / "real" / LEAF, binary)
    check(rc == 0, f"real mmseqs: worker exits 0 ({rc})")
    a = assignments(tmp / "real" / LEAF).set_index("name")
    check(set(a["status"]) == {"OK"} and set(a["n_clusters"]) == {"4"}, f"real mmseqs: 4 clusters ({set(a['n_clusters'])})")
    by_rep = a.groupby("cluster_rep").apply(lambda g: sorted({fam_of[n] for n in g.index}), include_groups=False)
    check(all(len(v) == 1 for v in by_rep) and len(by_rep) == 4,
          f"real mmseqs: every cluster is exactly one family ({by_rep.to_dict()})")
    check(set(a["cluster_size"]) == {"3"}, "real mmseqs: clusters of 3")
    up, _, _ = collect(table(dict(zip(names, seqs))), tmp / "real" / LEAF)
    check(all(v[f"{LEAF}_status"] == "OK" for v in up.values()), "real mmseqs: collector accepts the result")
    check(all(v[f"{LEAF}_mmseqs_version"] for v in up.values()), "real mmseqs: version recorded")


def main() -> None:
    tests = [test_layout_and_defaults]
    tmp_tests = [test_builder, test_worker, test_collector, test_task_script,
                 test_drivers_end_to_end, test_real_mmseqs]
    for t in tests:
        try:
            t()
        except Exception as e:  # noqa: BLE001
            FAILURES.append(f"{t.__name__} raised {type(e).__name__}: {e}")
    with tempfile.TemporaryDirectory() as d:
        for t in tmp_tests:
            sub = Path(d) / t.__name__
            sub.mkdir()
            try:
                t(sub)
            except Exception as e:  # noqa: BLE001
                import traceback
                traceback.print_exc()
                FAILURES.append(f"{t.__name__} raised {type(e).__name__}: {e}")
    print("=" * 70)
    for s in SKIPPED:
        print(f"SKIPPED: {s}")
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S) of {CHECKS} checks:")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print(f"all {CHECKS} seqcluster checks passed")


if __name__ == "__main__":
    main()
