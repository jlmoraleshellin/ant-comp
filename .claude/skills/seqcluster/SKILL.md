---
name: seqcluster
description: >-
  Which designs in a table are sequence-redundant? Clusters ONE single-chain sequence column with
  MMseqs2 easy-cluster in a SINGLE task over the whole table, and writes per design its cluster_id,
  representative, cluster size and is_rep. Load before composing a seqcluster run, before keeping one
  design per cluster or measuring sequence diversity, and before reading any seqcluster_ column. The
  traps: every label is a property of the SET (cluster_id is local to one run, never comparable
  across runs or dir labels); every run re-clusters the whole table, so --force is always on and a
  -f filter WIPES labels on the rows it excludes; is_rep is MMseqs2's pick, not the best design.
  Covers the required -i, the MMseqs2 defaults it passes through (min-seq-id 0.0, coverage 0.8),
  why -c is NOT coverage, the multi-chain refusal, and the run token that turns a stale result into
  an error.
---

# seqcluster

Sequence-identity clustering of a table's designs with MMseqs2 `easy-cluster`, one label set
per design. **`action: update`**: it annotates the table it runs on and mints no child. `-t`
and `-i` are required. CPU-only: the builder forces `gpus_per_task = 0`, so no `-g 0` is
needed.

> **Status: reviewed, not yet verified on real data.** Built 2026-10-10. The test suite passes
> (109 checks with a stub MMseqs2, 115 with the real binary on synthetic sequences).
> `tool-reviewer` (2026-10-10) verdict: **safe to proceed to real-data verification for a
> single run under a FRESH `-l` label with one `-i`** — exactly the shape of a first campaign
> use. It found two real risks that fire only on a RERUN under an EXISTING label (see "Known
> risks, not yet fixed" below) — read that before ever reusing a dir-label on this tool. Still
> needs a worker run under a `verification` label before a campaign uses its columns. See
> *What the tests do not establish* below.

## Known risks, not yet fixed (found by `tool-reviewer`, 2026-10-10)

Both require a **rerun under an already-used `-l` label** to fire — a first, single run under a
fresh label (the shape every verification and production use so far has taken) is not exposed to
either.

1. **A rerun with `-i` switched to a DIFFERENT column of the same row-count and same per-row
   lengths can silently cluster the OLD sequences under the NEW column's name.** Only the
   manifest (`run.json`'s token/params/input-column/count) is guaranteed fresh across a Modal
   task boundary — prosapia publishes just the manifest via `batch_upload`, and `members.tsv`
   is read by the task through its own Volume mount, which can be outrun by `spawn_map` and see
   a stale file. The collector currently checks names, count and `seq_len` against the
   workstation's copy, but **never the sequence content itself** — so a same-shape swap passes
   every check and reports `OK` with the new column's name while having clustered the old
   sequences. **Fix not yet applied**: hash `members.tsv` (sha256) into both the manifest and
   `run.json`, have the worker refuse on a mismatch and echo the hash it actually read.
2. **A row excluded by `-f` that lacks the new `-i` column only gets a printed WARNING, not a
   reset to `missing`**, when its labels are from an earlier run under the same `-l`. Two
   different clusterings can then sit under one set of columns: `cluster_id 3` from run 1 and
   `cluster_id 3` from run 2 both read `status == OK` and would merge in any `groupby` on
   `seqcluster_cluster_id`. **Fix not yet applied**: refuse outright, or explicitly reset those
   rows to `missing` in the collector factory, rather than warning and continuing.

**Until these are fixed: never rerun this tool under a dir-label it has already written to.**
Always use a fresh `-l` for every distinct clustering (which costs nothing — MMseqs2 reclusters
the whole set every time regardless).

## Reading the output: min-seq-id 0.0 answers a homology question, not a near-duplicate one

**Flagged by `tool-reviewer`, not a code defect — a reading caveat for whoever reports
`n_clusters`.** At the MMseqs2 defaults (`--min-seq-id 0.0`, 80% two-sided coverage, E-value
`1e-3` not exposed), `n_clusters` counts **families of detectable homology**, not near-duplicate
sequences. ProteinMPNN sequences redesigned on related backbones against one shared target can
share detectable homology at roughly 30–50% identity and collapse into a handful of clusters —
reading as "highly redundant" when the designs are not actually near-copies of each other.
**A diversity number from this tool is only meaningful alongside the `--min-seq-id` it was run
at.** Run at least one additional, higher threshold (e.g. 0.9, and optionally 0.5) under a
separate `-l` before reporting redundancy, and always quote the threshold next to the number.

## How it differs from every other tool here

**One task for the whole set.** A design's cluster depends on every other sequence it was
clustered with, so the builder hands MMseqs2 all designs at once. It returns exactly one
manifest row, and the driver prints **`Submitting 1 designs`** however many designs are
clustered. That line is correct. It is not a filter problem.

**Every run re-clusters the whole table, so `--force` is always on.** Every other tool
skips rows that are already `OK`. This one cannot: re-clustering only the new rows would
cluster a subset on its own and give them cluster IDs that collide with the earlier run's.
So:

- the **run** clusters every row of the (`-f`-filtered) table that has a sequence, `OK` or
  not, and prints `rows already OK from an earlier run will be RE-CLUSTERED` when there were
  some;
- the **collect** re-annotates every row with a sequence, whether `--force` is passed or not;
- to **keep** an earlier clustering, run the new one under a different `-l/--dir-label`.
  Nothing else preserves it.

## Premise, and when its numbers are meaningless

It assumes **one amino-acid sequence per row, one chain, one-letter code**, in a column of
**this** table. The answer describes **that set of sequences, and nothing else**.

| Situation | What it does |
| --- | --- |
| `-i` not given, or not a column of this table | **Refused at submit**, listing the table's sequence-like columns. There is **no lineage walk**: a sequence inherited from a parent would be shared by every child, so it would cluster them all together |
| A row holds a **multi-chain** string (`/` or `:`) | **Refused at submit**, naming the designs. ProteinMPNN joins the chains of a multi-chain design with `/`, so a complex column carries the target too, and a target shared by every design makes them all cluster together. Pass the binder-only column |
| A character outside `ACDEFGHIKLMNPQRSTVWY` + `XBZJUO` (a gap `-`, a stop `*`) | **Refused at submit**, naming the designs and characters |
| Lower-case letters | Upper-cased silently, then clustered |
| A row with no sequence (NaN or blank) | Not clustered. Collects as `missing` if the collector visits it (see the traps) |
| **The column is the wrong chain** (the target instead of the binder, in a table where each row has one chain) | **`OK`, with tidy clusters of the wrong protein.** The tool cannot detect this. Read `seq_len` |
| **A different set** (another `-f`, more designs added later) | **`OK`, with different labels.** Representatives, IDs and sizes can all change when one design is added or removed |
| A single design | `OK`: one cluster of size 1, `is_rep` True |

It always returns labels for a valid input. **There is no input for which it reports that
the clustering is meaningless for your question.**

## Trust columns: read these before the science

| Column | Expected value | What a surprise means |
| --- | --- | --- |
| `seqcluster_n_designs_clustered` | **The number of designs you meant to cluster**, the same on every member row | Smaller than you expected: a `-f` filter, or rows with a blank sequence. Record the expected count in the campaign log |
| `seqcluster_seq_len` | **The length of the chain you meant to cluster** (the binder length the generator was asked for) | A target-sized or complex-sized length means the wrong column. Kept on `error:` rows |
| `seqcluster_input_column` | The `-i` you passed | Provenance: the column the sequences came from. Kept on `error:` rows |
| `seqcluster_params` | `min_seq_id=0.0;c=0.8;cov_mode=0;cluster_mode=0` at the defaults | The exact MMseqs2 parameters. Two tables clustered with different values here are not comparable. Kept on `error:` rows |
| `seqcluster_mmseqs_version` | `8cc5ce367b5638c4306c2d7cfc652dd099a4643f` (release 18-8cc5c, pinned in the Modal image) | A different hash means a different binary |

## The question it answers, and the threshold

**How much sequence redundancy is in this set of designs, and which designs are
near-copies of each other?** Typical uses:

- keeping one design per sequence family before an expensive step (Boltz, AF3, wet lab);
- measuring diversity: `n_clusters` against `n_designs_clustered`;
- checking whether a sampling temperature or a bias actually spread the sequences out.

**This skill gives no threshold, deliberately.** The parameters default to **MMseqs2's own
`easy-cluster` defaults**, read from `mmseqs easy-cluster -h` of the pinned binary. They are
not a recommendation:

| Parameter | Default | What the default means |
| --- | --- | --- |
| `--min-seq-id` | **`0.0`** | **No identity floor.** Two sequences cluster whenever an alignment passes the coverage test and the E-value cutoff (`1e-3`, not exposed). Distantly related designs can share a cluster. **For a "near-copies" question, set this yourself (for example 0.9), and record the value and the reason in the campaign log** |
| `--coverage` | `0.8` | At least 80% of residues aligned, as counted by `--cov-mode` |
| `--cov-mode` | `0` | Coverage of both query and target |
| `--cluster-mode` | `0` | Set-Cover (greedy). `1` is connected components (transitive, so chains of neighbours merge into large clusters). `2`/`3` are CD-HIT-like |

## Invocation and flags

```bash
sapia run seqcluster <run_dir> -t table1 -i proteinmpnn_redesign_full_sequence --min-seq-id 0.9
# wait for the ONE task, then:
sapia collect seqcluster <run_dir> -t table1
```

Not yet run on Modal or on vib: this invocation is the intended one, **not a verified one**.

| Flag | Default | Meaning |
| --- | --- | --- |
| `-i/--input-column` | **none, required** | A single-chain amino-acid sequence column of **this** table |
| `--min-seq-id` | `0.0` | MMseqs2 `--min-seq-id`, 0.0–1.0 |
| `--coverage` | `0.8` | MMseqs2 `-c`, 0.0–1.0. **There is no `-c` short form here**: see the traps |
| `--cov-mode` | `0` | MMseqs2 `--cov-mode`, 0–5 |
| `--cluster-mode` | `0` | MMseqs2 `--cluster-mode`, 0–3 |
| `-l/--dir-label` | empty | Output dir / column suffix. **The only way to keep two clusterings side by side** |
| `-f/--filter` | none | Restricts the set that is clustered. Rows it excludes are **reset to `missing`** at collect |
| `--force` | — | **Has no effect.** Every run and every collect already covers the whole set |

Collect takes only `-t`, `-l` and `--force` (the last ignored).

Resources: 4 CPU, 8 GB, 1 h, set both in `modal_image.py` `RESOURCES` (Modal) and in the
`.sh` `#SBATCH` lines (SLURM). The worker uses `$SEQCLUSTER_THREADS` if set, else
`$SLURM_CPUS_PER_TASK`, else the CPUs the process can see. Thread count changes speed, not
the result: deterministic across `--threads 1` and `--threads 4` on a synthetic set.

**(Modal only)** The image is `debian_slim` plus the static MMseqs2 `18-8cc5c` **sse4.1**
build, pinned by URL and sha256 (`0ece6a6a…7b08`). There is no pip layer. sse4.1 rather
than avx2 so that an older host CPU cannot kill the task with an illegal instruction. The
build runs `mmseqs version` and a two-sequence `easy-cluster` as a smoke check. The first
submit builds the image inside the submit call.

**(vib)** **Not runnable.** There is no `SAPIA_ACTIVATE_SEQCLUSTER` and MMseqs2 is not
installed on the cluster. The tool registers in `sapia run --help` there, and every task
would then die at activation.

## Traps, each with its signature

- **`Submitting 1 designs` for a table of thousands.** This is correct: one task holds the
  whole set. The true count is in the run's own line `seqcluster: clustering N designs in
  ONE task` and, after collect, in `seqcluster_n_designs_clustered`.
- **`-c 0.9` sets CPUs, not coverage.** `-c` is prosapia's `--cpus-per-task` on every tool.
  Coverage is `--coverage` only. The signature is `seqcluster_params` still reading `c=0.8`
  after you meant to change it.
- **A `-f` filter wipes earlier labels on the rows it excludes.** Collect resets every row
  outside this run's set to `missing` with every column NA, deliberately, so that one
  clustering's labels never sit beside another's. The signature is rows that were `OK` now
  reading `missing`. To cluster a subset **and** keep the full-table labels, use a new `-l`.
- **Changing `-i` under the same dir label is refused at submit** when some row carries
  earlier labels but has no value in the new column: `N row(s) carry results in
  'seqcluster_status' from an earlier seqcluster run but have no '<col>'`. Collect only ever
  visits rows with a present input column, so those labels could not be reset. Use a new
  `-l`. (If the row is filtered out by `-f` the builder cannot see it. Collect then prints
  `WARNING: seqcluster: N row(s) keep labels from an EARLIER run` instead.)
- **Collecting before the task finishes** gives every member the **same** `error: no
  seqcluster_result/assignments.tsv -- the clustering task has not finished...`. It is not
  `missing`. Re-collect when the task is done: collect is always a full re-collect, so the
  error rows are replaced.
- **A result from an earlier submit is never collected as this one's.** The builder writes
  a fresh run token into `seqcluster_task/run.json`. The worker echoes it into every
  `assignments.tsv` row, and the collector compares the two. A mismatch gives every member
  `error: stale result -- assignments.tsv token ... is not this run's ...`.
- **The task is all-or-nothing.** An MMseqs2 failure, or a cluster TSV that is not a clean
  partition of the submitted set (a design missing, a design in two clusters, an unknown
  ID), gives **every** member the same `error:` and the task exits 1. The raw MMseqs2 log is
  `seqcluster_result/mmseqs.log`.
- **An unreadable member list makes collect raise** (`cannot collect -- this run's member
  list is unreadable`) rather than mark rows `missing`, because it no longer knows who was
  submitted. Re-run under the same `-l`.
- **Design names never reach MMseqs2.** FASTA headers are `>s<i> <name>`. MMseqs2 keys only
  on the first word, `s1, s2, …`, and the worker maps those back to table names from its own
  member list. Names with spaces or `|` are safe. `seqcluster_result/clusters_cluster.tsv`
  therefore holds `s<i>` IDs, not names; `input.fasta` beside it holds the mapping.

## What it collects

Bare names below are prefixed `seqcluster[_<label>]_`. `seqcluster_status == "OK"` is the
only proof a row was clustered.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `status` | `OK`, `error: ...` or `missing` | `OK`: a member of the validated clustering. `error:` is the **same on every member** of the run, with science NA and the trust columns kept. `missing`: **not in this run's set** (filtered out, a blank sequence, or no run under this label), with **every** column NA |
| `path` | `seqcluster_result/assignments.tsv` | One file for the whole run, the same on every member row. Empty on `missing` |
| `n_designs_clustered` | **Trust.** Designs in the one task | See the trust table |
| `seq_len` | **Trust.** This design's clustered length | See the trust table |
| `input_column` | **Trust.** Source column | See the trust table |
| `params` | **Trust.** MMseqs2 parameters | See the trust table |
| `mmseqs_version` | **Trust.** Binary build hash | See the trust table |
| `cluster_rep` | `name` of the cluster's representative | Use it to see who a design was grouped with. It is MMseqs2's set-cover centre, **not the best design** |
| `cluster_id` | 1-based integer per cluster | **Local to this one run**, numbered by first appearance in MMseqs2's output (**not** by size). Never compare across runs, tables or dir labels. Group on it within one run only |
| `cluster_size` | Members of this cluster within this run's set | `1` is a singleton: no other design within the thresholds |
| `is_rep` | True on the representative's own row | Exactly one per cluster. Reads back from the TSV as booleans, so `== True` works |
| `n_clusters` | Distinct clusters in this run's set | The same on every row. `n_clusters / n_designs_clustered` is the non-redundant fraction |

On disk, under `<run_dir>/<table>/seqcluster[_<label>]/`: `seqcluster_task/members.tsv`
(what was submitted), `seqcluster_task/run.json` (token, params, input column, count), and
`seqcluster_result/` (`assignments.tsv`, `input.fasta`, `clusters_cluster.tsv`,
`clusters_rep_seq.fasta`, `mmseqs.log`).

## Filter expressions

Keep MMseqs2's representative of every cluster:

```python
def apply_filter(df):
    return df[(df["seqcluster_status"] == "OK") & (df["seqcluster_is_rep"] == True)]  # noqa: E712
```

Keep the **best** design per cluster instead, ranked by a science column you trust. This is
usually what is wanted, because the representative is not chosen for quality:

```python
def apply_filter(df):
    ok = df[df["seqcluster_status"] == "OK"]
    return ok.sort_values("boltz_iptm", ascending=False).groupby("seqcluster_cluster_id").head(1)
```

Measure diversity, read in a notebook or by `table-analyst`, not as a filter:

```python
ok = df[df["seqcluster_status"] == "OK"]
ok.groupby("seqcluster_cluster_id").size().describe()
ok["seqcluster_n_clusters"].iloc[0] / ok["seqcluster_n_designs_clustered"].iloc[0]
```

All three were run against a table written by the real prosapia collect driver in the test
suite (with a stub MMseqs2).

## What it is blind to

- **Structure and fold.** Two designs in different clusters can share a fold, and two in
  one cluster can differ at every interface position that matters. For backbone
  similarity, use `usalign`.
- **Where on the target a design binds.** For the epitope, use `epitope`.
- **Quality.** `is_rep` and `cluster_rep` are MMseqs2's choice. Rank with a science column
  (Boltz confidence, `pyrosetta` `if_dG`, `cms`).
- **Anything outside the set.** It never compares against earlier rounds, other tables or
  natural sequences. To cluster across tables, build one table that holds them all.
- **Clustering errors that MMseqs2 itself permits.** With `--cluster-reassign` off (the
  MMseqs2 default, not exposed here), MMseqs2's own help says cascaded clustering *"can
  cluster sequence that do not fulfill the clustering criteria"*. A member may therefore sit
  slightly outside `--min-seq-id` of its representative. The tool checks that the result is
  a clean partition, not that every pair meets the thresholds.

## What the tests do not establish

`tools/seqcluster/test_seqcluster.py` exercises the builder, worker, collector and `.sh`,
and runs the real prosapia run/collect drivers with only the executor faked. With
`SEQCLUSTER_TEST_MMSEQS=/path/to/mmseqs` it also runs the real binary on 4 synthetic
families × 3 near-copies and gets 4 clusters of 3. It cannot establish:

- whether `proteinmpnn_redesign_full_sequence` (or any real input column) holds **one chain
  per row** in a real table. It is expected to be a single redesigned chain in the
  current campaign, but **no one has checked real data**. A `/` would be refused at submit;
- whether the pinned sse4.1 binary runs on the Modal host the task lands on (the build only
  smoke-tests it on the build host);
- whether the defaults group real ProteinMPNN redesigns the way the caller expects;
- whether `os.sched_getaffinity` in a Modal container reports the 4 allotted CPUs or the
  whole host (this affects speed only).
