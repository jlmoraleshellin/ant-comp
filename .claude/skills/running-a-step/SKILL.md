---
name: running-a-step
description: >-
  The executor-independent contract for running one prosapia step: the submit/wait/collect/verify loop, what counts as proof, the rerun-and-delete rule, what a worker may compute, and the mandatory report format with its Deviations block. Load in any *-worker agent before the first sapia command of a session.
---

# Running one step

This is what is true on **every** executor. Your agent definition holds the mechanics of yours — how a command reaches the compute, how you learn a task finished, how resources are requested. Where the two disagree, your definition wins and you should say so in your report.

You do not decide *what* to run. That comes from the caller. You decide *how*, run it, verify it, and report.

## The loop is four steps, not three

```
1 submit   →  capture the Output: line; the table is derived at submit time
2 wait     →  the executor's own completion signal, never the submit command
3 collect  →  -t <table from the Output: line>
4 verify   →  an invariant computed from raw files, NOT from the table
```

**Step 4 is not optional and it is the one that catches real corruption.** `Collected N row(s)` proves rows were written, not that the right work was done. Exit codes prove the tasks *ran*; some task scripts catch their own errors and still exit `0` (usalign and pyrosetta do). **`<leaf>_status == "OK"` is the only proof a design succeeded**, and even that is a claim by the collector.

So before reporting success, check something the tool did not compute:

| kind | example |
| --- | --- |
| counts | number of raw per-design result files == number of designs you meant to run |
| arithmetic | `n_res == target_len + binder_len`; `n_chains` after a merge; residue counts per entity |
| identity | a value matching an independent number — a length from the parent table, a residue name at a known position |

*Measured:* a template passed chain IDs, residue counts, byte-identical sequences, bit-identical coordinates and correct geometry while declaring 49 residues for its 41-residue chains — a phantom +8 shift that only an explicit count caught. Silent corruption has passed every status column and a full row count before now; only arithmetic caught it.

Say in your report **which invariant you checked and whether it held**. "I verified the output" is not a report; `48 result files for 48 designs; n_res 283 == 194 target + 89 binder on 3 spot-checked rows` is.

## `Submitting N designs` is not a design count

Several tools bin-pack before submitting, so N is the number of *manifest rows*:

| tool | what N counts |
| --- | --- |
| `boltz` | **shards** (`--shard-size`, default 10) — 2 designs print `Submitting 1 designs` |
| `proteinmpnn` | **parameter groups** — 8 backbones sharing flags print `Submitting 1 designs` |
| `cms` | **chunks** (`--designs-per-task`, default 100) |

N is therefore **not a usable guard against a filter failing open or shut**. Count the staged input files instead (`<out_dir>/boltz_inputs/*.yml`, `grp_*/inputs/*.pdb`, `cms_tasks/task_*.tsv`) and report that number. If a filter was used, report the count the filter itself printed as well.

**`No designs to submit.` exits 0.** If you do not see `Submitting N designs`, nothing was queued — usually a wrong `-i/--input-column`. Stop and report it; do not go on to collect.

**Staging dirs are not cleared between runs.** The manifest builders `mkdir(exist_ok=True)` and leave whatever was there, so a rerun in a run_dir that already held a pilot will re-process rows you meant to skip. Check the staged inputs are exactly the designs you intended before walking away.

## Rerun after a failure: delete the tool's output folder first

Delete `<run_dir>/<table>/<leaf>` — the tool's own output folder, nothing above it. **Never** delete the run_dir, `_registry.tsv`, or a table `.tsv`.

Why, each measured in this workspace:

- **Stale status artefacts read as success.** A rerun reuses the log dir, so a previous attempt's completion marker sits there looking like a finished task.
- **Different sharding conflicts.** Shard and task files are named by index (`shard_0.json`, `task_0.tsv`). Rerun with a different `--shard-size`, `--designs-per-task`, filter or design count and the new shards overwrite *some* of the old while orphans survive — a directory that is a silent mix of two runs.
- **Stale staged inputs get re-processed.**
- **Stale logs mislead.** Where log names carry a job id they sit *beside* the old ones rather than replacing them; always match on the id you just submitted.

**The exception — a partially successful run you are deliberately topping up.** Rerunning *without* `--force` resubmits only the non-`OK` rows, and that is the documented recovery. There, keeping the directory is the point.

> **failed attempt → delete. partial success you are topping up → keep and rerun without `--force`.**
> If you are unsure which you have, report the state and ask. Do not delete to find out.

## What you may compute, and what you may not

**If a measurement produces one value per design, it is a tool's job, not a script's.** When asked for one — "compute the interface contacts for each design", "score every backbone on X", "build a table of per-design distances" — **do not write it.** Reply with:

- that it should be a tool, because only a tool writes the number into the table where it can be filtered with `-f`, carries a `<leaf>_status`, and survives into child tables through lineage;
- the name of any existing tool that already produces it — read the **collector's column list**, do not guess (`cms` writes per-residue interface contributions; `pyrosetta` writes interface energetics; `usalign` writes TM/RMSD);
- then **stop**. You have no `Agent` tool, so you cannot commission one. The caller can.

**What you may still do:**

- **Read-only inspection** — row counts, file counts, reading a log, checking an invariant, printing a few columns, confirming a residue identity at a known position. *A fact about the run is yours; a number about a design is a column.*
- **Filter modules.** A filter selects rows, it does not measure them, and it must select on columns that already exist. Keep them in a subfolder of the run_dir (`<run_dir>/filters/`) so they don't clutter it, and report the **md5 and the row count the filter printed** — the filter's name alone does not say which version ran.
- **Batch-level statistics are not yours either**, but for the opposite reason: they belong to `table-analyst`, which the caller invokes. If asked for a correlation or a distribution, say so and hand the column names back.

## The report

Report every time, in this shape:

```
report file  <scratchpad>/report_boltz_table1.md   (first line of the hand-back)
run_dir      outputs/20260930_083738_binder_A
table        table1_biasC          (from the Output: line, never guessed)
command      sapia run boltz … (as sent, including every flag)
submitted    55 designs  (staged input files; `Submitting N` said 6 shards)
collected    55 rows      (expected 55)
status       boltz_cofold_status: OK 54, error 1, missing 1079
failures     1 task non-zero — <tail of its .err>
invariant    55 result dirs for 55 staged inputs; n_res 283 == 194+89 on 3 rows — HELD
odd          <anything that looked wrong even though it succeeded>
Deviations:  none
```

### Write the report to a file before you hand it back

**`Write` the full report to the session scratchpad first, then hand back — and make the file's absolute path the first line of the hand-back.** Name it `report_<tool>_<table>.md`, or `report_<tool>_<table>_2.md` if that name is taken, and never overwrite an earlier report.

This is not bookkeeping. **A hand-back carries exactly one report per worker run**: if you were resumed for a second mission, your second hand-back is refused as already delivered, and you have no `SendMessage` to fall back on. The caller then sees a bare "completed" with no numbers, and the only recovery is digging the text out of a huge JSONL transcript (measured 2026-10-12, on a worker resumed ~5 times). A report on disk survives that, and survives the caller's context being compacted.

So: **write the file even when the hand-back looks like it will work**, and if you are resumed for a second mission, say in your first sentence that the earlier report is at its own path — the caller may never have received it.

### The Deviations block is mandatory

End every report with `Deviations:` — either the literal word `none`, or one line per place you ran **less than, or other than, what the caller asked**. Subsampled, sharded differently, a filter you added, a flag you changed, a tool you substituted, rows you skipped, a step you could not reach.

This exists because the caller will turn your numbers into a scientific claim, and the difference between "no design passed" and "no design *of the 55 we ran* passed" is the difference between a finding and an error. A deviation you do not declare becomes a conclusion nobody can audit.

**State unknowns as unknowns.** "I could not determine the partition load" is useful; a confident guess is worse than nothing. Never describe a run as successful when tasks failed, the row count is short, or `<leaf>_status` is not `OK`.

### When something fails in a way you do not understand

**Stop and report the raw evidence** rather than patching a tool or retrying blind. Root-cause it if you can — re-run the task script by hand, read the source — but let the caller decide the fix. Then stop: **do not chain into the next tool unless you were asked to.** The caller decides what comes next.

## Two failure signatures worth recognising anywhere

**A task that fails with `.out` AND `.err` both empty died before its first statement produced output.** Look at the prelude and the variable assignments, not the worker. The classic cause is a local named after a bash special variable — `GROUPS=$(…)` fails with rc=1, is silently discarded, and `set -euo pipefail` then kills the task with no diagnostic anywhere. Same hazard: `UID`, `EUID`, `PPID`, `PIPESTATUS`, `SECONDS`, `RANDOM`, `LINENO`, `IFS`, `PATH`. Reproduce by re-running the task script by hand under `bash -x` with the real task environment. A missing activation script presents the same way.

**`collect` takes no filter**, so it walks the whole table and stamps `missing` on rows that were never submitted. Expected — but it means `<leaf>_status` is non-empty for rows you deliberately skipped, and a `missing` count is not a failure count.

## Skills

Load **`prosapia`** before your first `sapia` command (your definition may preload it). It is the workbench contract: tables and lineage, `create` vs `update` and how the output table is derived, the base run/collect flags, labels, the ready set, and the traps that make a run silently submit nothing.

Then load the **tool's own skill** before composing its flags. The tool skills were written against Modal; read them with that filter if you are on another executor — tool flags, input columns, collected columns, bin-packing and every scientific trap are still true; anything about Volumes, images, `modal-shell`, `--gpu-type`, build times or `.exit` files is not.

`sapia run <tool> --help`, run on your executor, is authoritative for flags. **Don't guess flag names.**
