---
name: running-a-step
description: >-
  The executor-independent contract for running one prosapia step: the submit/wait/collect/verify loop, what counts as proof, labels, filters and merging, the rerun-and-delete rule, what a worker may compute, and the mandatory report format with its Deviations block. Load in any *-worker agent before the first sapia command of a session.
---

# Running one step

True on every executor. Your agent definition holds the mechanics of yours; where the two disagree, it wins and you say so in your report. You do not decide what to run; you decide how, run it, verify it, and report.

## The loop

```
1 submit   →  capture the Output: line; the table is derived at submit time
2 wait     →  the executor's own completion signal, never the submit command
3 collect  →  -t <table from the Output: line>
4 verify   →  an invariant computed from raw files, NOT from the table
```

Step 4 is not optional. `Collected N row(s)` proves rows were written, not that the right work was done; exit codes prove tasks ran, and some task scripts catch their own errors and exit `0`; `<leaf>_status == "OK"` is the collector's claim. Silent corruption passes every status column and a full row count, and only arithmetic catches it. Check something the tool did not compute:

| kind | example |
| --- | --- |
| counts | raw per-design result files == designs you meant to run |
| arithmetic | `n_res == target_len + binder_len`; `n_chains` after a merge |
| identity | a length from the parent table, a residue name at a known position |

Say which invariant you checked and whether it held: `48 result files for 48 designs; n_res 283 == 194 + 89 on 3 spot-checked rows`, never "I verified the output".

## Submit

**`Submitting N designs` is not a design count**: several tools bin-pack, so N counts manifest rows (`boltz` and `rfdiffusion3`: shards of `--shard-size`, default 10; `proteinmpnn`: parameter groups; `cms`: chunks of `--designs-per-task`, default 100). N cannot guard against a filter failing open or shut. Count the staged inputs (`<out_dir>/boltz_inputs/*.yml`, `rfd3_inputs/*.json` keys, `grp_*/inputs/*.pdb`, `cms_tasks/task_*.tsv`) and report that number and the row count the filter printed.

- **`No designs to submit.` exits 0.** No `Submitting N designs` line means nothing was queued (usually a wrong `-i`): stop and report, do not collect.
- **Staging dirs are not cleared between runs**, so a rerun in a run_dir that held a pilot re-processes rows you meant to skip. Check the staged inputs are exactly the intended designs.
- **A chunked submit can die partway.** Submits split into arrays of 1000 tasks; if `sbatch` fails at chunk k, chunks 1 to k-1 stay queued. Do not cancel or resubmit blindly: let them drain, top up with a `-f` filter keeping only rows with no per-design output yet, collect once after both rounds, report both submits.
- **A long manifest build runs inside the submit call** (hundreds of CIF-to-PDB conversions take many minutes before `sbatch`). Run it detached with a walltime that covers the build, log to a file, find the job id from the scheduler (`squeue -u $USER`), then apply the array workarounds in `docs/vib.md`.
- **A failed or refused submit leaves debris.** The registry row (`create` tools) and the leaf dirs are written before the manifest is built, so an error leaves a registry row with no `<table>.tsv` (a phantom, not an empty table) and an empty leaf dir; an `update` run that queues nothing still rewrites its leaf's `.meta.json`. After any submit without a `Submitting N designs` line, compare `_registry.tsv`, the `*.tsv` files and the leaf dirs, report mismatches, and remove a phantom only with the caller's agreement.

## Collect

- **A collect rewrites the whole table**: it reads the full `.tsv` and writes it back. Copy the table aside first, keep one writer per table, and parse numbers rather than comparing strings (a column that gained blanks is re-serialised: `44` becomes `44.0`).
- **`collect` takes no filter**, so it stamps `missing` on every row with an input that was never submitted. `missing` is not a failure count; a row with a blank input cell is not touched.

## Labels, filters and row selection

- **One label per parameter set.** A rerun of the same recipe keeps its label (with `--force` to redo `OK` rows); a different input column or flag set gets a new label; a verification run gets its own. Each label adds a leaf dir and a full column set, and a reused label silently overwrites that leaf's `.meta.json`, so list the labels under the table before any run. An `update` tool on a filtered table needs no label unless the same rows are rerun with different parameters. `--table-label` is the same rule for forking a `create` lineage.
- **`-f` or `-w`.** `-f <module.py>` runs `apply_filter(df) -> df`; `-w column=value` (repeatable, after `-f`) keeps rows whose column equals the value as a string, so a re-serialised number stops matching. A filter should assert the expected row count and print an md5 of the selected names; report both. A missing or empty filter file errors at submit after the table is registered (`FileNotFoundError`, or `does not define an 'apply_filter' function`): check it is non-empty right after upload and read the submit output.
- **Pooled rows from several sources.** Before splitting an `update` tool into runs per source, check how it treats blank input cells: the ready set skips them, and a tool that resolves its input up the lineage may fill them from an ancestor. Count rows per source in the staged inputs.

## Merging

Check `sapia merge -h` and the installed source, not this text. The installed prosapia merges siblings only: sources must share a parent (otherwise it refuses), row names get the source's label as suffix (an unlabelled source keeps bare names; a collision refuses), columns are the union, and rows are copied, so annotations on the merged copy never reach the sources. Each row carries `merge_source`. `sapia merge` takes whole tables and has no row filter, so merging after selecting is not possible: select at run time with `-w merge_source=<table>` or `-f`. Names are unique only inside a run_dir; a set pooled across runs is keyed on (run, name).

Pattern: fork with `--table-label`, merge the siblings right after the fork, run each `update` tool once on the merged table, use `-w merge_source=...` for per-fork variants. A merge is a ledger step; read the merged header and row count before building on it, and re-check any name-suffix test afterwards.

## Rerun after a failure

Delete `<run_dir>/<table>/<leaf>`, the tool's own folder and nothing above it; never the run_dir, `_registry.tsv` or a table `.tsv`. A rerun reuses the log dir, so an old completion marker reads as success; shard and task files are named by index, so a different `--shard-size`, filter or design count leaves orphans mixed with the new run; stale staged inputs get re-processed; and logs carrying a job id sit beside the old ones, so match on the id you just submitted.

> failed attempt → delete. partial success you are topping up → keep and rerun without `--force`. If unsure which, report the state and ask; do not delete to find out.

- **The top-up is for `update` tools only.** A `create` tool never skips finished rows: a rerun under the same `--table-label` resubmits everything and reuses `grp_N` / `task_N` names, overwriting finished outputs. Clean top-up: a new `--table-label` with a `-f` filter selecting only the failed rows (verified equal to the staged inputs of the failed groups), then `sapia merge` of the two sibling tables.
- **A first slice of a `create` run may not be completable under its label.** Drivers that cut groups from 0 on every submit and do not prune put the slice's stems in two groups on a full run (double collect, or `error:` rows); some refuse such a leaf instead, and the tool's skill says which. Run the whole set at once when cheap, or have the caller delete the slice leaf and its child `.tsv`.
- **A `create` tool on a table that already carries its own leaf's `<leaf>_status` = `OK` submits nothing** (the ready set drops `OK` rows unless `--force`): it prints `No designs to submit.`, exits 0, and still registers the child table. Pass `--force` with the `-f` filter and keep the label, since the phantom row is already that label's table.
- A create leaf's `.meta.json` records input column, label, filter and `--where`, not the tool parameters; those live in the staged `<leaf>_tasks/task_N.tsv` and the ledger command, so write the command into the ledger first.

## What you may compute

**A measurement with one value per design is a tool's job, not a script's.** Reply that it should be a tool (only a tool writes the number into the table, filterable with `-f`, with a `<leaf>_status`, surviving into child tables), name any existing tool that already produces it (read the collector's column list, do not guess), and stop; the caller commissions tools. You may do read-only inspection (counts, a log, an invariant, a few columns); a fact about the run is yours, a number about a design is a column. Filter modules select on existing columns and never measure: keep them in `<run_dir>/filters/` and report md5 and printed row count. Batch statistics belong to `table-analyst`: hand the column names back.

## The report

```
report file  <scratchpad>/report_<tool>_<table>.md   (first line of the hand-back)
run_dir      <as printed by new_run>
table        <table>          (from the Output: line, never guessed)
command      sapia run <tool> … (as sent, every flag)
submitted    55 designs  (staged input files; `Submitting N` said 6 shards)
collected    55 rows      (expected 55)
status       <leaf>_status: OK 54, error 1, missing 1079
failures     1 task non-zero: <tail of its .err>
invariant    55 result dirs for 55 staged inputs; n_res 283 == 194+89 on 3 rows: HELD
odd          <anything that looked wrong even though it succeeded>
Deviations:  none
```

- **Write the full report to the session scratchpad first** (`report_<tool>_<table>.md`, `_2` if taken, never overwrite) and make its absolute path the first line of the hand-back: a hand-back carries one report per worker run, so a second mission's is refused as delivered, and the file survives that and context compaction. If resumed, say in your first sentence that the earlier report is at its path. If the harness refuses the `Write`, return the full report as text, keep command files in the scratchpad, and say in the first line that no report file exists.
- **Every worker also writes the ledger-brief file its agent definition names** (`.claude/agents/*-worker.md` give its path and content); hand back with both.
- **`Deviations:` is mandatory**: `none`, or one line per place you ran less than, or other than, asked (subsampled, sharded differently, a filter or flag you added, a tool substituted, rows skipped, a step unreachable). The caller turns your numbers into a claim, and "no design passed" differs from "no design of the 55 we ran passed". State unknowns as unknowns; never call a run successful when tasks failed, rows are short, or `<leaf>_status` is not `OK`.
- **When a failure is not understood**, stop and report the raw evidence rather than patching a tool or retrying blind; root-cause it if you can (re-run the task script by hand, read the source), let the caller decide the fix, and do not chain into the next tool unless asked.

## One failure signature

A task with `.out` and `.err` both empty died before its first statement: look at the prelude, variable assignments and activation script, not the worker. The classic cause is a local named after a bash special variable (`GROUPS=$(…)` fails silently under `set -euo pipefail`; the list is in `authoring-a-tool`). Re-run the task script by hand under `bash -x` with the real task environment.

## Skills

Load `prosapia` before your first `sapia` command, then the tool's own skill before composing its flags. Tool skills were written against Modal: flags, columns, bin-packing and scientific traps hold everywhere; Volumes, images, `modal-shell`, `--gpu-type`, build times and `.exit` files do not. `sapia run <tool> --help` on your executor is authoritative for flags; do not guess flag names.
