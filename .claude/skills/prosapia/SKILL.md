---
name: prosapia
description: >-
  The prosapia workbench model and the `sapia` CLI contract: run_dirs, tables and lineage, create vs update, the base run/collect flags. Load before composing any sapia command, on any executor. Trap: several mistakes make a run exit 0 having submitted nothing. Executor mechanics live in the worker agents, per-tool flags in the tool skills.
---

# prosapia

`prosapia` (CLI `sapia`) is a **workbench, not a pipeline**. There is no fixed order of steps: there
is one tabular data format and tools that consume and produce it. Your job is to compose one step
correctly, run it, and report what landed in the table.

Two principles define how tools meet the data:

1. **A design is a row; a generation of designs is a table.**
2. **When a protein diverges in sequence or structure it is no longer the same protein** — it is a
   child of a parent, and therefore needs a new table.

**This page is executor-neutral.** How a command reaches the compute, how you learn a task finished,
and how resources are requested belong to your worker agent and to `running-a-step`. Nothing here
assumes Modal or SLURM.

Authoritative references: **the installed package does not ship `docs/`** — only `src/`. Query
<https://github.com/jlmoraleshellin/prosapia/tree/dev> (the `dev` branch specifically), plus
`sapia run <tool> --help` **on your executor**, which is always authoritative for a tool's own flags.
**Don't guess flag names.**

| Question | Page |
| --- | --- |
| Base `sapia run` flags, executors, the task contract | `docs/running-a-tool.md` |
| `sapia collect` | `docs/collecting-a-tool.md` |
| `create` vs `update`, roots, lineage, `lookup` | `docs/lineage-and-tables.md` |
| `-l/--dir-label` vs `--table-label` | `docs/using-labels.md` |
| The workstation, tool images, `.exit` files | `docs/running-on-modal.md` |
| `-f/--filter` | `docs/writing-a-filter-function.md` |

## The data model

```
run_dir/
├── _registry.tsv      catalog of every table: gen, table_label, parent_table, tool, created_at
├── table0.tsv         one row per design, keyed by `name`
├── table1.tsv         a child table (the next generation)
├── .manifests/        transient per-run manifests
└── table0/            per-table, per-tool outputs on disk
    └── rfdiffusion3/  = run_dir/<table>/<leaf>/
```

- **Each tool contributes columns to a row**, never files to another tool.
- **Columns are leaf-prefixed**: the leaf is `<tool>`, or `<tool>_<dir_label>` if you passed `-l`.
  Every tool writes `<leaf>_status` and usually `<leaf>_path`.
  **`<leaf>_status == "OK"` is the only proof a design succeeded.**
- **Row lineage** lives in `parent_name`, `parent_table`, `gen`. A root table's `parent_table` is the
  sentinel `root`.
- **Only `sapia new_run` mints a run_dir.** No tool ever does.
- **`_registry.tsv` is the authority on which tables exist** and how they relate. A remembered
  lineage tree drifts; read the registry.

## The verbs

`sapia` has exactly these: `new_run`, `run`, `collect`, `merge`, `modal-shell`, `init`, `fork-tool`.
(`modal-shell` is the Modal executor's workstation entry point and is meaningless elsewhere.)

`merge` — "Pool sibling tables (same parent) into one new child table" **Siblings means same parent**: it is for pooling tables that branched off one parent, not for pooling the ends of separate lineages. However, although not recommended, prosapia does allow pooling tables whose `parent_table` values differ (cousin tables). This feature is available for backward compatibility. Run `sapia merge --help` in case of doubt.

```bash
sapia new_run --label <label>        # prints the run_dir path and nothing else
sapia run     <tool> <run_dir> [-t <table>] [flags]
sapia collect <tool> <run_dir> -t <table> [-l <dir_label>]
```

**There is no verb for reading a table.** Read TSVs with shell or pandas, on your executor — see
*Reading tables* below.

## Two phases, and where output lands

`run` fans tasks out — each writes only its own files. `collect` folds those files into the table.
That split is what keeps table writes safe under parallelism and reruns idempotent.

**You never name the output table.** The driver derives it at submit time from the tool's `action`:

- **`create`** reserves a **child table**, `table<gen+1>[_<label>]`, registers it before any task
  runs, and links each new row to its parent. With no `-t`, it starts a fresh root, `table0`.
- **`update`** writes back into the table it read. **`-t` is required**; without it the run errors.

So for a `create` tool the `Output:` line that `sapia run` prints is how you learn the child table's
name — you need it for `collect -t`. **Never guess it**; `_registry.tsv` also has it.

```
Submitting 5 designs
Output:  outputs/2026…_run/table1/proteinmpnn          <- <run_dir>/<table>/<leaf>
Logs:    outputs/2026…_run/table1/proteinmpnn/proteinmpnn_logs
```

## Base run flags

Shared by every tool; each tool adds its own on top.

| Flag | Default | Meaning |
| --- | --- | --- |
| `run_dir` (positional) | — | From `sapia new_run`. Pass it exactly as printed. |
| `-t`, `--table` | — | Source table, name without `.tsv`. Required for `update`; omit on a `create` to start a root. |
| `-i`, `--input-column` | the tool's `default_input_column` | Which column feeds the tool. **Check it every time** — see the traps. |
| `-l`, `--dir-label` | `""` | Same-tool variants in one table. Leaf becomes `<tool>_<label>`. Must be repeated on `collect`. |
| `--table-label` | `""` | `create` only. Names the child table `table<gen>_<label>`; labels accumulate down generations. |
| `-f`, `--filter` | none | **Path to a `.py` file** defining `apply_filter(df) -> df`, applied to the source frame before the manifest. Not an inline expression — see the traps. |
| `--force` | off | Re-submit designs this tool already finished. |
| `-e`, `--executor` | `$SAPIA_EXECUTOR` | Leave alone unless your worker says otherwise. |
| `-C`, `--max-concurrent` | `40` | Tasks running at once. |
| `-g`, `--gpus-per-task` | `1` | **`-g 0` for CPU-only tools.** |
| `-c` / `--mem` / `-T` | tool's `RESOURCES` | CPUs / memory / wall time per task. |

Executor-specific resource flags (`--modal-gpu`, `--partitions`, `-a/--account`, …) belong to the
worker agent for that executor, not here.

### The ready set

You never subset by hand. The driver submits rows with a **present `--input-column`**,
minus those already `OK` for this leaf, unless `--force`. A `-f` filter is applied first;
the framework's own filtering applies on top of whatever it returns. Because the module
sees the whole frame before the manifest is built, it is also the supported way to run a
tool on **one named design** or a small subset — `return df[df["name"] == "design_3"]` —
which is how you buy a cheap probe before committing a batch.

The already-`OK` skip only fires for **`update`** tools, whose status column lives in the table being
read. For a **`create`** tool that column lives in the *child* table, so **every ready row is
resubmitted on a rerun.**

## Collect

```bash
sapia collect <tool> <run_dir> -t <table> [-l <dir_label>] [--force]
```

`-t` is always required and names the table the **run wrote to**: the derived child for a `create`,
the same table for an `update`. There is no `-i` — collect reuses the column the run recorded in its
`.meta.json` sidecar, so the phases cannot disagree. `-l` must match the run's. It prints
`Collected N row(s) into <table>.`; rows already `OK` are skipped unless `--force`.

- `create` stamps lineage on each new row and **fails** if a child has no real parent row.
- `update` only **warns** if it produces a row that isn't in the table, and marks a ready design with
  no output on disk as `missing`.

**`Collected N row(s)` is not proof the right work was done** — see `running-a-step` for the
independent shape check that is.

## Labels

| | `-l` / `--dir-label` | `--table-label` |
| --- | --- | --- |
| Scopes | output dir + columns *within* a table | the child *table* name |
| Phases | run **and** collect (must match) | run only |
| Tools | any | `create` only |
| Use for | same-tool variants (seeds, params, two different comparisons) on the same designs | forks of the lineage |

Collect has no `--table-label` — the label is part of the table name, so pass that name to `-t`
(e.g. `-t table1_lowT`).

## Filters

Filters are **Python modules** that define `apply_filter(df) -> df`. They are applied to the source table before the manifest is built, so they can see the whole frame and return a subset. They must live where the tasks can see them, not on your own machine. Create a `filters/` subdir in the run_dir to keep them together. A filter is **not** a pandas-style string expression — see the traps.

## Reading tables

Run these **on your executor**, wherever the run_dir actually lives (the Modal workstation, or inside
an allocation on a cluster) — never on a machine that cannot see it. Both environments have pandas.

```bash
# what tables exist, and how they relate
cat <run_dir>/_registry.tsv

# columns of a table
head -1 <run_dir>/table1.tsv | tr '\t' '\n'

# how many designs actually succeeded
python -c "import pandas as pd; d=pd.read_csv('<run_dir>/table1.tsv',sep='\t'); print(d['boltz_status'].value_counts())"

# selected columns, ordered by the one you are gating on
python -c "import pandas as pd; d=pd.read_csv('<run_dir>/table1.tsv',sep='\t'); \
print(d[['name','boltz_status','boltz_confidence_score']].sort_values('boltz_confidence_score',ascending=False).head(10).to_string())"
```

**Return aggregates and bounded slices, never whole tables.** A `value_counts()`, a count above a
threshold, or a head of ≤10 rows ordered by the gating column answers the question; a 200-row dump
spends the caller's context on rows nobody will read. Say what you withheld.

Reading a `_status` column: `OK` means that design succeeded for that leaf. `missing` means the row
was never submitted — expected whenever a filter was used, and **not a failure count**. Anything else
is an error whose text is usually the diagnosis.

## Traps

Each of these produces a *successful-looking* command that did the wrong thing.

- **`No designs to submit.` exits 0.** A wrong `-i`, an empty filter or an all-`OK` table
  prints that one line and returns success. Always read the submit output: if you don't
  see `Submitting N designs`, nothing was queued — say so instead of moving on to collect.
- **Default input columns rarely match.** `proteinmpnn`'s is `rfdiffusion_path`
  (RFdiffusion, *not* rfdiffusion3) — from an rfd3 table you must pass
  `-i rfdiffusion3_path`. The predictors default to `proteinmpnn_sequence`. Check the
  parent table's header before composing the command.
- **A `create` rerun from the same parent reuses the same child table name.** `table1` is
  derived from the parent's `gen`, not from what already exists, and re-registering an
  identical lineage is a no-op. So a second rfd3 run on `table0` writes into the same
  `table1` — and a *different* `create` tool on `table0` also derives `table1` and shares
  it. Pass `--table-label` whenever you branch, or you will merge two experiments.
- **A `--dir-label` mismatch fails at collect**, with `Output dir not found` — collect
  looks under the labelled leaf.
- **`-f` is a module *path*, not an expression.** It is declared `type=Path` and loaded
  with `importlib.util.spec_from_file_location`; the module must define
  `apply_filter(df) -> df`, or you get
  `AttributeError: … does not define an 'apply_filter' function`. Passing a pandas-style
  string fails at submit with `ImportError: Could not load module from <your expression>`
  — measured, and it costs a submit.
- **A `-f` filter module must live on the Volume.** The workstation image carries only
  prosapia's source, your tools dirs and `.env`; a local path won't exist there. Write the
  filter into the run_dir through the workstation (heredoc via `--cmd`) and pass that path.
  `modal volume put` is unreliable on this network.
- **Not every tool runs on Modal.** A tool needs a `modal_image.py`. Present for
  rfdiffusion3, proteinmpnn, boltz, alphafold3, colabfold, usalign. **Absent** for
  rfdiffusion (v1), openfold3, make_symmdef, align_symm_axis — those cannot run here; say
  so rather than trying.
- **Run `sapia` from `/runs`, never `cd` into a run_dir.** Stored paths are relative to the
  cwd, and tasks start at `/runs`; a table collected from inside a run_dir points at files
  no later task can resolve.
- **`sapia` never runs locally.** There is no `/runs` on this machine.

## Tool roster

| Tool | action | default `-i` | Modal | Notes |
| --- | --- | --- | --- | --- |
| `rfdiffusion3` | create | `pdb_path` | yes | Root run takes no `-t`. Skill: `rfdiffusion3`. |
| `proteinmpnn` | create | `rfdiffusion_path` | yes | CPU → `-g 0`. Skill: `proteinmpnn`. |
| `boltz` | update | `proteinmpnn_sequence` | yes | Skill: `boltz`. |
| `alphafold3` | update | `proteinmpnn_sequence` | yes | Needs its own image/params env. |
| `colabfold` | update | `proteinmpnn_sequence` | yes | |
| `usalign` | update | *n/a* | yes | No input column: selects on `--col-a`, plus `--col-b` (resolved up the lineage) or `--ref`. Forces its own `-g 0`. Skill: `usalign`. |
| `rfdiffusion` | create | `pdb_path` | **no** | |
| `openfold3` | update | `proteinmpnn_sequence` | **no** | |
| `make_symmdef` | update | `pdb_path` | **no** | |
| `align_symm_axis` | update | `make_symmdef_path` | **no** | |

Load the per-tool skill (`rfdiffusion3`, `proteinmpnn`, `boltz`, `usalign`) before
composing that tool's flags; this page only covers what every tool shares.

## One caveat this page's status logic inherits

A task script may catch its own errors and still exit `0` (`usalign` and `pyrosetta` do). **Exit
codes prove the tasks *ran*; the `<leaf>_status` column after collect is what proves they
*worked*.** Check it before calling a step successful.
