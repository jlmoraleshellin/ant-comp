---
name: authoring-a-tool
description: >-
  How to add a new tool to this pipeline: the spec.py descriptor, run_<name>.py manifest builder, collect_<name>.py collector, .sh task script and optional modal_image.py, wired into the `sapia` CLI. Use when asked to create, add, scaffold or write a new tool, build a manifest or collect function, wrap a binary as a pipeline step, or pick a tool `action` (create / update).
---

# Authoring a pipeline tool

A **tool** wraps one protein-design program as a pipeline step: it fans designs out as tasks (a SLURM array or Modal containers, chosen by `--executor`), runs the program per task, and collects the results back into a TSV table inside a run_dir. This skill is the contract for writing one so it matches the existing tools exactly.

A tool carries **no orchestration of its own**. It is a declarative `Tool` descriptor (in `spec.py`) plus two pure hooks and one `.sh` task script; the shared drivers (`run_from_args` / `collect_from_args`) supply the flow, and the `sapia` CLI discovers every tool and dispatches to it. There are **no per-tool `__main__` entry points** — the run/collect files are plain modules that export functions.

**Before writing anything, read the real code** — it is the source of truth, not this file:

- Drivers & types: `src/prosapia/core/tool.py`, `base_run.py`, `base_collect.py`, `naming.py`, `base_parser.py`, `data_manager.py`, `tool_registry.py`.
- Executors: `src/prosapia/core/executors/` (`slurm.py`, `modal.py`, `volume_path`).
- The `{expr}` mini-language (no docs page exists): `src/prosapia/utils/expr.py`, `src/prosapia/utils/positions.py`.
- CLI dispatcher: `src/prosapia/cli/cli.py` (builds `sapia {run,collect} <tool>`).
- Task prelude sourced by every `.sh` task script: `src/prosapia/core/scripts/sapia_task_prelude.sh`.
- User docs: `docs/writing-a-tool.md`, `docs/writing-a-build-manifest-function.md`, `docs/writing-a-collect-function.md`, `docs/running-on-modal.md`.
- Full example (create + extra args + Modal image): `src/prosapia/tools/proteinmpnn/` — `spec.py`, `run_proteinmpnn.py`, `proteinmpnn.sh`, `collect_proteinmpnn.py`, `modal_image.py`.
- Minimal example (update, no worker, CPU-only): `src/prosapia/tools/make_symmdef/`.
- Worker example (a Python pre-step): `src/prosapia/tools/align_symm_axis/`.
- Root-create example (no `--table`): `src/prosapia/tools/rfdiffusion3/`.
- Collect-flags example: `src/prosapia/tools/usalign/`.

Then find the **closest existing tool** to what's being asked and mirror its structure, naming, and idioms. Prefer copying a neighbour over inventing.

## Anatomy of a tool

A tool lives in `src/prosapia/tools/<name>/` and consists of:

| File | Role | Required? |
|---|---|---|
| `spec.py` | Assembles the `Tool(...)` descriptor as a module global `TOOL`, importing the two hooks from the sibling modules. This is the single wiring point; the CLI discovers it. | Always |
| `run_<name>.py` | Pure module. Exports `build_<name>_manifest(ctx)` (and any extra-args fn / args subclass). No `__main__`. | Always |
| `<name>.sh` | The per-task script, shared by every executor. Sources the prelude, reads one manifest line, runs the program, writes per-design output. | Always |
| `collect_<name>.py` | Pure module. Exports `collect_<name>(ctx)` — a per-design collector factory (see below). No `__main__`. | Always |
| `modal_image.py` | The tool's Modal image (`image()`, optional `RESOURCES` / `volumes()`), used by `--executor modal`. | Only to run on Modal |
| `<name>_worker.py` | A thin Python step **before** the program when extra pre-processing is needed. | Only when needed |

Tools are discovered by `tool_registry.discover`, which loads every `tools/*/spec.py`. Built-ins live under `src/prosapia/tools/`; user tool dirs can be added via `$PROSAPIA_TOOLS_DIR` (later dirs shadow built-ins). To start from a built-in, `sapia fork-tool <name>` copies it into your tools dir; to change only a field or two, reuse the built-in via `get_builtin(<name>).with_overrides(...)` in your own `spec.py`.

**Do not add a worker for a simple shell step.** If the per-design work is a straightforward command, put it directly in the `.sh` task script — see `make_symmdef.sh`. Add a `_worker.py` only when the per-design step needs real Python: structure parsing, format conversion, a multi-step subprocess chain (e.g. `align_symm_axis_worker.py` uses gemmi + numpy).

## Step 0 — pick the `action`

`Tool.action` sets the table contract (`src/prosapia/core/tool.py`). There are **two** actions:

- **`create`** — reserves a new table that `collect` fills. Covers two cases, distinguished only by whether `--table` is given at run time:
  - **child** (with `--table`): spawns a child table from the parent. Every collected row MUST carry a parent (`Collected.parent` → `parent_name`, a row in the parent table); the driver stamps `parent_table` + `gen`.
  - **root** (no `--table`): spawns a root table with no parent; the driver stamps `parent_name = NA` and `gen 0`. Use for a tool that starts a fresh lineage (from scratch, or from a direct input like `--input-pdb`). `--table` is optional for a create tool.
- **`update`** — annotates the table it ran against, in place. No new rows; adds columns to existing rows. Use when the tool measures/derives properties of designs that already exist (a metric, a relaxed structure path). `--table` is required.

The driver enforces the contract: a `create` collect whose rows lack a valid parent on a child table raises; an `update` collect that invents new rows warns.

## The `spec.py` descriptor

`spec.py` binds the whole tool together and exposes it as `TOOL`. Mirror a neighbour — `make_symmdef/spec.py` (minimal) or `proteinmpnn/spec.py` / `usalign/spec.py` (with args):

```python
from pathlib import Path

from prosapia.core import Tool

from .collect_make_symmdef import collect_make_symmdef
from .run_make_symmdef import build_make_symmdef_manifest

TOOL = Tool(
    name="make_symmdef",
    action="update",                                   # or "create"
    description="Make symmetry definition files.",     # shown in run/collect --help
    default_script=str(Path(__file__).parent / "make_symmdef.sh"),
    default_input_column="pdb_path",                   # table column the run reads
    build_manifest_fn=build_make_symmdef_manifest,
    collect_fn=collect_make_symmdef,
    # Optional, when the tool adds its own flags:
    # add_run_args_fn=add_run_proteinmpnn_args,
    # add_collect_args_fn=add_collect_usalign_args,
)
```

`default_input_column` is the table column the run iterates over by default (overridable at run time with `--input-column`); `collect` reads back whichever column the run actually used, via the sidecar (see below). `Tool.metadata` derives the `ToolMetadata` (name, action, description, default_input_column) the drivers consume.

## The run module (`run_<name>.py`)

Export a `build_<name>_manifest(ctx) -> list[tuple[str, ...]]` — one tab-separated row per task. The driver handles table resolution, output-dir creation, `--filter`, writing the manifest, and submission via the chosen executor. The builder only returns rows; it never reads driver state back or submits anything itself.

```python
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.utils import ensure_pdb


def build_make_symmdef_manifest(ctx: ManifestCtx[CommonArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0          # CPU-only tool → no GPU request

    ready = ctx.ready                   # designs to submit (see below)
    manifest_rows: list[tuple[str, ...]] = []
    for name in ready.index:
        name = cast(str, name)
        src = Path(str(ready.at[name, ctx.args.input_column]))
        # CIF→PDB up front, cached under run_dir/.cif_to_pdb:
        input_pdb = ensure_pdb(src, ctx.args.run_dir) if src.exists() else src
        manifest_rows.append((name, str(input_pdb)))
    return manifest_rows
```

`ManifestCtx` gives the builder: `df` (source table frame, already filtered; empty for a root create), `args`, `out_dir` (the resolved output dir), `lookup(name, column)` (read-only lineage walk up parent tables — never mutate through it), `write_meta` (record extra run params into the sidecar — needed by root creates, see below), and the `ready` property.

- **`ctx.ready`** is the set of designs this run should submit: rows with a present `--input-column`, minus those this tool already finished (`<leaf>_status == "OK"`) unless `--force`. Iterate it instead of hand-filtering. The already-OK skip (resume-on-rerun) and `--force` are framework behaviour — don't re-implement them per tool.
- **Root create (no `--table`)**: `ctx.df`/`ctx.ready` are empty, so the builder synthesizes its design group(s) from `args` instead. It **must** record their names with `ctx.write_meta(root_designs=[...])` — a root run has no parent table for the collect phase to iterate, so this sidecar record is the only thing collect can key off. (`write_meta` merges into the sidecar; it won't clobber base run params.)
- **CIF→PDB**: if the program needs a PDB, convert up front with `ensure_pdb(src, args.run_dir)` (from `prosapia.utils`).
- **Paths for Modal**: any path written into the manifest must resolve inside the task container. Wrap paths the task reads with `volume_path(...)` (from `prosapia.core.executors`) so they are in the `/runs` mount's form, as `run_proteinmpnn.py` does; off Modal it is just `abspath`.
- **Batching**: a manifest row is a task, not necessarily a design. A tool may pack several designs per task (e.g. proteinmpnn writes one sub-manifest per task and emits a 1-field row pointing at it).
- `run_dir` is a positional arg and must already exist — mint one with `sapia new_run --label <label>`; the driver raises if it is missing.

### Extra CLI flags

If the tool needs its own flags, add an `add_run_<name>_args(parser)` function wired via `add_run_args_fn` in `spec.py`, plus a `CommonArgs` subclass used only as the type parameter of `ManifestCtx[...]` for type-checking (see `proteinmpnn`):

```python
from argparse import ArgumentParser

from prosapia.core import CommonArgs, ManifestCtx


class ProteinMPNNArgs(CommonArgs):
    num_seq_per_target: int
    # ...


def build_proteinmpnn_manifest(ctx: ManifestCtx[ProteinMPNNArgs]) -> list[tuple[str, ...]]:
    ...


def add_run_proteinmpnn_args(parser: ArgumentParser) -> None:
    parser.add_argument("--num-seq-per-target", type=int, default=2)
```

`run_dir`, `-t/--table`, `--force`, `-i/--input-column`, `-l/--dir-label`, `-f/--filter`, `-s/--script`, `-e/--executor`, `--modal-gpu`, the resource flags (`-g/--gpus-per-task`, `-c/--cpus-per-task`, `-T/--time`, `--mem`, `-C/--max-concurrent`), the SLURM flags (`-a/--account`, `-p/--partitions`, `--max-gpu-fraction`), and (for create tools) `--table-label` are all supplied by the base parser — don't redeclare them.

### The `{expr}` mini-language — per-design values inside a flag

A flag is **one string for the whole run**, but designs differ: the motif ends at a different residue in every row. The mini-language closes that gap, so a user writes one flag that means something different per design instead of splitting the run or writing a filter per value.

Source (there is **no `docs/` page** for this — the modules are the reference): `src/prosapia/utils/expr.py` and `src/prosapia/utils/positions.py`.

**It is opt-in per flag.** Nothing resolves unless your builder calls the resolver. Two layers, both exported from `prosapia.utils`:

| Layer | Call | Returns | Grammar |
| --- | --- | --- | --- |
| Expression | `resolve_expr(expr, lookup, name)` | `int` | One expression, no braces |
| Template | `resolve_template(s, lookup, name)` | `str` | Substitutes every `{expr}` island in `s`; text outside braces is **left verbatim** |
| Positions | `parse_positions(spec, lookup, name, lengths=…)` | `list[list[int]]` | Template first, then `/` chain breaks, `,` fragments, `start:end` inclusive, optional outer `[…]` |

#### The division of labour — this is the whole idea

> **Outside the braces is your tool's native language, and prosapia never parses it. Inside the braces is prosapia's.**

RFdiffusion contigs, ProteinMPNN position lists, a YAML spec field — prosapia passes them through untouched and only substitutes the islands. So **never invent a syntax for "take this number from the table."** Accept the wrapped program's own syntax verbatim and run it through `resolve_template`:

```python
from prosapia.utils import resolve_template

for name in ctx.ready.index:
    contig = resolve_template(ctx.args.contigs, ctx.lookup, name)
    # 'A1-{motif_end},30'  ->  'A1-131,30' for this design
```

`ctx.lookup` is already on `ManifestCtx`; that call is the entire integration.

#### What the grammar allows

Integer literals, **bare column names**, and `+ - * //`. Nothing else — `_safe_eval` walks the AST and rejects every other node type.

```
{motif_end}                 # column, resolved up the lineage for this design
{motif_end - 1}             # arithmetic on it
{binder_len // 2}
{n_res * 2 + 10}
```

- **`/` is not an operator.** `ast.Div` is absent from `_BIN_OPS`, so `{a/b}` raises. Use `//`. (True division would not return an int, and `/` is already the chain separator in the position language.)
- **Integers only.** A column resolving to a non-integral float raises `column 'x'=0.93 is not an integer` — `{boltz_plddt}` will never work. This is a geometry/indexing language, not a general templating engine.
- **Lineage walk, nearest non-null wins.** `lookup` reads the column on this row, and on a miss or a null follows `parent_name`/`parent_table` upward (with a cycle guard). A child's own value shadows its parent's.

#### Resolve at manifest-build time, never in the task script

The manifest carries **already-resolved** values. The task script gets a literal contig and a literal position list; it has no table, no lineage and no resolver. A per-row failure therefore raises during the build, **before anything is submitted** — the cheap place to fail. Let it.

#### Structured flags: walk the tree

For a flag that is a YAML/JSON file rather than a string, resolve recursively — strings and dict *keys* go through `resolve_template`, containers are walked, and `int`/`float`/`bool`/`None` pass through so native types survive. Copy `_resolve_tree` in `run_rfdiffusion3.py:249`; that is what lets `--extra-spec` carry `{expr}` in its values.

#### Traps

- **A root `create` has no lineage.** `ctx.df` is empty, there is no row to resolve against, so `{expr}` cannot work. Root runs must use literal values — the rfd3 skill states this as "`{expr}` placeholders … need `-t`".
- **Braces collide with native syntax that uses them.** Any YAML/JSON flow mapping in a flag value will be read as an island and fail to parse as an expression. Avoid literal braces in `--extra-spec`, or your users will hit a confusing error.
- **An unresolvable column is a per-row `ValueError` naming the design.** Decide deliberately whether your builder lets it abort the run or skips that row — silently skipping means a user sees a smaller `Submitting N designs` than they expected, which reads like a filter problem rather than a typo.
- **The position parser's error is migration-friendly on purpose** — a bare column name outside braces reports `wrap table column expressions in braces, e.g. '{motif_end}'` rather than a raw `int()` failure. Match that tone if you write your own parser.
- **Positions are 1-indexed per chain, order preserved, not de-duplicated.** Open-ended ranges (`10:`, `:50`, `:`) need a `lengths` argument; without one they raise. ProteinMPNN never passes lengths, so open ranges simply don't work there.

## The `.sh` task script

Source the shared prelude, which loads `.env` and exports `MANIFEST` (`$1`), `OUT_DIR` (`$2`), `SAPIA_SCHEDULER` (`slurm` | `modal`), `SAPIA_TASK_ID` (this task's 1-based index) and `SAPIA_LINE` (this task's manifest line — cut your own fields from it), then call `sapia_activate SAPIA_ACTIVATE_<NAME>` (a no-op under Modal, where the image provides the environment). Manifests are **tab-separated**. The same script runs under every executor, so use `SAPIA_TASK_ID`, never `SLURM_ARRAY_TASK_ID`. The executor also exports `SAPIA_TOOL_DIR` (the dir holding the `.sh`), so address sibling files like a worker as `"${SAPIA_TOOL_DIR:?}/<name>_worker.py"`.

```bash
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --job-name=make_symmdef

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_SYMMDEF

NAME=$(echo "$SAPIA_LINE" | cut -f1)
SRC=$(echo "$SAPIA_LINE" | cut -f2)

# ... run the program on $SRC, write a per-design result the collector can read.
```

Write a **per-design result the collector can find** — either program output in a known layout under `$OUT_DIR/<name>/`, or a small per-design TSV (`$OUT_DIR/<name>.tsv` with a `status` column). Record errors as data (a `status` of `error: ...`) rather than only failing the task, so partial runs still collect.

Set `#SBATCH` resources for the tool (GPU vs CPU — and set `ctx.args.gpus_per_task = 0` in the builder for CPU-only tools); the `-c/-T/--mem` flags override them at run time. If a step needs a specific interpreter (a gemmi worker, the program's own Python), address it through an overridable variable rather than bare `python` — e.g. `align_symm_axis.sh` uses `PY=${PIPELINE_PYTHON:-python}`, `proteinmpnn.sh` uses `PROTEIN_MPNN_PYTHON=${PROTEIN_MPNN_PYTHON:-python}`.

## The Modal image (`modal_image.py`, optional)

A tool runs under `--executor modal` only if it ships a `modal_image.py` next to its `.sh` (see `docs/running-on-modal.md` and `proteinmpnn/modal_image.py`):

```python
import modal

RESOURCES = {"gpu": "L4", "cpu": 8, "memory": "8G", "timeout": "01:00:00"}  # optional


def image() -> modal.Image:                              # required
    return (
        modal.Image.debian_slim(python_version="3.12")
        .pip_install("torch", "numpy")
        .env({"PROTEIN_MPNN": "/opt/ProteinMPNN"})
    )


# def volumes() -> dict[str, modal.Volume]: ...          # optional: weights, DBs, caches
```

The image provides everything `sapia_activate` would have sourced on SLURM, and should export the same env vars the `.sh` expects. Its Python version need not match the workstation's; an image with no Python (e.g. `Image.from_registry`) needs `add_python=...`. The `.sh`, its sibling files and the prelude are added to the image automatically. `RESOURCES` is overridden by `--modal-gpu` and the resource flags.

## The collect module (`collect_<name>.py`)

Export `collect_<name>(ctx)` — a **per-design collector factory**. It runs once per collect (do any directory scan here), then returns a per-design function `one(design) -> Iterable[Collected]`. The driver iterates the ready designs, calls `one` for each, and stamps `<leaf>_status` / `<leaf>_path` / `parent_name` from every `Collected` you yield — so a tool never touches those column names itself.

Full guide with worked update + create examples: `docs/writing-a-collect-function.md`.

```python
from typing import Iterable

import pandas as pd

from prosapia.core import Collected, CollectCtx, CollectEach, DesignCtx


def collect_make_symmdef(ctx: CollectCtx) -> CollectEach:
    # (setup once — e.g. scan ctx.out_dir into an index — then return `one`)
    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="")
            return
        row = pd.read_csv(tsv_path, sep="\t").iloc[0]
        yield Collected(status=str(row["status"]), path=str(row["symm_path"]))

    return one
```

`Collected` carries a tool's output for one row: `data` (tool-specific columns, as **bare names** — the driver prefixes each key to `<leaf>_<key>` so same-tool variants never collide), `path` (→ `<leaf>_path`), `status` (→ `<leaf>_status`, default `"OK"`), `name` (override the row key — create tools minting child rows), `parent` (→ `parent_name`). Yield **one** `Collected` for an update, **several** (with `name`/`parent`) for a create. Yield **nothing** when a design has no output: an update marks the row `missing`; a create skips it. For the rare tool whose status is per-comparison (one out_dir hosting several named comparisons), carry each status as its own `data` column and set `status=None` to suppress the leaf stamp.

`DesignCtx` (the per-design handle) gives `name`, `out_dir`, `lookup`. The setup factory receives the full `CollectCtx` — `df`, `args`, `table_name`, `out_dir`, `parent_table`, `parent_df`, `lookup`, `creates_table`, `default_input_column`, and `ready` — so capture whatever the closure needs (indices, `args.*`). `status_col` / `path_col` still exist on `CollectCtx` but a tool that yields `Collected` never needs them.

- **`ctx.ready`** are the designs to collect: for a create tool, the parent table's ready rows (or, for a **root** create, the group names the run recorded in the sidecar via `ctx.write_meta(root_designs=…)` — there is no parent table); for an update tool, this table's ready rows minus already-collected (`OK`) ones unless `--force`. It reads the input column the *run* recorded in the out_dir sidecar (`.meta.json`), so run and collect can't disagree — which is why **collect takes no `--input-column`**.
- **`create`**: set `Collected.parent` to a row present in `parent_df` (child table). The driver validates the edge and stamps `parent_table`/`gen`; don't set those. Don't propagate ancestor values — leave them to `lookup` later. (Root create: no parent to set.)
- **`update`**: yield rows for designs already in `df`. A new name (via `Collected.name`) triggers a warning (it means the tool should probably be `create`).

Optional collect flags follow the same pattern as run flags: an `add_collect_<name>_args(parser)` fn wired via `add_collect_args_fn` in `spec.py`, plus a `CollectArgs` subclass for typing `CollectCtx[...]` (see `usalign`). `run_dir`, `-t/--table`, `-l/--dir-label` and `--force` come from the base parser.

## Running a tool

Everything dispatches through the `sapia` console script (no per-tool entry points):

```
sapia new_run --label <label>                              # mint a run_dir
sapia run <name> <run_dir> --table <table> [flags]         # submit the tasks (-e slurm|modal)
sapia collect <name> <run_dir> --table <out_table>         # collect outputs into <out_table>
```

The run reserves the table (`create`) or targets it (`update`); collect only fills it. For an update, pass the same table to both. For a create, collect takes the *new* table the run reserved (named `table<gen>[_<parent_label>][_<table_label>]`; see `docs/using-labels.md`). A root create takes no `--table` at run time. On Modal, run both from `sapia modal-shell` with the run_dir under `/runs`.

## Checklist for a new tool `<name>`

1. `mkdir src/prosapia/tools/<name>/`.
2. `run_<name>.py`: `build_<name>_manifest(ctx)` using `ctx.ready`; optional args subclass + `add_run_<name>_args`. No `__main__`. For any flag naming a residue, length or index, pass it through `resolve_template(s, ctx.lookup, name)` so `{expr}` works — accept the wrapped program's native syntax, don't invent your own.
3. `collect_<name>.py`: `collect_<name>(ctx)` returning a per-design `one(design)` that yields `Collected(...)`; honour the `create`/`update` contract. No `__main__`. See `docs/writing-a-collect-function.md`.
4. `<name>.sh`: `source "${SAPIA_PRELUDE:?}"`, `sapia_activate SAPIA_ACTIVATE_<NAME>`, cut fields from `$SAPIA_LINE`, run the program, write per-design output; set `#SBATCH` resources (and set `ctx.args.gpus_per_task = 0` in the builder for CPU-only tools).
5. `modal_image.py` if the tool should run on Modal; wrap manifest paths with `volume_path`.
6. `<name>_worker.py` **only** if the per-design step needs real Python; call it via `$SAPIA_TOOL_DIR`.
7. `spec.py`: assemble `TOOL = Tool(name, action, description, default_script, default_input_column, build_manifest_fn, collect_fn, [add_run_args_fn, add_collect_args_fn])`.
8. Match the neighbours for imports (`from prosapia.core import ...`, `from prosapia.utils import ...`), docstring style, and error-as-data reporting.
9. Verify discovery and flags: `uv run sapia run <name> -h` and `uv run sapia collect <name> -h`. Sanity-check manifest/collect logic locally with `uv run`. Full runs execute on the HPC/SLURM cluster or on Modal.

## Trap: never name a task-script variable after a bash special variable

Cost four silent task failures on 2026-09-28.

A `<tool>.sh` typically unpacks its manifest line into locals:

```bash
NAME=$(echo "$SAPIA_LINE" | cut -f1)
GROUPS=$(echo "$SAPIA_LINE" | cut -f3)   # <-- BROKEN
```

`GROUPS` is **pre-set by bash** as an indexed array of the user's group IDs. Assigning a
command substitution to it **fails with rc=1 and the assignment is silently discarded**:

```
$ bash -c 'GROUPS=$(echo hello); echo rc=$?; echo val=[$GROUPS]'
rc=1
val=[0]
```

Task scripts run under `set -euo pipefail`, so that rc=1 **kills the task instantly** —
before any `echo`, before the worker is ever invoked. The result is the worst possible
failure signature:

- `.exit` is non-zero (so it is not mistaken for success), but
- **`.out` and `.err` are both 0 bytes**, and
- `modal app logs <app_id>` returns nothing.

There is no diagnostic anywhere. The only way to find it is to re-run the task script by
hand with `bash -x` under the real task environment.

**Avoid these names for locals** (non-exhaustive): `GROUPS`, `UID`, `EUID`, `PPID`,
`PIPESTATUS`, `SECONDS`, `RANDOM`, `LINENO`, `BASH_*`, `FUNCNAME`, `HOSTNAME`, `IFS`,
`OPTARG`, `OPTIND`, `REPLY`, `SHLVL`, `PWD`, `OLDPWD`, `PATH`, `HOME`.

Prefix manifest-field locals (`SEL_GROUPS`, `TOOL_UID`) or check first with
`bash -c 'declare -p NAME' 2>/dev/null` — if it prints, pick another name.

**Debugging rule this implies:** an `.exit` that is non-zero with *empty* `.out` **and**
`.err` means the script died before its first statement produced output — look at variable
assignments and the prelude, not at the worker.
