---
name: authoring-a-tool
description: >-
  How to add a new tool to this pipeline: the spec.py descriptor, run_<name>.py manifest builder, collect_<name>.py collector, .sh task script and optional modal_image.py, wired into the `sapia` CLI, plus the naming, column and error conventions that outlive the code. Use when asked to create, add, scaffold or write a new tool, build a manifest or collect function, wrap a binary as a pipeline step, or pick a tool `action` (create / update).
---

# Authoring a pipeline tool

A tool wraps one program as a pipeline step: it fans designs out as tasks (a SLURM array or Modal containers, chosen by `--executor`), runs the program per task, and collects results into a TSV table in a run_dir. It carries no orchestration of its own: a declarative `Tool` descriptor (`spec.py`), two pure hooks and one `.sh` task script; the shared drivers (`run_from_args` / `collect_from_args`) supply the flow and the `sapia` CLI discovers and dispatches every tool. There are no per-tool `__main__` entry points.

Read the real code first, it is the source of truth: `core/tool.py`, `base_run.py`, `base_collect.py`, `naming.py`, `base_parser.py`, `data_manager.py`, `tool_registry.py`, `core/executors/` (`slurm.py`, `modal.py`, `volume_path`), `utils/expr.py` and `utils/positions.py` (the `{expr}` language), `cli/cli.py`, `core/scripts/sapia_task_prelude.sh`, and the prose docs in a prosapia checkout (`docs/writing-a-tool.md`, `writing-a-build-manifest-function.md`, `writing-a-collect-function.md`, `running-on-modal.md`). Examples to mirror, in the bundled tools: `proteinmpnn` (create, extra args, Modal image), `make_symmdef` (minimal update, CPU only), `align_symm_axis` (a Python worker), `rfdiffusion3` (root create), `usalign` (collect flags). Find the closest existing tool and copy its structure, naming and idioms.

## Anatomy

| File | Role | Required |
|---|---|---|
| `spec.py` | Assembles `TOOL = Tool(...)`; the single wiring point the CLI discovers | always |
| `run_<name>.py` | Pure module exporting `build_<name>_manifest(ctx)` (and any args function) | always |
| `<name>.sh` | Per-task script shared by every executor: sources the prelude, reads one manifest line, runs the program, writes per-design output | always |
| `collect_<name>.py` | Pure module exporting `collect_<name>(ctx)`, a per-design collector factory | always |
| `modal_image.py` | `image()`, optional `RESOURCES` / `volumes()`, for `--executor modal` | only to run on Modal |
| `<name>_worker.py` | A Python step before the program (parsing, format conversion, a subprocess chain) | only when needed |

Do not add a worker for a simple shell step: put a straightforward command directly in the `.sh`. Discovery loads every `tools/*/spec.py`; built-ins come first, then each dir on `$PROSAPIA_TOOLS_DIR`, later dirs shadowing earlier ones by `Tool.name`. To start from a built-in use `sapia fork-tool`; to change a field or two use `with_overrides` (see `editing-a-tool`). **To change a bundled tool in this workspace, fork it under `tools/` and keep its `name`**: the columns keep their names and downstream steps are unaffected.

## Pick the `action`

- **`create`** reserves a new table that `collect` fills. With `--table` it is a **child** (every collected row must carry a parent row via `Collected.parent`; the driver stamps `parent_table` and `gen`); without `--table` it is a **root** (`parent_name = NA`, gen 0; the builder must record the group names with `ctx.write_meta(root_designs=[...])`, the only thing collect can key off).
- **`update`** annotates the table it ran against in place: no new rows, `--table` required. Use it when the tool derives properties of designs that exist.

The driver enforces the contract: a create collect without a valid parent raises; an update collect that invents rows warns.

## `spec.py`

```python
TOOL = Tool(
    name="make_symmdef",
    action="update",                                   # or "create"
    description="Make symmetry definition files.",     # shown in run/collect --help
    default_script=str(Path(__file__).parent / "make_symmdef.sh"),
    default_input_column="pdb_path",                   # overridable with --input-column
    build_manifest_fn=build_make_symmdef_manifest,
    collect_fn=collect_make_symmdef,
    # add_run_args_fn=..., add_collect_args_fn=...     # only when the tool has its own flags
)
```

Collect reads back whichever column the run used from the `.meta.json` sidecar, so collect takes no `--input-column`.

## The run module

`build_<name>_manifest(ctx) -> list[tuple[str, ...]]` returns one tab-separated row per task; the driver handles table resolution, output dirs, `--filter`, the manifest and submission, and the builder never submits or reads driver state back. `ManifestCtx` gives `df` (source frame, already filtered; empty for a root create), `args`, `out_dir`, `lookup(name, column)` (read-only lineage walk), `write_meta`, and `ready`.

```python
def build_make_symmdef_manifest(ctx):
    ctx.args.gpus_per_task = 0                  # CPU-only tool: no GPU request
    rows = []
    for name in ctx.ready.index:                # designs to submit
        src = Path(str(ctx.ready.at[name, ctx.args.input_column]))
        rows.append((name, str(ensure_pdb(src, ctx.args.run_dir) if src.exists() else src)))
    return rows
```

- **`ctx.ready`** is rows with a present `--input-column` minus those this leaf already finished (`<leaf>_status == "OK"`) unless `--force`; the skip and `--force` are framework behaviour, never re-implement them. The skip fires only for `update` tools (the status column lives in the table read); a `create` tool resubmits every ready row.
- **CIF to PDB** up front with `ensure_pdb(src, run_dir)` (cached under `run_dir/.cif_to_pdb`). **Modal paths**: wrap every path the task reads in `volume_path(...)`; off Modal it is `abspath`.
- **A manifest row is a task, not necessarily a design**: a builder may pack several designs per task (one sub-manifest per task, a 1-field row pointing at it).
- Extra flags: an `add_run_<name>_args(parser)` wired via `add_run_args_fn`, plus a `CommonArgs` subclass used only to type `ManifestCtx[...]`. The base parser already supplies `run_dir`, `-t`, `--force`, `-i`, `-l`, `-f`, `-w`, `-s`, `-e`, `--modal-gpu`, the resource flags, the SLURM flags and (create) `--table-label`; do not redeclare them.
- `run_dir` is positional and must already exist (`sapia new_run`).

### `{expr}`: per-design values inside a flag

A flag is one string for the whole run but designs differ; `{expr}` islands close the gap. Outside the braces is the wrapped program's native language and prosapia never parses it; inside is prosapia's. Accept the program's own syntax and pass it through `resolve_template(flag, ctx.lookup, name)` (`'A1-{motif_end},30'` becomes `'A1-131,30'`); `parse_positions(spec, lookup, name, lengths=…)` adds `/` chain breaks, `,` fragments and `start:end` ranges. It is opt-in per flag: nothing resolves unless the builder calls it.

- Grammar: integer literals, bare column names, `+ - * //` (and unary sign); `/` is not an operator. Integers only: `{boltz_plddt}` raises `not an integer`. A column resolves up the lineage, nearest non-null wins.
- Resolve at manifest-build time, never in the task script: the manifest carries literal values, and a per-row failure raises before anything is submitted. Decide deliberately whether an unresolvable column aborts the run or skips that row (a silent skip reads as a filter problem).
- A root create has no lineage, so `{expr}` cannot work there. Braces in a native flag (a YAML or JSON flow mapping) are read as islands and fail to parse. Positions are 1-indexed per chain, order preserved, not de-duplicated; open-ended ranges need `lengths`. For a structured flag (a YAML or JSON file) walk the tree and resolve strings and dict keys (`_resolve_tree` in `run_rfdiffusion3.py`).

## The `.sh` task script

Source the prelude, which loads `.env` and exports `MANIFEST` (`$1`), `OUT_DIR` (`$2`), `SAPIA_SCHEDULER` (`slurm` | `modal`), `SAPIA_TASK_ID` (1-based) and `SAPIA_LINE` (this task's tab-separated manifest line), then call `sapia_activate SAPIA_ACTIVATE_<NAME>` (a no-op under Modal). Use `SAPIA_TASK_ID`, never `SLURM_ARRAY_TASK_ID`; address sibling files as `"${SAPIA_TOOL_DIR:?}/<name>_worker.py"`; address interpreters through an overridable variable (`PY=${PIPELINE_PYTHON:-python}`), never bare `python`.

```bash
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --job-name=make_symmdef
set -euo pipefail
source "${SAPIA_PRELUDE:?}"
sapia_activate SAPIA_ACTIVATE_SYMMDEF
NAME=$(echo "$SAPIA_LINE" | cut -f1)
SRC=$(echo "$SAPIA_LINE" | cut -f2)
# run the program on $SRC, write a per-design result the collector can read
```

Write a per-design result the collector can find (program output in a known layout under `$OUT_DIR/<name>/`, or a small `$OUT_DIR/<name>.tsv` with a `status` column), and record errors as data (`status` of `error: ...`) rather than only failing the task, so partial runs still collect. `#SBATCH` sets resources; `-c/-T/--mem` override them. **Never name a local after a bash special variable** (`GROUPS`, `UID`, `EUID`, `PPID`, `PIPESTATUS`, `SECONDS`, `RANDOM`, `LINENO`, `BASH_*`, `FUNCNAME`, `HOSTNAME`, `IFS`, `OPTARG`, `OPTIND`, `REPLY`, `SHLVL`, `PWD`, `OLDPWD`, `PATH`, `HOME`): `GROUPS=$(…)` fails with rc 1, the assignment is discarded, and `set -euo pipefail` kills the task with `.out` and `.err` both empty and no diagnostic anywhere. Prefix manifest locals (`SEL_GROUPS`) or test with `bash -c 'declare -p NAME'`; debug by re-running the script under `bash -x` with the real task environment.

## The Modal image (optional)

A tool runs under `--executor modal` only if `modal_image.py` sits next to its `.sh`: `image()` (required, returns a `modal.Image` providing everything `sapia_activate` would have sourced, with the env vars the `.sh` expects), optional `RESOURCES = {"gpu": "L4", "cpu": 8, "memory": "8G", "timeout": "01:00:00"}` (overridden by `--modal-gpu` and the resource flags) and `volumes()` for weights and caches. The `.sh`, its siblings and the prelude are added to the image automatically; an image with no Python needs `add_python=...`.

## The collect module

`collect_<name>(ctx)` is a per-design collector factory: it runs once per collect (do any directory scan here), then returns `one(design) -> Iterable[Collected]`. The driver iterates `ctx.ready`, calls `one` per design, and stamps `<leaf>_status`, `<leaf>_path` and `parent_name` from every `Collected`, so a tool never writes those names itself.

```python
def collect_make_symmdef(ctx):
    def one(d):
        tsv = ctx.out_dir / f"{d.name}.tsv"
        if not tsv.is_file():
            yield Collected(status="missing", path="")
            return
        row = pd.read_csv(tsv, sep="\t").iloc[0]
        yield Collected(status=str(row["status"]), path=str(row["symm_path"]))
    return one
```

`Collected` carries `data` (bare column names; the driver prefixes `<leaf>_`), `path`, `status` (default `"OK"`), `name` (create tools minting child rows) and `parent`. Yield one for an update, several (with `name`/`parent`) for a create; yield nothing when a design has no output (an update marks the row `missing`, a create skips it). A tool whose status is per comparison carries each as its own `data` column and sets `status=None`. `CollectCtx` exposes `df`, `args`, `table_name`, `out_dir`, `parent_table`, `parent_df`, `lookup`, `creates_table`, `default_input_column`, `ready`. For a create, set `parent` to a row present in `parent_df` and do not propagate ancestor values (leave them to `lookup`); for an update, a new name warns that the tool should probably be a create. Optional collect flags follow the run-flag pattern (`add_collect_<name>_args`, a `CollectArgs` subclass).

## Conventions that outlive the code

- **Name columns with one convention** (one spelling, one position for qualifiers such as a nitrogen or a side, one prefix style) and document every irregularity in the tool's skill: names become permanent once a table exists, because renaming strands old rows, filters and downstream `-i` flags.
- **A hard-coded column list lives in several places**: the worker, the collector's typed lists, and the error-row list. Update all of them together; the error row must carry every science column as `NA`, never `0`.
- **Additive blocks are non-fatal**: a block added to a working tool (a second model, a stability term) must fail into an `<block>_error` column with the row's status left `OK` and its columns `NA`, and collect must print one count line, including zero (`N of M OK rows carry a <block>_error`), so an `NA` is not read as a zero.
- **State the sign and direction of every column** (higher is better, negative is wanted, units) in the docstring and in the skill, verified against the source, not recalled; a mirrored sign convention is the classic silent error. Say what each column is blind to and which tool answers that.
- **Trust columns exist** for anything that maps, aligns or renumbers; emit them on error rows too (they diagnose the failure), and make the premise failure an error row.
- **Errors are data, and NA is not 0**: a deliberately bad input must produce an error status (see `verifying-a-tool`).

## Checklist

1. Create the tool folder; write `run_<name>.py` with `ctx.ready`, optional args, no `__main__`; pass any residue, length or index flag through `resolve_template`.
2. Write `collect_<name>.py` yielding `Collected(...)` and honouring the create/update contract.
3. Write `<name>.sh` (prelude, `sapia_activate`, cut fields, per-design output, `#SBATCH` resources, `gpus_per_task = 0` in the builder for CPU-only tools); add `modal_image.py` if it should run on Modal (wrap manifest paths in `volume_path`); add a worker only if the step needs real Python.
4. Assemble `spec.py`; match neighbours for imports, docstring style and error-as-data reporting.
5. Verify discovery and flags (`sapia run <name> -h`, `sapia collect <name> -h`), sanity-check manifest and collect locally, then run the checklist in `verifying-a-tool`.
6. Write the skill per `writing-a-tool-skill`, and tell the user to commit the tool.
