---
name: editing-a-tool
description: >-
  How to change an existing (bundled or custom) tool: picking the lightest mechanism that works (run flags, a `--filter` module, a custom `--script`, the activation script, `with_overrides`, `sapia fork-tool`, or editing in place) and the naming/discovery rules that decide whether a change shadows the original or registers alongside it. Use when asked to edit, customize, override, fork or shadow a tool, or to change how it is activated. For a brand-new tool use authoring-a-tool.
---

# Editing a pipeline tool

The bundled tools are intentionally general, so most workflows need to bend one at some point. There are several ways to do it, from lightest to heaviest; **always pick the lightest one that does the job**, because every step up the ladder is more code the user now owns and must keep in sync with upstream. The anatomy of a tool (`spec.py`, `run_<name>.py`, `collect_<name>.py`, `<name>.sh`, optional `modal_image.py` / `_worker.py`) and the hook contracts are described in the authoring-a-tool skill; this skill is about changing one that already exists.

**Before changing anything, read the real code and docs** — they are the source of truth, not this file:

- The tool being edited: `src/prosapia/tools/<name>/` (all of it — `spec.py`, the run/collect modules, the `.sh`, `modal_image.py`), plus its user doc under `docs/tools/<name>.md` if one exists.
- `docs/writing-a-tool.md` ("Customizing bundled tools", "How tools are discovered"), `docs/configuration.md` (activation binding, tool discovery), `docs/running-a-tool.md` (run flags, `-s/--script`), `docs/writing-a-filter-function.md`, `docs/using-labels.md`.
- Mechanics: `src/prosapia/core/tool.py` (`Tool`, `ToolOverrides`, `with_overrides`), `src/prosapia/core/tool_registry.py` (`discover`, `get_builtin`), `src/prosapia/cli/fork_tool.py`, `src/prosapia/core/naming.py` (output dir / column leaf).

## Step 0 — find the lightest mechanism

Work out *what* has to change, then walk down this ladder and stop at the first rung that covers it:

| What has to change | Mechanism | Code the user owns |
|---|---|---|
| Parameters, resources, input column, a variant run | Run/collect flags (`--set`, `-i`, `-l`, `-c/-T/--mem`, `--modal-gpu`, ...) | None |
| Which designs are submitted | A `-f/--filter` module with `apply_filter(df) -> df` | One small `.py` |
| How the binary is found / its environment on SLURM | The user's activation script (`SAPIA_ACTIVATE_<NAME>`) | Their activation `.sh` |
| The per-task shell commands, for one run | `-s/--script path/to/my.sh` | One `.sh` |
| One hook or field, permanently | `get_builtin(<name>).with_overrides(...)` in a user `spec.py` | `spec.py` + the replaced hook |
| Several pieces, or the pieces are tightly coupled | `sapia fork-tool <name> [dest]`, then edit the copy | The whole tool folder |
| The built-in itself is wrong or missing a general feature | Edit `src/prosapia/tools/<name>/` in place (contributing to prosapia) | Nothing extra; it ships |

If the change is general-purpose (a missing flag, a bug, a Modal image), prefer fixing the built-in in place over teaching every user to override it. If it is site- or project-specific, keep it out of `src/prosapia/tools/` and use a user-side mechanism.

## Rung 1 — flags, no code

Many "edits" are already a flag. Read the tool's `add_run_<name>_args` / `add_collect_<name>_args` and the base flags (`sapia run <name> -h`) before writing code. Common ones: `-i/--input-column` (feed a different column), `-l/--dir-label` (a same-tool variant with its own output dir and columns, e.g. a second seed), `--table-label` (label the child table of a create run), the resource flags `-g/-c/-T/--mem/-C` (override `#SBATCH` and Modal `RESOURCES`), `--modal-gpu`, and passthroughs like proteinmpnn's `--set '--ca_only'` that forward raw flags to the binary. `--force` re-submits/re-collects designs that are already `OK`.

## Rung 2 — a filter module

When the need is "run the tool on only some designs" (a threshold, a sample, a subset by name), write a filter instead of touching the tool. A filter is a module defining exactly `apply_filter(df: DataFrame) -> DataFrame`; the driver applies it to the source table at submit time, before `ctx.ready`'s input-column and resume checks. The signature is fixed, so parameters come from environment variables (set in `.env` or inline). See `docs/writing-a-filter-function.md` and `examples/filters/test_filter.py`.

```python
from pandas import DataFrame


def apply_filter(df: DataFrame) -> DataFrame:
    return df[df["boltz_complex_plddt"] > 0.8]
```

## Rung 3 — the activation script

How a tool's binary and environment are found on SLURM is not in the tool at all: each `.sh` calls `sapia_activate SAPIA_ACTIVATE_<NAME>`, which sources the user's activation script named in `.env`. Changing conda env, module, checkout path, interpreter pins (`RFDIFFUSION_PYTHON`, `PROTEIN_MPNN_PYTHON`, `USALIGN_BIN`, `PIPELINE_PYTHON`), input paths (`PROTEIN_MPNN`, `ROSETTA`, `AF3_*`) or per-tool runtime setup (framework caches on node-local scratch) all belong there — never hard-code them into the tool's `.sh`. `sapia init --config` scaffolds templates under `activation/`; `docs/configuration.md` lists each tool's variables and whether they go in the activation script or in `.env` (values read at submit time in Python must be in `.env`). Under `--executor modal` activation is skipped and the tool's `modal_image.py` provides the environment instead.

## Rung 4 — a custom task script for one run

To change the per-task shell commands without touching Python, copy the tool's `.sh`, edit it, and pass it with `-s/--script path/to/my.sh`. The manifest is unchanged, so the copy must keep reading the same `$SAPIA_LINE` fields and writing output where the collector expects it. The script's directory becomes `SAPIA_TOOL_DIR`, so sibling files the original calls through `"${SAPIA_TOOL_DIR}/..."` (e.g. a `_worker.py`) must sit next to the copy; likewise `--executor modal` looks for `modal_image.py` next to the script. The log folder is named after the script stem, so a custom script gets its own `<stem>_logs/`. If the custom script should become the default, move up to rung 5 and override `default_script`.

## Rung 5 — `with_overrides` in a user `spec.py`

`Tool` is a frozen dataclass; `get_builtin(name)` returns a bundled one and `with_overrides(**fields)` returns a copy with some fields swapped. In a user tool folder (on `$PROSAPIA_TOOLS_DIR`, default `./tools`), write a `spec.py` that reuses everything and replaces only what differs:

```python
from pathlib import Path

from prosapia.core import get_builtin

from .collect_my_rfdiffusion import collect_my_rfdiffusion

TOOL = get_builtin("rfdiffusion").with_overrides(
    collect_fn=collect_my_rfdiffusion,
    # default_script=str(Path(__file__).parent / "rfdiffusion.sh"),
)
```

Overridable fields (`ToolOverrides` in `core/tool.py`): `name`, `action`, `description`, `default_script`, `default_input_column`, `build_manifest_fn`, `collect_fn`, `add_run_args_fn`, `add_collect_args_fn`. `spec.py` is loaded as part of a synthetic package, so relative imports of sibling modules work. A replacement hook must honour the same contract as the original: a new `build_manifest_fn` must emit the fields the `.sh` cuts (or ship a matching `default_script`), a new `collect_fn` must read the layout the `.sh` writes, and a new `add_run_args_fn` replaces the original flags wholesale — if the built-in manifest builder reads `ctx.args.<flag>`, either keep that flag or replace the builder too. To extend rather than replace, import the original's pieces from the built-in module path (e.g. `from prosapia.tools.proteinmpnn.run_proteinmpnn import add_run_proteinmpnn_args`) and wrap them.

Reach for this when one or two pieces change and the rest should keep tracking upstream. Once you are replacing most of the hooks or copying the `.sh` and its siblings, fork instead.

## Rung 6 — `sapia fork-tool`

`sapia fork-tool <name> [dest] [--tools-dir DIR]` copies the whole built-in folder (minus `__pycache__`) into the first entry of `$PROSAPIA_TOOLS_DIR` (default `./tools/<folder>`), after which it is an ordinary user tool: edit its `spec.py`, `run_*`, `collect_*`, `.sh`, `modal_image.py` freely, following the authoring-a-tool conventions. `default_script` is built from `Path(__file__).parent`, so the fork automatically runs its own `.sh`, worker and Modal image. A fork no longer receives upstream fixes, so keep the diff to the original small and note what changed in the module docstrings.

## Naming: shadow or register alongside

Discovery (`tool_registry.discover`) loads the built-ins first, then every dir on `$PROSAPIA_TOOLS_DIR` in order, and keys tools by `Tool.name` — not the folder name. Later dirs win. So a user `spec.py` (rung 5 or 6) either:

- **keeps `name`** → it **shadows** the built-in: `sapia run <name>` now uses the user version everywhere that tools dir is on the path; or
- **changes `name`** (via `with_overrides(name=...)` or editing the fork's `spec.py`) → it registers a **new tool alongside** the built-in.

The name is also the **leaf** of the output dir and every column the tool writes (`run_dir/<table>/<name>[_<dir_label>]/`, `<name>_status`, `<name>_path`, `<name>_<key>`; see `core/naming.py`). Renaming therefore changes column names: downstream tools whose `default_input_column` is the original (`rfdiffusion_path`, `boltz_path`, ...) will need `-i <newname>_path`, and a run and its collect must use the same tool. Shadowing keeps the columns identical, so downstream steps are unaffected. For "the same tool with different settings" in one run_dir, prefer `--dir-label` over a rename.

## Rung 7 — editing the built-in in place

When working on prosapia itself and the change belongs to everyone, edit `src/prosapia/tools/<name>/` directly, following the authoring-a-tool conventions. Keep it backward compatible where users depend on it: new flags get defaults that preserve current behaviour, manifest field order and the `.sh` must change together, and output/column names should not change silently (they are the table contract for downstream tools and existing run_dirs). Update `docs/tools/<name>.md` and `docs/configuration.md` when flags, env variables or activation needs change, and remember that user forks and `with_overrides` wrappers import these modules and hook names.

## Checklist

1. Read the tool's folder, its docs, and `sapia run <name> -h` / `sapia collect <name> -h`; state what exactly must change.
2. Pick the lightest rung that covers it (flags → filter → activation → `-s` → `with_overrides` → fork → in place).
3. For rungs 5–6, decide shadow vs new name and check the column-name consequences for downstream `--input-column`s.
4. Keep manifest ↔ `.sh` ↔ collector consistent: the fields the builder emits are the fields the script cuts, and the files the script writes are the files the collector reads.
5. If the tool runs on Modal, make sure a `modal_image.py` sits next to whichever `.sh` will run, and that manifest paths still go through `volume_path`.
6. Verify: `uv run sapia run <name> -h` shows the expected tool and flags (a user tool must be on `$PROSAPIA_TOOLS_DIR`), then sanity-check the manifest/collect logic locally with `uv run` before a full run on SLURM or Modal.
