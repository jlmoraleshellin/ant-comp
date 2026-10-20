---
name: editing-a-tool
description: >-
  How to change an existing (bundled or custom) tool: picking the lightest mechanism that works (run flags, a `--filter` module, a custom `--script`, the activation script, `with_overrides`, `sapia fork-tool`, or editing in place) and the naming/discovery rules that decide whether a change shadows the original or registers alongside it. Use when asked to edit, customize, override, fork or shadow a tool, or to change how it is activated. For a brand-new tool use authoring-a-tool.
---

# Editing a pipeline tool

Pick the lightest mechanism that does the job: every rung down is more code the user owns and must keep in sync with upstream. The anatomy of a tool and the hook contracts are in `authoring-a-tool`; this skill is about changing one that exists. Read the real code first, it is the source of truth: the tool's folder (`spec.py`, run/collect modules, the `.sh`, `modal_image.py`), `core/tool.py` (`Tool`, `ToolOverrides`, `with_overrides`), `core/tool_registry.py` (`discover`, `get_builtin`), `cli/fork_tool.py`, `core/naming.py`. The bundled tools are in the installed package; the prose docs (`docs/writing-a-tool.md`, `docs/configuration.md`, `docs/writing-a-filter-function.md`) exist only in a prosapia checkout.

## The ladder

| What has to change | Mechanism | Code the user owns |
|---|---|---|
| Parameters, resources, input column, a variant run | Run/collect flags (`--set`, `-i`, `-l`, `-c/-T/--mem`, ...) | None |
| Which designs are submitted | `-f` filter module with `apply_filter(df) -> df`, or `-w column=value` | One small `.py` |
| How the binary is found / its environment on SLURM | The activation script (`SAPIA_ACTIVATE_<NAME>`) | Their activation `.sh` |
| The per-task shell commands, for one run | `-s/--script my.sh` | One `.sh` |
| One hook or field, permanently | `get_builtin(<name>).with_overrides(...)` in a user `spec.py` | `spec.py` + the replaced hook |
| Several pieces, or tightly coupled pieces | `sapia fork-tool <name> [dest]`, then edit the copy | The whole tool folder |
| The built-in itself is wrong | Edit the bundled tool in the prosapia checkout | Nothing extra; it ships |

- **Flags.** Read `sapia run <name> -h` and the tool's `add_run_*_args` / `add_collect_*_args` before writing code: `-i` feeds another column, `-l` gives a same-tool variant its own leaf and columns, `--table-label` labels a `create` child, `--force` redoes `OK` rows, and passthroughs like proteinmpnn's `--set` forward raw flags.
- **Filter.** The signature is fixed, so parameters come from environment variables; the driver applies it to the source table before the input-column and resume checks.
- **Activation script.** Conda env, module, checkout path, interpreter pins and per-tool runtime setup belong in the activation script named in `.env`, never hard-coded in the tool's `.sh`; values read at submit time in Python must be in `.env`. Under `--executor modal` activation is skipped and `modal_image.py` provides the environment.
- **`-s`.** The copy must keep reading the same `$SAPIA_LINE` fields and writing where the collector expects; its directory becomes `SAPIA_TOOL_DIR`, so siblings the original calls (a `_worker.py`) and `modal_image.py` must sit next to it; its logs go to `<stem>_logs/`.
- **`with_overrides`.** Overridable: `name`, `action`, `description`, `default_script`, `default_input_column`, `build_manifest_fn`, `collect_fn`, `add_run_args_fn`, `add_collect_args_fn`. A replacement hook honours the original's contract: a new builder emits the fields the `.sh` cuts, a new collector reads the layout the `.sh` writes, and a new `add_run_args_fn` replaces the flags wholesale, so keep any flag the builder reads. Once most hooks are replaced, fork instead.
- **Fork.** `fork-tool` copies the folder into the first entry of `$PROSAPIA_TOOLS_DIR` (default `./tools`); `default_script` is built from `Path(__file__).parent`, so the copy runs its own `.sh`, worker and image. A fork stops receiving upstream fixes: keep the diff small and note changes in the module docstrings.

## Shadow or register alongside

Discovery loads the built-ins, then each dir on `$PROSAPIA_TOOLS_DIR` in order, keyed by `Tool.name` (not the folder name); later dirs win.

- **Keep `name`: the fork shadows the built-in.** `sapia run <name>` uses it everywhere that tools dir is on the path, and every column name stays identical, so downstream steps are unaffected. This is the default way to change a bundled tool in this workspace: the install is not editable, so an edit to the bundled source takes effect only after prosapia's `dev` is pushed, the pin bumped and the cluster venv re-bootstrapped, while a fork under `tools/` is live on the next run.
- **Change `name`: it registers a new tool alongside.** The name is also the leaf of the output dir and of every column (`<name>_status`, `<name>_path`), so renaming changes column names, downstream `default_input_column`s need `-i <newname>_path`, and a run and its collect must use the same tool. For the same tool with other settings in one run_dir use `-l`, not a rename.

## Hard-won rules

- **A hard-coded column list lives in several places.** Adding or renaming a column means updating the worker, the collector, and the error-row column list together; a column known to one list and not the others is silently dropped or gets the wrong dtype, and the error row must carry it as NA.
- **Manifest, `.sh` and collector move together**: the fields the builder emits are the fields the script cuts, and the files the script writes are the files the collector reads. On Modal a `modal_image.py` must sit next to whichever `.sh` runs, and manifest paths must go through `volume_path`.
- **Column names are a contract.** A bundled tool's names are what downstream tools and existing run_dirs expect; new flags need defaults that preserve current behaviour. Say which existing columns change meaning, and give the edited tool a new `-l` label (see `verifying-a-tool`).
- **The edit is not real until it is committed.** Tell the user to commit it the moment it is verified, say in the hand-back whether it is committed, and edit the tool's skill only together with the commit; a working-tree edit that a skill documents is a trap the next session cannot see.

## Checklist

1. Read the tool's folder and `sapia run <name> -h` / `sapia collect <name> -h`; state what exactly must change.
2. Pick the lightest rung (flags, filter, activation, `-s`, `with_overrides`, fork, in place).
3. For overrides and forks, decide shadow vs new name and check the column-name consequences for downstream `-i`.
4. Keep manifest, `.sh` and collector consistent; update every hard-coded column list.
5. Check `sapia run <name> -h` shows the expected tool and flags (a user tool must be on `$PROSAPIA_TOOLS_DIR`), sanity-check manifest and collect locally, then verify per `verifying-a-tool`.
6. Make the edit durable: commit, then update the skill.
