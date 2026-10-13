# ant-comp — a protein-design workspace driven by Claude Code

A workspace where you describe a protein-design problem in plain language and Claude plans it, runs it, and reads the results back. Nothing is installed locally: every tool runs on **Modal** (a container per task) or on the **VIB DataCore** (SLURM jobs). Built on [`prosapia`](https://github.com/jlmoraleshellin/prosapia) (CLI: `sapia`).

---

## Install

**You need:** Python 3.13, [`uv`](https://docs.astral.sh/uv/), [Claude Code](https://claude.com/claude-code), and an account on at least one backend (Modal, VIB DataCore, or both).

**1. Clone and install**

```bash
git clone git@github.com:jlmoraleshellin/ant-comp.git
cd ant-comp
uv sync
```

`uv sync` pulls `prosapia` from GitHub at the commit pinned in `uv.lock`; it needs read access to that repo.

**2. Set up Modal** — skip if you only use vib

```bash
uv run modal setup
uv run sapia modal-shell --cmd 'sapia run --help'   # should list the tools
```

Set your personal volume values in `.env`:
```bash
SAPIA_MODAL_RUNS_VOLUME="sapia-runs-<your-name>"
```

If you don't do this, your results will appear in mine.

**3. Set up the VIB DataCore** — skip if you only use Modal

Add the cluster to `~/.ssh/config` as `vib`, and check that `ssh vib true` returns silently:

```
Host vib
    HostName <login-node-hostname>
    User <your-username>
```

Set your personal values in section 1 of the vib part of `.env`:

```bash
SAPIA_VIB_HOST="vib"
SAPIA_VIB_WORKSPACE="/data/groups/csb/<...>/<you>/ant-comp"   # yours; anywhere on group storage
SAPIA_VIB_ACCOUNT="anastassia_vorobieva"                      # required on every job
```

Then build the workspace — creates it, syncs this repo into it, and builds its venv from the `Miniconda3` module on a compute node:

```bash
bash scripts/vib.sh bootstrap    # ~5 min
```

For vib alone you need only `bash`, `ssh`, `rsync`, `timeout` and `python3` on your machine — `uv` is needed for Modal, never for vib. Everything else about the vib connection is in [`docs/vib.md`](docs/vib.md).

Auth is a short-lived SSH certificate via a Microsoft browser sign-in. When it lapses the next `ssh` opens a sign-in page and waits — so if Claude reports that ssh hung, that is your cue to sign in.

**4. Start working**

```bash
claude --agent thinker
```

It will ask which backend to use, and about the design decisions it foresees, before running anything.

---

## How it works

`prosapia` is a workbench, not a pipeline: **a design is a row, a generation of designs is a table**, and each tool adds columns. A `create` tool (rfdiffusion3, proteinmpnn) mints a child table; an `update` tool (boltz, usalign, cms) annotates the table in place. Columns are prefixed per tool, and `<leaf>_status == "OK"` is the only success signal — a blank means "not applicable", never zero.

The corollary worth internalising: **if a measurement produces one value per design, it belongs in a table column, not an ad-hoc script.** Columns can be filtered on, carry a status, and survive into child tables. A script's output is a file the next session cannot see.

Claude Code subagents cannot spawn subagents, so the setup is two layers:

```
you ── the science ──▶ thinker            (main session, Opus)
                          │ one step at a time
              ┌───────────┴───────────┐
              ▼                       ▼
         modal-worker             vib-worker        (Sonnet)
              ▼                       ▼
       Modal containers        SLURM array jobs
```

The `thinker` owns the goal and reads the tables but never runs `sapia`. The workers execute, report, and **stop** — neither chains into the next tool on its own. `tool-creator` builds a new tool when no existing one answers a measurement. Per-tool skills in `.claude/skills/` hold the flags and traps.

### The two backends

Same CLI, same tables, same skills. A `run_dir` lives on **one** backend — nothing syncs them.

| | Modal | VIB DataCore |
| --- | --- | --- |
| Environment | `modal_image.py` per tool | activation script per tool |
| Data | the `sapia-runs` Volume | `$SAPIA_VIB_WORKSPACE/outputs` |
| Completion signal | `.exit` files | `squeue` / `sacct` |
| Waiting | no queue; first run builds an image (10–15 min) | queue waits, minutes to hours |
| Cost | billed per second | free, but nodes are shared with the lab |
| Custom tools | **all six** | `mkcomplex`, `chainsel`, `ringfit` only |

On vib the workspace is an rsync'd copy of this repo with its own venv, resynced before every submit — so **a local edit to a tool or activation script is live on the next run, with no commit and no push**. That is the vib equivalent of Modal shipping your working tree into the image.

---

## What is in the repo

```
CLAUDE.md          standing context, read at every session start
.claude/agents/    thinker, modal-worker, vib-worker, tool-creator
.claude/skills/    one SKILL.md per tool, plus binder-campaign and the process skills
tools/             custom prosapia tools
activation/vib/    the SLURM half of each custom tool (what modal_image.py is on Modal)
scripts/           vib.sh — the one way to reach the cluster
docs/              vib.md
campaigns/         campaign reports — the scientific record
```

**Custom tools:** `cms` (interface contact surface + shape complementarity), `chainsel` (extract or merge chains), `mkcomplex` (rebuild a complex around a designed binder), `ringfit` (does a binder straddle two protomers of an oligomer), `bindcraft2` (a whole binder campaign in one step), `atomium` (ProteinMPNN-like design with a private model).

Ask the thinker to write a campaign report into `campaigns/` before closing a session. Skills accumulate tool knowledge; nothing accumulates scientific knowledge unless you write it down. See `campaigns/20260928_7ojg_slyb_binder.md` for the shape.

---

## Day to day

Pull before you start. If you have uncommitted edits, park them: `git stash`, `git pull`, `git stash pop` (`pop` applies and deletes; `git stash apply` keeps a copy until you `git stash drop`). A conflicting `pop` leaves the stash in place — resolve, `git add`, then drop it manually.

Share improvements back on a branch. Skills and campaign reports are the point of the shared repo: a trap you wrote down is a trap nobody else hits. Keep run outputs, loose `.cif`/`.pdb` dumps and `.venv` out of commits.

To move the prosapia pin: `uv lock --upgrade-package prosapia && uv sync`, test both backends, then commit `uv.lock`. **This does not update vib by itself** — re-run `bash scripts/vib.sh bootstrap`, or the cluster keeps running the old library silently (`bash scripts/vib.sh status` shows the mismatch).

---

## Things that will bite you

**Both backends**

- **Never run `sapia` on your own machine** — neither `/runs` nor the cluster filesystem exists locally.
- **Never `cd` into a run_dir** — stored paths become relative to it and break for every later tool.
- **`No designs to submit.` exits 0.** A run that queues nothing looks like success; it is almost always a wrong input column (ProteinMPNN defaults to `rfdiffusion_path`, so from an rfd3 table you must pass `-i rfdiffusion3_path`).
- **Exit 0 does not mean it worked.** Several tools write errors as data and still exit 0. `<leaf>_status` is the proof.
- **High confidence is not proof.** A confident prediction does not mean the sequence folds to its designed backbone — that is what `usalign` self-consistency is for. For binders, fold and pose are separate questions.
- **The broken predictions produce the best-looking numbers.** A huge interface is evidence of a collapsed target until the target RMSD says otherwise.

**Modal only**

- A tool's first run builds its image inside the submit call — allow 15+ minutes before assuming it is stuck.
- Never delete a weight Volume to tidy up; `rfd3-checkpoints` has to be populated by hand.

**VIB DataCore only**

- **Every job needs `-A`.** `sbatch`, `srun` and `salloc` are all rejected without it despite your default association, and prosapia surfaces this only as `sbatch exited 1`.
- **Use `--partitions <name>:<gpu_count>`, never the bare name** — prosapia's GPU-count query cannot parse this cluster's `gpu:h100:4(S:0-1)` and dies with `ValueError: invalid literal for int()`.
- **GPU jobs need `--partitions` at all**: the default partition has none. Where a GPU type has a `_co_pi` partition, use only that one (the plain ones belong to other groups), and stay under half its GPUs (`--max-gpu-fraction` defaults to 0.5).
- **Do real work in an allocation, not on the login node** — it is throttled, and a venv build there wedges partway through.
- **`-g 0` for CPU-only tools**, or they request a GPU on a partition that has none and pend forever.
- **Never add `--delete-excluded` to the sync** — `outputs/` and `.venv/` live inside the workspace and that flag deletes them outright.

**Network (imec):** TLS interception breaks `modal volume get`/`put` while `modal volume ls` works. Data lives on the Volume anyway, so normal work is unaffected. Do not weaken TLS verification to work around it.

---

## Extending

Ask the thinker for a measurement in scientific terms; if no tool covers it, `tool-creator` scaffolds one (`spec.py`, `run_*.py`, `collect_*.py`, the `.sh`, and `modal_image.py`). That gives a Modal-only tool — running it on vib additionally needs an activation script in `activation/vib/` and a `SAPIA_ACTIVATE_<NAME>` entry in `.env`.

Teach the workspace what you learn by editing the relevant `SKILL.md`. Measured facts with numbers beat general advice.

For your own agent variant, copy rather than edit: `cp .claude/agents/thinker.md .claude/agents/thinker-membrane.md`, change the `name:` in the frontmatter to match, and run `claude --agent thinker-membrane`. Editing the shared file gives every colleague a merge conflict and imposes your preferences on their campaigns. Corrections and newly-found traps are the exception — those belong in the shared file.
