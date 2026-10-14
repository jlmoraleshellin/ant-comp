# What this project is

A **protein-design workspace** built on [`prosapia`](https://github.com/jlmoraleshellin/prosapia)
(CLI: `sapia`), running on one of **two backends**: **Modal** (a container per task) or the
**VIB DataCore** (SLURM array jobs over ssh). No design software is installed locally —
every tool runs remotely, and the data stays there too, on a Modal Volume or on cluster
group storage.

The intended way of working: **Claude drives the campaign, the backend does the compute.**
Claude never holds the data; it submits steps, waits for them, collects them, and reads the
tables.

The science is identical either way — same CLI, same tables, same lineage, same skills.
What differs is the execution machinery. **A single `run_dir` lives on one backend** — the
Volume and the cluster filesystem are separate worlds, and nothing syncs them.

**The vib workspace.** `$SAPIA_VIB_WORKSPACE` is an rsync'd copy of this repo on cluster
storage with its own `.venv`, and the worker resyncs it before every submit. That is
the vib equivalent of Modal shipping your working tree into the image: edit a tool or an
activation script locally and it is live on the next run, **with no commit and no push**
(verified end to end). All 17 tools register there.

Why the venv must live inside the workspace: prosapia calls `load_dotenv()` at import
(`core/base_run.py`) and `find_dotenv` walks up from `base_run.py`'s own path in
`site-packages`, **not** from the cwd. So the venv's location decides which `.env` supplies
submit-time config, while the task prelude reads the `.env` in the submit cwd. With the venv
inside the workspace both resolve to the same file — this repo's `.env`, which is therefore
the single place tool activation is configured for vib.

**A tool needs two things on vib**: its code in the synced `tools/` (automatic), and a
`SAPIA_ACTIVATE_<NAME>` entry. Registered is not the same as runnable — without the second,
a tool appears in `sapia run --help` and then every task dies at activation. Current state:
`mkcomplex`/`chainsel`/`ringfit` runnable (they need nothing beyond prosapia's own deps);
`cms`/`atomium`/`bindcraft2` not ported (they need real software installed on the cluster);
`pyrosetta`/`openfold3` registered with no activation script anywhere.

## prosapia in one page

A **workbench, not a pipeline.** There is no fixed order of steps — there is a shared
tabular data format and tools that consume and produce it. You compose a workflow
dynamically, forking and back-tracking as the science demands.

Two principles:

1. **A design is a row; a generation of designs is a table.** Each row is keyed by a design
   `name`; each tool contributes columns (a structure path, a sequence, a pLDDT, an RMSD).
2. **When a protein diverges in sequence or structure it is a child, not the same protein**
   — so it needs a new table.

Everything happens inside a **`run_dir`**, minted by `sapia new_run` and never by a tool:

```
run_dir/
├── _registry.tsv      catalog of every table + its lineage (parent, gen)
├── table0.tsv         one row per design, keyed by `name`
├── table1.tsv         a child table (another generation)
├── .manifests/        transient per-run manifests
└── table0/            per-table, per-tool outputs on disk
    └── rfdiffusion3/
```

**Two phases per tool.** `sapia run` fans tasks out (each writes only its own files);
`sapia collect` folds those files into the table. This keeps table writes safe under heavy
parallelism and repeatable on reruns.

**Two kinds of tool.** The `action` decides where output lands, and **you never name the
output table** — the driver derives it:

- **`create`** (rfdiffusion3, proteinmpnn) mints a **child table**, `gen+1`, and links each
  new row to its parent. Produces new entities.
- **`update`** (boltz, alphafold3, usalign) annotates the **same table in place**, adding
  columns. Derives a property of designs that already exist.

Columns are leaf-prefixed per tool (`boltz_ptm`, `proteinmpnn_score`), and
`<leaf>_status == "OK"` marks a design that actually succeeded. Reruns skip already-`OK`
rows unless `--force`.

## How work gets done here: the agents

Claude Code subagents **cannot spawn further subagents**, so the setup is two layers, not
three:

```
main session = thinker                    (claude --agent thinker; Opus)
   ├─ modal-worker subagent               (Sonnet; Bash/Read/Skill)   → Modal
   └─ vib-worker subagent                 (Sonnet; Bash/Read/Skill)   → VIB DataCore
         └─ reads .claude/skills/<tool>/SKILL.md on demand
```

- **`thinker`** (`.claude/agents/thinker.md`) — owns the scientific problem: goals, what to
  try next, reading result tables, what to keep. **Never runs `sapia`, `modal` or `ssh`
  itself.** Delegates intent ("20 backbones, length 90–110") and requires the run_dir,
  table, row count and failures back. **Asks the user which backend** at the start of a
  session if they haven't said, then uses that one worker throughout.
- **`modal-worker`** (`.claude/agents/modal-worker.md`) — executes on Modal: the
  workstation, tool images, the `.exit` wait loop.
- **`vib-worker`** (`.claude/agents/vib-worker.md`) — executes on the VIB DataCore: ssh to
  the login node, `sbatch`, partitions, the `squeue`/`sacct` wait loop.
- Both workers report back and **stop**; neither chains into the next tool on its own.
- **Per-tool skills** (`.claude/skills/<tool>/SKILL.md`) — flags,
  verified invocations, collected columns and the specific traps of each tool. Loaded by
  the worker instead of re-reading the full docs. **They were written against
  Modal**: the tool flags, input/output columns and scientific traps hold everywhere, but
  anything about Volumes, images, `modal-shell`, `--modal-gpu` or `.exit` files does not
  apply on vib.

Start a campaign with `claude --agent thinker`.

## Execution model A: Modal

`run_dir`s live **only on a Modal Volume**, mounted at `/runs`. `sapia` runs next to it in
a small container, the **workstation**, not on this machine:

```bash
sapia modal-shell                       # interactive shell; cwd is /runs
sapia modal-shell --cmd '<command>'     # one command — ALWAYS exits 0, see below
```

Non-negotiables, each learned the hard way:

- **`sapia modal-shell --cmd` ALWAYS returns 0. Never test `$?` after it.** Measured
  2026-10-02: `--cmd 'exit 3'` → `rc=0`; `--cmd 'ls /nonexistent'` → `rc=0` with the real
  `No such file or directory` on stderr, i.e. the command ran, failed, and the code was lost
  on the way out. prosapia is not at fault — `modal_shell_from_args` is literally
  `raise SystemExit(subprocess.run(...).returncode)`; the loss is in `modal shell --cmd`
  itself or in the base64 `bash <(...)` wrapper prosapia uses to survive quoting. **This
  silently hides a failed submit**: a tool that raises while building its manifest (a bad
  column, a refused input) queues nothing, and the caller sees success. Detect it properly:
  * **authoritative, on disk:** `<out_dir>/<script>_logs/<script>_modal.json` exists and
    holds `{"app_id", "n_tasks"}`. No file ⇒ nothing was queued.
  * **secondary, in the captured stdout:** the line `Submitting N designs`. It is printed
    only, never written to disk, so it is checkable only in the call's own output.
  Treat `No designs to submit.` and a Python traceback as failures, however the shell exits.
- **Never run `sapia` outside the workstation.** Locally there is no `/runs`, so the run
  fails or, worse, writes paths no task container can resolve.
- **Run `sapia` from `/runs`** (the workstation's cwd) and use the relative `run_dir` that
  `new_run` printed. Never `cd` into a run_dir — stored paths become relative to it and
  break for every later tool.
- **`sapia run` is detached.** It returns once tasks are queued. The `.exit` files below are
  the only reliable completion signal.
- **A tool's first run builds its Modal image inside the submit call** — Boltz took over
  10 minutes. Allow 15+ minutes before assuming a submit is stuck. Later runs are seconds.
- **`-g 0` for CPU-only tools.** `--gpus-per-task` defaults to 1.

### Task status: the `.exit` loop

`<run_dir>/<table>/<leaf>/<script>_logs/` holds, per task, `<script>_<id>.out`,
`.err`, `.exit` (the exit code — written even on failure, `255` if the wrapper itself
failed), plus `<script>_modal.json` = `{"app_id", "n_tasks"}`.

| State | How to tell |
| --- | --- |
| done | `n_tasks` `.exit` files, all `0` |
| failed | an `.exit` that isn't `0` → read the matching `.err` |
| running | `.exit` files missing, app still running |
| killed | `.exit` files missing, app `stopped` (timeout/OOM) → `modal app logs <app_id>` |

The loop is **submit → poll `.exit` → check codes → collect**. Poll every 30–60 s (each
workstation call costs 5–10 s of cold start).

`modal` CLI commands (`app list`, `app logs`, `volume ls`) run **locally**, not in the
workstation. Set `NO_COLOR=1` before parsing their output — ANSI codes break JSON parsing.

## Execution model B: the VIB DataCore (SLURM)

No containers here. `sapia` runs from the workspace and submits SLURM array jobs; each
tool's environment comes from an activation script. Everything goes through **one script,
`scripts/vib.sh`** — the only sanctioned way to reach the cluster; full contract in
`docs/vib.md`. Its personal inputs are the three `SAPIA_VIB_*` values in section 1 of
**this repo's `.env`** (host, workspace, account); the `SAPIA_ACTIVATE_*` entries below them
are read by prosapia itself on the compute node.

```bash
bash scripts/vib.sh exec --sync --on-node '<command>'   # push the working tree, then run in a small srun (use for `sapia run`)
bash scripts/vib.sh exec --on-node '<command>'          # `sapia collect`, anything reading tables (2 cpus, 8G, 30 min)
bash scripts/vib.sh exec '<command>'                    # login node: squeue / sacct / sinfo polls and logs only
bash scripts/vib.sh status                    # marker, venv, prosapia commit vs pin, uv

# anything heavier than a scheduler query goes in an allocation:
#   bash scripts/vib.sh exec 'srun -A "$SAPIA_VIB_ACCOUNT" -c 8 -t 1:0:0 bash -lc "<command>"'
```

Non-negotiables, each learned the hard way:

- **Sync before every submit** (`exec --sync`). Only `tools/ activation/ utils/ .env`
  travel, and `--delete` acts only inside those directories, so `outputs/` and `.venv/`
  cannot be reached. The sync refuses unless the remote workspace holds a
  `.sapia-workspace` marker. Not needed before poll calls.
- **Do real work in an allocation, never on the login node.** `srun -A "$SAPIA_VIB_ACCOUNT"
  -c 8 -t 1:0:0 bash -lc '<cmd>'` (or `salloc -c 8 -t 12:0:0` for a human shell). `-A` is
  required on `srun`/`salloc` just as on `sbatch`, and `bash -lc` is what makes `module`
  available. Measured: a venv build on the login node wedged twice, and took minutes on a
  compute node. Scheduler queries (`squeue`/`sacct`/`sinfo`) are fine on the login node.
- **No `uv` on vib.** The venv is built by `bash scripts/vib.sh bootstrap` from the
  `Miniconda3` module's Python (3.12) with stdlib `venv`/`pip`, pinned to the prosapia
  commit read from `uv.lock`. On this cluster `uv` exists only under one user's home, so it
  cannot be assumed. rsync never touches the venv, so after a pin bump the cluster silently
  runs the old library until someone re-runs `bootstrap` — `vib.sh status` shows the
  mismatch. If `uv.lock` points at the editable `../prosapia`, `bootstrap` stops and needs
  `--commit <pushed sha>`.
- **Commands always start in `$SAPIA_VIB_WORKSPACE`, and never `cd`.** `vib.sh exec`
  does the `cd` itself and refuses `cd`, `..`, `~`, absolute paths outside the workspace and
  an `rm` outside the workspace (an accident guard, not a sandbox; removals inside it are
  allowed). Why the cwd matters: a SLURM task
  starts in the submit directory and its prelude sources `.env` relative to it, and
  `PROSAPIA_TOOLS_DIR` defaults to the literal `tools` resolved against the cwd. Submit from
  anywhere else and the custom tools vanish and every task dies with
  `sapia: set SAPIA_ACTIVATE_<TOOL> in your .env`.
- **`--partitions <name>:<gpu_count>`, never the bare name.** prosapia's GPU-count query
  splits `sinfo -o %G` on `:` and this cluster reports `gpu:h100:4(S:0-1)`, so the bare
  form dies with `ValueError: invalid literal for int() with base 10: '0-1)'`.
- **Run_dirs are absolute**, minted with `--base "$SAPIA_VIB_WORKSPACE/outputs"`. Group
  storage is shared by login and compute nodes, so an absolute path resolves wherever a
  task lands. Still never `cd` into a run_dir.
- **GPU jobs need `--partitions`.** The default partition `gp_64C_128T_512GB` has no GPUs,
  so a GPU task submitted without one never runs. Check `sinfo -o '%P %G'` and `sinfo -s`;
  don't pick a partition unilaterally for a large batch. `--modal-gpu` is Modal-only.
- **Stay under half a partition's GPUs.** Each holds 4–16 in total across its nodes;
  `--max-gpu-fraction` defaults to `0.5` and caps concurrency accordingly. Leave it there
  unless the caller has a reason — the capacity is the lab's, not ours.
- **Where a GPU type has a `_co_pi` partition, only the `_co_pi` one may be used** — it is
  this group's entitlement, and the plain `gpu_h100_*` / `gpu_b300_*` belong to other groups
  even when idle. a100, l40s, `gpu_ds` and `gpu_short` have no `_co_pi` form and are fine
  as they are. Any GPU type, b300 included, is allowed through its `_co_pi` partition.
- **Every `sapia run` needs `-a $SAPIA_VIB_ACCOUNT`.** A default association is not enough:
  the job_submit plugin rejects the job with `! Missing slurm account`, and prosapia
  reports only `RuntimeError: sbatch exited 1` — the real reason is on stderr above the
  traceback. If the variable is empty, ask the user — never guess a value.
- **Certificates expire.** Auth is a short-lived SSH cert from a Smallstep CA via an Azure
  AD browser sign-in. If an ssh call hangs or says `Permission denied`, **stop and ask the
  user to sign in** — never retry in a loop. `vib.sh` wraps every ssh in `timeout`.
- **The login node and filesystem are shared with the lab.** Only scheduler queries and
  quick inspection belong there; everything else goes in an allocation. Never touch
  another user's directories, or anything under the shared `$SAPIA_VIB_SHARED_ACTIVATION`.
- **`-g 0` for CPU-only tools**, same as Modal.

### Task status: the SLURM loop

**There are no `.exit` files.** The scheduler is the source of truth; the log dir
(`<script>_<jobid>_<taskidx>.out` / `.err`) is the evidence. Record every job ID —
`--partitions` or >1000 tasks produce several arrays.

```bash
squeue -h -j <jobid> -o '%i %T %R'
sacct -n -P -X -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

| State | Meaning |
| --- | --- |
| `PENDING` | queued; `%R` gives the reason. `ReqNodeNotAvail` / `PartitionConfig` will never start — stop and report |
| `RUNNING` | keep polling |
| `COMPLETED` `0:0` | done — collect only when **every** array task is here |
| `FAILED` | read that task's `.err` |
| `TIMEOUT`, `OUT_OF_MEMORY`, `CANCELLED` | killed; report `Elapsed`/`MaxRSS` |

The loop is **submit → poll `squeue`/`sacct` → check states → collect**. Poll every
60–120 s; queue waits run from minutes to hours, and each call is a fresh ssh connection.

A task that fails with **empty `.out` and `.err`** died before its first statement — look
at the activation script and the prelude, not the tool.

## Local configuration

`prosapia` is installed **from GitHub, branch `dev`**, pinned to a commit in `uv.lock`:

```toml
[tool.uv.sources]
prosapia = { git = "https://github.com/jlmoraleshellin/prosapia", branch = "dev" }
```

Consequences, both of which matter:

- **The install is not editable.** Edits in a local `../prosapia` checkout have **no effect**
  on what runs — neither locally nor in the workstation image, which mounts the installed
  package from `.venv/lib/python3.*/site-packages/prosapia/`. To pick up new commits on
  `dev`: `uv lock --upgrade-package prosapia && uv sync`. To develop the library, switch the
  source back to `{ path = "../prosapia", editable = true }` for the duration.
- **`docs/` is not shipped.** The wheel carries `src/prosapia/` only. The library *source* is
  therefore readable at `.venv/lib/python3.*/site-packages/prosapia/`, but the prose docs
  that several skills cite (`docs/running-on-modal.md`, `docs/lineage-and-tables.md`, …) exist
  only in a checkout of the repo. Keeping a sibling `../prosapia` clone on `dev` is still
  worthwhile for that reason alone — it is just no longer the dependency.

This pin governs **Modal only**. The cluster env is installed and updated separately by
whoever maintains it; a `uv.lock` bump here does not move it, and the two can drift.

`dev` keeps moving. Read the source before trusting a detail, and prefer
`sapia run <tool> --help` (in the workstation, or on vib) over memory.

`.env` (no secrets; Modal auth lives in `~/.modal.toml`, vib auth in a short-lived SSH
cert). Note `SAPIA_EXECUTOR` is **not** set here: the Modal workstation exports
`SAPIA_EXECUTOR=modal` itself, and the cluster's own `.env` sets `slurm`.

| Variable | Value here | Meaning |
| --- | --- | --- |
| `SAPIA_MODAL_RUNS_VOLUME` | `sapia-runs` | The runs Volume. **Required.** Created on first use. |
| `SAPIA_MODAL_VOLUME_RFD3_CKPT` | `rfd3-checkpoints` | rfd3 checkpoints, mounted at `/checkpoints`. |
| `SAPIA_MODAL_VOLUME_BOLTZ_CACHE` | `boltz-cache` | Boltz weights + CCD, mounted at `/boltz_cache`. |
| `SAPIA_MODAL_VOLUME_BINDCRAFT_CACHE` | `bindcraft-cache` | BindCraft2's AlphaFold params (~5.3 GB). |
| `SAPIA_VIB_HOST` | `vib` | ssh alias for the DataCore login node (`~/.ssh/config`). |
| `SAPIA_VIB_WORKSPACE` | see `.env` | The rsync'd copy of this repo, holding its own `.venv` and `outputs/`. Every remote command runs from here. Per user. |
| `SAPIA_VIB_ACCOUNT` | `anastassia_vorobieva` | SLURM account for `-a`. **Required** — the cluster rejects jobs without an explicit `-A`. |
| `SAPIA_VIB_SHARED_ACTIVATION` | see `.env` | The lab's activation scripts for the built-in tools, read-only. The software they activate is installed on the cluster, not in this repo. |
| `SAPIA_ACTIVATE_*` | absolute paths | Read by prosapia on the compute node. Inert under Modal. |
| `CONDA_PREFIX` | see `.env` | Needed by the shared activation scripts; copied from the shared env's `.env`. |
| `RFD3_CKPT` | not in `.env` | Exported by the rfd3 activation script; `rfdiffusion3.sh` passes it as `ckpt_path` at task time. Never in `.env`: Modal tasks receive it too. |

The weight Volumes were deliberately renamed **without** a `sapia-` prefix so teammates who
don't use prosapia can share them. Note that prosapia's own defaults are still
`sapia-rfd3-checkpoints` / `sapia-boltz-cache`, and `get_named_volume` uses
`create_if_missing=True` — so a teammate missing these `.env` lines silently gets a new
**empty** Volume rather than an error. Keep these lines when copying the project.

**Never delete a weight Volume to tidy up.** `rfd3-checkpoints` is 2.5 GiB and must be
populated by hand (see below); `boltz-cache` re-downloads on first use.

## Validated pipeline

The full chain has been run end to end **on Modal**. Reference run: `outputs/20260927_211428_rfd3_denovo`
on `sapia-runs-test` (the previous test Volume; the current one is `sapia-runs`). The
timings below are Modal's; on vib, add the queue wait and drop the image-build time.
**The chain has not yet been validated end to end on vib** — treat a first run there as a
test, and start small.

| Step | Command | Result |
| --- | --- | --- |
| rfd3 de novo | `sapia run rfdiffusion3 <run_dir> --length 80-120 --num-designs 5` | 5 backbones → `table0`, ~3 min |
| ProteinMPNN | `sapia run proteinmpnn <run_dir> -t table0 -i rfdiffusion3_path --num-seq-per-target 2` | 10 sequences → `table1`, ~1.5 min |
| Boltz | `sapia run boltz <run_dir> -t table1` | 10 predictions annotated onto `table1`, ~16 min (incl. first weight download) |
| PyRosetta | `sapia run pyrosetta <run_dir> -t table1` | FastRelax (1 cycle) + ref2015 metrics onto `table1`, ~1.5 min for 2 designs (`outputs/20260928_113443_pyrosetta_test` on `sapia-runs`) |

All 10 Boltz predictions came back confident (confidence 0.91–0.97, pLDDT 0.93–0.97). Note
that high confidence is **not** proof the sequence folds to its designed backbone — the real
test is self-consistency: compare each prediction back to its parent backbone with
`usalign` (an `update` tool). That step has not been run yet.

## Known gaps

- **Tool weights are installed by hand.** There is no `sapia` command for it. rfd3's
  checkpoint was populated with a one-off `foundry install rfd3 --checkpoint-dir /checkpoints`
  from a container of the tool's image. A `setup()` hook in each `modal_image.py` plus a
  `sapia modal-setup <tool>` verb would let an agent do this itself.
- **`--length min-max` does not vary length within a batch.** rfd3 draws one length per
  batch, so `--num-designs 5` gave 5 backbones of the *same* length (103 aa). For a spread,
  use several batches (`--set n_batches=5`) or several design keys.
- **ProteinMPNN's `default_input_column` is `rfdiffusion_path`** — RFdiffusion, not
  rfdiffusion3. Coming from an rfd3 table you must pass `-i rfdiffusion3_path` or the run
  silently submits nothing.
- **Untested:** polling from inside a single long-running workstation container (unclear
  whether its Volume mount refreshes to show task commits). The worker therefore uses
  repeated short `modal-shell` calls, which is what was actually verified.
- **The custom tools are Modal-only** (see the top of this file). Giving one to vib means
  writing an activation script there and a `SAPIA_ACTIVATE_<NAME>` entry in the cluster's
  `.env` — which lives in the shared env dir and is not ours to edit unilaterally.
  `tool-creator` scaffolds a `modal_image.py`, not an activation script.
- **A failed `sapia run` leaves a phantom table in `_registry.tsv`.** Verified in
  `core/base_run.py`: inside the `DataManager` context the driver calls
  `resolve_output_table()` (which does `registry.register_table()`), then creates `out_dir`
  and `<leaf>_logs/`, then writes `.meta.json` — and only *after* all that calls
  `build_manifest_fn`. So any error a tool raises while building its manifest (rpxdock's
  input guard, a bad `--allowed-residues`, a missing column) arrives too late: the registry
  row and the directories are already committed, and no `<table>.tsv` is ever written.
  **A registry row with no matching `.tsv` is a refused or abandoned submit, not a table
  with zero rows** — remember that when auditing lineage. A tool cannot clean this up
  itself: `ManifestCtx` carries no handle on the registry. Fixing it properly means
  reordering upstream (build the manifest before reserving the table, or roll back the
  registration on exception). Until then, delete the leaf dir and the registry row by hand.
- **Nothing syncs the two backends.** No command moves a `run_dir` between the Modal Volume
  and cluster storage; a campaign started on one finishes on that one.
- **Untested on vib:** the full chain, and every custom tool by definition.

## Network note (imec)

This network intercepts TLS to Modal's blob storage, so `modal volume get` / `put` can fail
with a certificate error while `modal volume ls` works. Because data stays on the Volume and
`sapia` runs in the workstation, normal work is unaffected — only local up/downloads are.
Do not work around it by weakening TLS verification; move the data through the workstation
or raise it with IT.
