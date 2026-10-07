---
name: laproteina
description: How to run the custom laproteina tool — La-Proteina Complexa (LPC), NVIDIA's fully atomistic flow-matching binder generator. Covers why it is root-only (the target comes from LPC's own registry, not a table column), the --num-jobs x --nsamples cost model, why only the generate stage is run and what that gives up, the two Volumes it needs (checkpoints and the target registry), the silent-skip trap in LPC's generate.py, and the columns it collects. Load before composing a laproteina run or interpreting its columns.
---

# laproteina

**Custom tool** (lives in `tools/laproteina/`, not bundled with prosapia).
**`action: create`** — mints a child table, one row per generated backbone.
**Root-only**: it takes no `-t/--table`, and the builder raises if you pass one.

Upstream: [NVIDIA-BioNeMo/Proteina-Complexa](https://github.com/NVIDIA-BioNeMo/Proteina-Complexa),
default branch `dev`. A partially latent flow-matching model that samples backbone
geometry, side chains and sequence **jointly**, conditioned on a target structure
and a hotspot list. NVIDIA Open Model License — commercially usable, outputs not
claimed by NVIDIA.

> **Status: validated end to end on Modal, 2026-10-07.** One job, 2 designs
> against `90_TNFa_HUMAN_xreact`: exit 0, both rows `OK`, **~2 min on one A100**
> including model load. Reference run `outputs/20261007_123651_tnf_smoke` on
> `sapia-runs-luca`. Everything below marked "measured" comes from that run.

## Verified invocation

```bash
sapia new_run --label lpc_tnf
sapia run laproteina <run_dir> --task-name 90_TNFa_HUMAN_xreact \
    --num-jobs 1 --nsamples 2 --design-prefix tnfh
sapia collect laproteina <run_dir> -t <the table the run reserved>
```

## The target is a name, not a path

This is the thing that surprises people. LPC resolves its target out of **its own
registry** — `configs/targets/targets_dict.yaml` inside the LPC checkout — which
carries the structure path, the chain/residue ranges the model sees
(`target_input`), the hotspots and the binder-length range.

So there is no input column, nothing for `-t` to supply, and `--hotspots` does not
exist as a flag here. To design against a new epitope you **edit the registry entry**,
not the `sapia` command line. `complexa target list` shows what is registered.

That also means `{expr}` placeholders are meaningless for this tool — a root create
has no lineage to resolve against, and the values live in LPC's config anyway.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--task-name` | **required** | The key in LPC's `targets_dict.yaml`. |
| `--num-jobs` | 1 | Array tasks to **split** the run across. Parallelism, not extra work. Also sets LPC's `gen_njobs`, which is mandatory — see the trap below. |
| `--nsamples` | config's (4) | **Total** binders for the whole run, divided across the jobs. `--num-jobs 4 --nsamples 20` is 20 designs, 5 per job — not 80. |
| `--seed` | config's | Job *j* runs with `seed + j`, so jobs don't redraw the same binders. Leave unset and the jobs share a seed, separated only by `--job-id`. |
| `--design-prefix` | `--task-name` | Prefix for every row name. The default is descriptive but long; `tnfh` beats `90_TNFa_HUMAN_xreact`. |
| `--config` | `configs/search_binder_local_pipeline.yaml` | Only its `generation.*` section is ever used. |
| `--set KEY=VALUE` | none | Raw Hydra override, `++`-prefixed for you. Repeatable. The escape hatch. |

Default Modal resources: **A100**, 8 CPU, 32 GiB, 4 h timeout — a guard, not an
estimate. Re-size with `-T/--time` once a real run has been measured.

## Only the `generate` stage runs — and what that gives up

LPC's own `complexa design` is a four-stage pipeline: **generate → filter →
evaluate → analyze**. This tool runs **generate alone**.

What that skips, deliberately:

- its own **SolubleMPNN** sequence design — because the campaign designs sequences
  with `atomium` instead, which is where tied-positions multi-target design happens;
- its own **AF2 / ESMFold / RF3 refolding** — because `boltz` is the independent
  validator here, and LPC's AF2 scores are not independent anyway when AF2 also
  steered generation;
- **inference-time search** (beam search, MCTS, best-of-N) driven by reward models.

That last one is the real cost. Search is what the upstream paper is named after
and where much of its reported quality comes from. The `--minimal` image omits JAX,
ColabFold and tmol, so the reward models are not even installed. **This is the cheap
mode.** If the budget allows turning search back on later, that is a different image
and a different cost model, not a flag.

Same division of labour as `bindcraft2 --trajectory-only`.

## Two Volumes, and why the second one exists

| Mount | Volume (`.env` key, default) | Holds |
| --- | --- | --- |
| `/ckpts` | `SAPIA_MODAL_VOLUME_LPC_CKPT`, `lpc-checkpoints` | `complexa.ckpt` (2.7 GB) + `complexa_ae.ckpt` (3.8 GB) |
| `/lpc_assets` | `SAPIA_MODAL_VOLUME_LPC_ASSETS`, `lpc-assets` | `configs/targets/targets_dict.yaml` + `data/target_data/` |

The checkpoints are public on NGC — plain `wget`, **no account and no API key**
(verified: the endpoint 302s to a signed URL and returns 200 unauthenticated). But
6.5 GB per container is not something to re-download per task, so populate the
Volume once.

The assets Volume is the less obvious one. **The Modal image clones LPC fresh from
GitHub, so a target added locally does not exist in it.** The task script overlays
`/lpc_assets` onto the checkout before running; without it, `--task-name
90_TNFa_HUMAN_xreact` dies at config resolution in every container. Whenever the
epitope or the target structure changes, the assets Volume must be updated too —
editing the local checkout alone changes nothing on Modal.

## Traps

### A stale results CSV makes generate exit 0 having done nothing

`generate.py` returns immediately if `./inference/results_<config_name>_<job_id>.csv`
exists:

```
Results already exist at ... Exiting generate.py.
```

That filename is keyed by **job id, not run name**, so a later run reusing the same
`--job-id` would silently produce nothing and still exit 0. `laproteina.sh` deletes
this job's marker before running. If you ever invoke `complexa generate` by hand in
the same checkout, know that this is why a job "ran" in two seconds.

### The collected structure is the complex, not the binder

`laproteina_path` is the **binder plus the target chains LPC was conditioned on**.
Anything that needs the binder alone — self-consistency against its own backbone,
an interface metric — must run `chainsel` first. Same contract as `bindcraft2_path`.

### `gen_njobs` must match the task count, or every job but the first does nothing

`generate.py` reads `njobs = cfg.get("gen_njobs", 1)` and hands it to
`split_by_job(cfg, job_id, njobs)`, which slices the sample budget: job *j* takes
`ceil(nsamples/njobs)` samples starting at `j * that`. If `gen_njobs` stays at the
config's default of **1** while tasks run with `--job-id 1, 2, …`, the slice check
`nsamples_per_split * job_id >= nsamples` is true immediately and the job logs

```
Job id 2 get 0 samples. Finishing job...
```

and calls `sys.exit(0)` — **exit 0, no structures, nothing to collect**. The tool
sets `gen_njobs` from `--num-jobs` for exactly this reason. Never override it by
hand with `--set`.

A related consequence: `--num-jobs` greater than `--nsamples` leaves the surplus
tasks with empty slices. The builder warns; prefer fewer, fuller jobs.

### The nsamples override path is `nres.nsamples`

The count lives on the `nres` sampler
(`generation.dataloader.dataset.nres.nsamples`), not on `dataset` itself. Overriding
`dataset.nsamples` with Hydra's `++` quietly **adds a dead key** and leaves the real
value at its default of 4 — a run that looks configured and is not.

### Row names get long fast

Without `--design-prefix`, a row is named
`90_TNFa_HUMAN_xreact_j0_n_87_id_0`. Pass a short prefix.

## What it collects

Child rows named `<group>_n_<length>_id_<idx>`, each linked to its job group.

Columns (leaf-prefixed `laproteina_`): `job`, `sample`, `sample_dir`,
`task_name`, `lpc_n`, plus the chain shape measured from the structure —
`chains`, `n_chains`, `binder_chain`, `binder_length`, `target_length`,
`complex_length` — plus `_status` and `_path` (the complex PDB). If the job wrote
a `rewards_*.csv`, any column containing reward/score/plddt/ptm/iptm/pae/rmsd is
lifted in alongside.

**Filter on `binder_length`, never on `lpc_n`.** LPC's own `n_<N>` in the sample
directory name is the **complex** length, target included. Measured: a 1TNF
trimer target gave `lpc_n` 519 and 535 for binders of **63 and 79** residues
(456 of target + the binder). A length filter against `lpc_n` would be comparing
a binder-length threshold to a number four hundred residues too big.

**The structures are the contract; the rewards CSV is a bonus.** A job that wrote
structures but no CSV still collects cleanly, with the reward columns simply absent.

### Ranking

LPC's own reward columns — when present — come from the model that generated the
designs, so they are **not** independent validation, for the same reason bindcraft2's
confidence scores are not. Rank by what comes after: `boltz` on the complex, and
self-consistency of the binder back to its generated backbone.

## Wiring to the next step

| Going to | Pass |
| --- | --- |
| `atomium` / `proteinmpnn` | `-i laproteina_path` — the backbone is the **complex**, so also pass `--chains-to-design <binder chain>` |
| `chainsel` / `cms` / `usalign` | `-i laproteina_path`, and `chainsel` the binder out first |
| `boltz` | not directly — design sequences first |

**Measured on the TNF trimer: the binder is chain `D`.** LPC writes the target
chains it was conditioned on followed by the binder, so the binder is the last
chain — A/B/C of 152 residues each (the trimer) then D. The collector records
this per design in `binder_chain`, so read it from the table rather than
assuming: a target with a different chain count moves the letter.

## Image pins — do not loosen these

LPC installs several dependencies unpinned, and three of them have drifted since
its 1.0.0 release. Each one breaks the build or the first import, and each pin
below exists because it cost a cycle:

| Pin | Why |
| --- | --- |
| `atomworks==2.2.1` | 3.0.0 (2026-10-05) moved `AtomSelectionStack` from `io/utils/selection.py` to `io/utils/query.py`; LPC's `pdb_utils.py:51` still imports the old path. |
| `dm-haiku==0.0.13` | 0.0.14+ calls `jax.core.take_current_trace`, absent from the jax 0.4.29 this stack pins. Comes in unpinned via colabdesign. |
| `libxrender1`, `libxext6` | From LPC's own Dockerfile. openbabel `dlopen`s its format plugins; without libXrender **every** plugin fails, and openbabel then mis-parses its own error output and raises `ValueError: not enough values to unpack`. Surfaces as a Hydra `Error locating target 'gen_dataset.collate_fn'` — a missing system library, three frames from where it shows. |

Two things **not** to add:

- **`rc-foundry[all]`** (build_uv_env.sh step 8) pulls torch 2.14, replacing the
  pinned 2.7.0; `torch_cluster`/`torch_geometric` are compiled against 2.7.0 and
  the container dies with a bare `Bus error`. It is not needed.
- **`--minimal`.** The flag suggests JAX and ColabDesign are optional. They are
  not: `search/__init__.py` imports every algorithm eagerly, including
  `sequence_hallucination`, whose module body does `import jax` and
  `from colabdesign import mk_afdesign_model`. The import fires even with
  `search.algorithm=single-pass` and `reward_model=null`. The **libraries** are
  required; the AF2 **weights** are not.

The image ends with two import checks (`proteinfoundation.search`, and
`gen_dataset.collate_fn`). Keep them. Both failures above first appeared on a GPU
task; as build steps they fail in minutes on a CPU builder instead.

## Submit trap: "Submitting N task(s)" is printed even when the build failed

prosapia prints

```
Submitting 1 task(s) to Modal app 'sapia-laproteina' (gpu=A100, ...)
```

**after** a `modal.exception.ImageBuildError`, because the two come from
different code paths and stdout/stderr interleave. Do not treat that line as
proof a task exists. Check for `ImageBuildError` in the submit output, or that
the log dir's `.exit` is newer than the submit. Cost 30 minutes of polling a task
that was never created — twice.
