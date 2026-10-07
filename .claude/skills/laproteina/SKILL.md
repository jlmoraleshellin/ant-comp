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

> **Status: the code is tested, a real run is not.** The manifest builder and the
> collector have unit tests against a mock output tree (3 jobs → 3 tasks with
> distinct seeds; 2 jobs × 2 samples → 4 rows, rewards joined, empty job handled).
> Nothing has yet run on Modal. Treat the first run as a test and start with
> `--num-jobs 1 --nsamples 2`.

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

Columns (leaf-prefixed `laproteina_`): `job`, `length`, `sample`, `sample_dir`,
`task_name`, plus `_status` and `_path` (the complex PDB). If the job wrote a
`rewards_*.csv`, any column whose name contains reward/score/plddt/ptm/iptm/pae/rmsd
is lifted in alongside.

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

Confirm the binder's chain id from a collected PDB before composing the atomium
run; LPC writes the target chains it was given plus the binder, and the binder's
letter depends on the target's.
