---
name: bindcraft2
description: How to run the custom bindcraft2 tool on Modal — BindCraft2 binder-design campaigns, where one task is a whole campaign rather than one design. Covers the campaign cost model (--num-designs vs --max-trajectories), root vs child runs, hotspots with {expr}, --targets for TRUE multi-specificity (one binder optimised jointly against several targets, with detarget counter-selection) and the per-target column layout it collects, --trajectory-only for handing backbones to atomium/proteinmpnn instead of BC2's own MPNN, --reuse-campaigns for collecting one campaign into two tables (how to compare BC2's own designs against atomium's on the same backbones), the table/dir labels that fork needs and the silent table collision when you omit them, the AlphaFold-parameter cache Volume, the three collect stages, --copies for a homo-oligomeric (dimeric, trimeric) binder and why its lengths are per copy, how to get the binder's chain LETTER, the secondary-structure presets (mixed_topology is the ANTI-helix one), and why its own confidence scores are not independent validation. Load before composing a bindcraft2 run or reading its columns.
---

# bindcraft2

**Custom tool** (lives in `tools/bindcraft2/`, not bundled with prosapia).
**`action: create`** — mints a child table, one row per binder design a campaign produced.

Upstream: [PacesaLab/BindCraft2](https://github.com/PacesaLab/BindCraft2), pinned to tag
`v1.0.3` in `modal_image.py`. Its `docs/source/reference.md` is authoritative for
settings, `docs/source/outputs.md` for columns.

> **Not yet run end to end.** Everything below is read off the upstream source and this
> tool's own self-test (`uv run python tools/bindcraft2/test_bindcraft2.py` — manifest,
> settings files and the collector against fabricated campaign folders, no GPU). The
> first real submit should be one small campaign, not a fan-out — see *First run* below.
> **The multi-target path has never been run against real BindCraft2 output**; the
> per-target layout is derived from `campaign_output.py` at the pinned tag.

## The thing to understand first: a task is a campaign

Every other design tool here maps one design to one unit of work. BindCraft2 does not.
One `bindcraft design` process takes **one target** and loops: hallucinate a backbone
with AF2, redesign its sequence with ProteinMPNN, refold it, filter it, keep or discard —
until `--num-designs` have been **accepted** or `--max-trajectories` attempts are spent.

Consequences that shape how you drive it:

- **One manifest row = one campaign = one target.** `-C/--max-concurrent` throttles
  targets, not designs.
- **The row count is unknown until collect.** A campaign can accept fewer designs than
  you asked for, or none at all. That is a normal outcome, not a failure.
- **`--num-designs` is a stopping condition, not a batch size.** The real cost knob is
  `--max-trajectories`: a campaign runs until one of the two fires.
- **It already contains the rfd3 → mpnn → predict → filter chain.** Do not compose those
  tools around it; compose *independent* checks after it (see the last section).

## Two shapes of run

- **Root run (no `-t`)** — starts a fresh lineage in `table0`, one campaign. The target is
  `--target-pdb <file>` (design group `<stem>_bc2`) or `--shipped-target <name>` (group
  `<name>_bc2`), for a target BindCraft2 ships: `hPDL1`, `hPD1`, `mPDL1`, `hIL2R`,
  `hIL7RA`, `dynorphin_a`.
- **Child run (`-t <table>`)** — one campaign per ready row, targets from
  `-i/--input-column` (default `pdb_path` — set it to whatever column actually holds your
  target, e.g. `mkcomplex_path`, `chainsel_path`). Mints a child table.
- **Multi-target run (`--targets <file>`)** — one campaign, several targets, **one
  binder** optimised against all of them. A different experiment from either of the
  above; see the next section.

`{expr}` placeholders resolve up the lineage, so they need `-t`. A root run must use
literal values, and says so rather than resolving to nonsense.

## Several targets, one binder: `--targets`

This is the flag for "design a binder that binds **both** human and mouse EGFR", or "binds EGFR and **not** ERBB2". It is not a loop and it is not a post-hoc merge of two campaigns: BindCraft2 v1.0.3 puts one loss instance per (loss × target) into a single weighted sum (`loss.py`), `multitarget_merged_gradients` (default true) refuses a sequence update until *every* target has contributed a gradient, and `multitarget_tied_redesign` (default true) ties the ProteinMPNN redesign across targets. The binder that comes out was selected to satisfy all of them at once.

`--targets` takes a YAML/JSON file — a list of target mappings, or `{targets: [...]}` so you can lift the block straight out of upstream's `examples/pdl1_crossreactive_detarget.json`. Each entry carries **only** BindCraft2's seven per-target keys (`settings.py:38`): `name`, `target_path`, `chains`, `hotspots`, `coldspots`, `weight`, `objective`.

```yaml
# targets/egfr_pair.yaml
- name: hEGFR
  target_path: targets/hEGFR_d3.pdb
  chains: A
  hotspots: A355,A356,A440,A441
  weight: 1.0
- name: mEGFR
  target_path: targets/mEGFR_d3.pdb
  chains: A
  hotspots: A355,A356,A440,A441
  weight: 1.0
- name: hERBB2          # off-target: counter-select AGAINST this one
  target_path: targets/hERBB2.pdb
  chains: A
  weight: -0.5
```

```bash
sapia run bindcraft2 <run_dir> --targets targets/egfr_pair.yaml \
    --binder-lengths 60-100 --num-designs 10 --max-trajectories 300
sapia collect bindcraft2 <run_dir> -t <the table the run reserved>
```

**Counter-selection.** A negative `weight` (or `objective: detarget`, which forces the weight negative) makes a target something the campaign is pushed *away* from. Upstream excludes detargets from its own ranking, and so does this tool — see the column layout below.

**Rules this tool enforces at submit time**, so a typo costs a second and not a container-hour:

| Refused | Because |
| --- | --- |
| an unknown per-target key | upstream rejects it too, but only after the image is built. You get a `did you mean 'hotspots'?` here. |
| `--targets` with `--target-pdb` / `--shipped-target` | three ways to name a target, one campaign. |
| `--targets` with `--chains` / `--hotspots` / `--coldspots` | **in this mode those are per-target sub-keys.** Move them into the entry of the target they describe. Silently applying one epitope to all three targets is exactly the plausible-looking wrong answer to avoid. |
| a `name` with `.`, `;`, `/` or a space | the name is a metric suffix (`i_pTM.hEGFR`), a filename suffix, an entry in the `;`-joined `targets` cell, and a table-column fragment. |
| a `name` starting with `off_` | the collector uses that infix to mark an off-target's columns. |
| duplicate names, a bad `objective`, a non-numeric `weight` | |
| a list of nothing but detargets | the campaign would have nothing to bind. |
| a `target_path` that does not exist | checked per campaign, **after** `{expr}` resolution. |

**Relative `target_path`** is read against the submit cwd first (i.e. `/runs`, like every other path here), then against the `--targets` file's own directory — so a self-contained targets folder also works. Both attempts are named when neither exists, and the chosen path is written into the settings file absolute.

**Design-group naming.** Without `-t`, the group is `<file stem>_bc2` (`egfr_pair.yaml` → `egfr_pair_bc2`), matching the `--target-pdb` convention. With `-t`, the group is the **row name**, one campaign per ready row, all against the same target list — and the row's `--input-column` is then *not* a target; it only decides which rows are ready and supplies the lineage `{expr}` resolves against. The tool prints that in so many words when it happens.

**Shipped targets can go multi too** — upstream accepts `"target": ["hPDL1", "mPDL1"]` and accumulates each preset's own `targets` block (`settings.py:requested_preset_names`). `--shipped-target` here takes **one** name only; the pair form is not wired up. Use `--targets` with explicit paths instead.

## Invocation

```bash
# root: one campaign against a target not in any table yet
sapia run bindcraft2 <run_dir> \
    --target-pdb targets/PDL1.pdb --chains A \
    --hotspots 'A54,A56,A66,A115' \
    --binder-lengths 60-100 --num-designs 10 --max-trajectories 200

# child: one campaign per target row, epitope read from the table
sapia run bindcraft2 <run_dir> -t table0 -i pdb_path \
    --hotspots 'A{epitope_start}-{epitope_end}' --num-designs 5

sapia collect bindcraft2 <run_dir> -t <the table the run reserved>
```

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--targets FILE` | none | YAML/JSON list of targets **one** binder is optimised against jointly. Negative `weight` = counter-select. Excludes `--target-pdb`/`--shipped-target` and the top-level `--chains`/`--hotspots`/`--coldspots`. See above. |
| `--trajectory-only` | off | Stop at **backbones** — no MPNN, no refold, no filter, nothing accepted. The handoff to this workspace's own designers; see below. |
| `--reuse-campaigns` | none | `TABLE[:LABEL]`. Submit nothing; reserve a second table over an earlier run's campaigns, to collect a second stage of the same GPU hours. See below. |
| `--max-trajectories` | BC2's own | **The cost knob.** Attempts spent before giving up. Set it or a hard target can burn the whole timeout. With `--trajectory-only` it is the *only* budget (BC2 defaults it to 100). |
| `--num-designs` | BC2's own | Accepted designs to stop at (`number_of_final_designs`). A target, not a guarantee. |
| `--hotspots` | none | Target residues the binder should contact, BC2 syntax (`A54,A56,A66-70`). Takes `{expr}`. Omit to let BC2 pick the epitope. |
| `--coldspots` | none | Regions to avoid, same syntax. |
| `--chains` | all | Target chains to design against (`A`, `A,B`). |
| `--binder-lengths` | modality's | `80` or `60-100`. Takes `{expr}`. **Per copy** when `--copies` is above 1. |
| `--copies` | BC2's own (1) | Identical chains the binder is built from — `3` for a homo-trimer. Above 1 switches on BC2's multi-chain-binder feature. See below. |
| `--modality` | `binder` | `binder`, `VHH`, `peptide`, `cyclic_peptide`, `ARP`, `scFv`, `Fab`, `large_binder`, `homo_oligomer`, `multidomain`, `induced_fit`, `fold_switch`. Comma-separated to combine. |
| `--property` | none | Repeatable preset: `humanize`, `protease_stable`, `disulfide_staple`, `forced_targeting`, `initial_guess`, `mixed_topology`, `termini_together`, `termini_accessible`, `bigbang`. Validated at submit time. **Read the secondary-structure note below before reaching for `mixed_topology`.** |
| `--core` | none | Core profile under every preset; `benchmark` for a reproducible run. Pair with `--campaign-seed`. |
| `--extra-settings` | none | YAML/JSON merged into every campaign's settings (`objective`, `aa_bias`, `min_iptm_final`, `save_design_trajectory`, …). Values take `{expr}`. |
| `--set KEY=VALUE` | none | Verbatim `bindcraft design --set`, applied over the generated file. **No spaces** — the task script word-splits it. Use `--extra-settings` for anything richer. |
| `--no-resume` | off | See below. |

`bindcraft design --list-settings` names every setting `--extra-settings` and `--set`
accept. This tool writes the settings file itself; a key set by both a flag and
`--extra-settings` is refused at submit time rather than silently resolved.

Default Modal resources: **A100**, 8 CPU, 32 GiB, **12 h** timeout. Size the timeout to
the campaign with `-T`, and remember a killed container leaves no `.exit` file.

### A homo-oligomeric binder: `--copies`

`--copies N` writes BindCraft2's `copies` campaign setting, the number of **identical** chains the binder is built from. It is independent of the target's chain count, which is `--chains`. A trimeric target and a trimeric binder are two separate decisions, and either can be used without the other.

```bash
# a homo-trimeric binder against a trimeric target
sapia run bindcraft2 <run_dir> --target-pdb targets/trimer.pdb \
    --chains A,B,C --copies 3 --modality homo_oligomer \
    --binder-lengths 50-80 --hotspots 'A54,B210-214'
```

Above 1, upstream switches on its *multi-chain binder* feature (`settings.py`, `CAMPAIGN_FEATURES`): validation moves to the **multimer** model, protomer-scoped losses go **per protomer**, and an `Oligomer_Symmetry_RMSD` check is installed. The copies are tied symmetrically unless you pass `--set oligomer_tie=none`.

Four things to get right:

- **`--binder-lengths` is per copy.** `--copies 3 --binder-lengths 50-80` is a 150–240 residue assembly, not an 80-residue one split three ways.
- **`--copies` beats `--modality homo_oligomer`.** That preset sets `copies: 2` itself, so without the flag you get a **dimer**. Upstream layers the campaign settings file over every preset (`campaign_over_presets`), and this tool writes `copies` into that file — so pass both to get the oligomer filters *and* your own chain count.
- **The binder is several chain letters.** Behind a 3-chain target a trimeric binder is `DEF`. `bindcraft2_binder_chain` carries them all, but **only** when `bindcraft2_binder_chain_src` reads `stamp`. The `last_chain` fallback returns one letter and is an under-count — and that is exactly the stage `--trajectory-only` lands in, since upstream stamps accepted designs only. Check the `_src` column before passing letters to `chainsel`, `cms` or `ssprofile`.
- **It cannot be combined with a multi-chain `binder_scaffold`.** Upstream refuses it: a multi-chain scaffold is one binder spanning its chains, not copies of one.

A target chain count has no such flag and needs none — the wrapper never opens the target structure, so a trimer, a pentamer or a fused receptor passes through on `--chains` alone.

### Secondary structure: `mixed_topology` is the ANTI-helix preset

**Do not use `--property mixed_topology` to ask for helical binders.** Its own description in `settings/property/mixed_topology.json` is *"At most 50% helix and at least 20% beta sheet"*. It does three things, all of them away from helix:

```jsonc
{ "weights_non_helical": 1.0,          // rewards NON-helical contacts
  "weights_binder_helicity": 0,        // zeroes the default helicity reward
  "max_helix_fraction_final": 0.5,     // rejects anything over 50% helix
  "filters": { "Binder_BetaSheet_Fraction": { "threshold": 0.2, "higher": true } } }
```

The knob you actually want is **`weights_binder_helicity`, which defaults to `-0.3`** (`settings/core/default.json`). BindCraft2 minimises a weighted sum of losses, so a **negative weight favours helix** and a more negative one favours it harder. There is no dedicated flag; reach it through `--extra-settings`:

```yaml
# helical.yaml
weights_binder_helicity: -0.8     # more helical than the -0.3 default
```

```bash
sapia run bindcraft2 <run_dir> --targets targets/egfr_pair.yaml \
    --extra-settings helical.yaml
```

So: BC2's **default is already mildly helix-favouring**; `mixed_topology` turns that off and pushes the other way. If you want helix, either leave the default alone or make `weights_binder_helicity` more negative — and verify the outcome with `ssprofile` rather than trusting the weight, since a loss weight is a preference, not a constraint.

### `resume` is on by default

This tool writes `"resume": true` into every settings file, so re-submitting continues
the campaigns already started rather than refusing or restarting. It has to: the
framework's own resume filter reads `<leaf>_status`, which for a `create` tool lives in
the *child* table, not the one the run reads — so `sapia run` cannot skip finished
campaigns on its own. `--no-resume` forces a clean start.

## Backbones only: `--trajectory-only`

BC2's ProteinMPNN is **not** swappable by configuration — it is a JAX/Haiku
reimplementation inside the package, `jit`-compiled, fed in-memory arrays. `mpnn_model`
and `mpnn_variant` only choose among `.npz` checkpoints of that same architecture
(`neutral` / `negative` / `positive`). Nothing reaches a PyTorch model.

What BC2 *does* ship is a clean seam. With `--trajectory-only` it builds the backbone
and stops — `mpnn_model` is never constructed, no design is ever accepted — so you can
redesign the sequences with `atomium` or `proteinmpnn` and finish the campaign with this
workspace's own tools:

```bash
sapia run bindcraft2 <run_dir> -t table0 -l traj \
    --trajectory-only --max-trajectories 40 --hotspots 'A54,A56'
sapia collect bindcraft2 <run_dir> -t <the table the run reserved> -l traj
```

`sapia collect` needs no `--stage`: the run records `trajectory_only` in the out_dir
sidecar and collect follows it, the same contract that keeps `input_column` consistent
across the two phases.

### Two things to get right

- **Filter on `bindcraft2_completed`.** `!_Trajectories.csv` has one row per attempt,
  *terminated ones included*. Those are collected, not silently dropped — a terminated
  attempt that wrote no structure yields no row, but one that died late still can. BC2's
  own `terminated` column is blank when the attempt *succeeded*, which inverts prosapia's
  "blank means not applicable", so the collector adds an explicit boolean. Filter on it
  before spending sequence design on junk backbones.
- **Use `bindcraft2_binder_chain` for the chain letter — do not read `binder_chain`.**
  The collected structure is the **complex**, and the binder's chain *letter* is now a
  collected column. See *The binder's chain letter* below for why the setting of the same
  name is the wrong thing to read. Pass that letter to atomium's `--chains-to-design`. Do
  *not* `chainsel` the binder out first here — you want it redesigned in the target's
  context. (`chainsel` comes later, before `usalign`/`cms`, which do want the binder
  alone.)

Table naming and labels for this fork are in *Two tables from one campaign* below.

### What you give up

BC2's redesign call is not a plain MPNN pass. `mark_redesign_residues` in `MPNN_stage.py`
holds interface contacts fixed when `redesign_interface` is set, ties chains for
oligomers, masks fold-conditioning scaffold and inter-domain linkers, and applies the
campaign's `aa_bias`. A vanilla `atomium` run reproduces none of that unless you rebuild
it with `--fixed-positions`. That is a real difference in *what is being designed*, not
just in which model designs it — decide it deliberately rather than discovering it in the
results.

You also lose the whole filter battery, which is the point: you are replacing it with
`boltz` + `usalign` + `cms`, which are independent of the AF2 that built the backbone.

## Two tables from one campaign, and how they are named

One run reserves **one** table, and its output dir is derived from it
(`run_dir/<the table it reserved>/bindcraft2[_label]/`). So a second
`sapia collect -t <other table>` looks in a directory nothing ever wrote to and finds
nothing — even though the campaign it wants is right there, since one campaign folder
holds all three stages.

That bites the moment you want to **compare BindCraft2's own designs against atomium's
on the same backbones**. Both branches have to come out of one campaign, or they are
built on different trajectories and the comparison is confounded.

`--reuse-campaigns` closes it: a run that submits nothing, reserves a table, and records
in its sidecar that its campaigns live in an earlier run's out_dir. Collect follows the
record instead of assuming the two coincide.

```bash
# table0 already holds the target rows. One FULL campaign per ready row.

# 1. the campaigns, reserving table1. Collect their backbones.
sapia run bindcraft2 R -t table0 -l traj --hotspots 'A54,A56' --num-designs 10
sapia collect bindcraft2 R -t table1 -l traj --stage trajectories

# 2. the SAME campaigns' accepted designs, into table1_bc2. Submits nothing, no GPU.
sapia run bindcraft2 R -t table0 -l bc2 --table-label bc2 --reuse-campaigns table1:traj
sapia collect bindcraft2 R -t table1_bc2 -l bc2 --stage ranked

# 3. the atomium branch off the BACKBONES table, reserving table2
sapia run atomium R -t table1 -i bindcraft2_traj_path --chains-to-design <binder chain>
```

**Two different tables are named here, and confusing them is the trap.** `-t` is the
table the run *reads ready rows from* — the reuse run passes the **same `-t` the original
run used** (here `table0`), because the collector maps each ready design to a campaign
folder by name and those folders are named after the *original* run's design groups.
The `--reuse-campaigns` argument is a different thing: `TABLE[:LABEL]` where **TABLE is
the table the earlier run COLLECTED INTO** (`table1`) and **LABEL is that earlier run's
`-l/--dir-label`** (`traj`) — *not* its `--table-label`. A root campaign needs no `-t` at
all; the reuse run reads the group names off disk.

**The reuse run needs its own `-l`** (`-l bc2` above), or its columns land as plain
`bindcraft2_*` and the `bindcraft2_traj_hash` / `bindcraft2_bc2_hash` join described
below does not exist.

> **Corrected 2026-10-06.** This example previously showed a root run (no `-t`) while its
> own prose and the lineage diagram described a child run, and it passed
> `--reuse-campaigns table1:traj` against a run that had reserved `table0`. It also
> omitted `-l bc2` from the reuse run, and branched atomium off `table0` (the targets)
> rather than `table1` (the backbones). Verified against `sapia run bindcraft2 --help`
> and `collect_bindcraft2.py`.

### The two labels are different things, and one of them is mandatory

| Flag | Names | Omit it and… |
| --- | --- | --- |
| `-l/--dir-label` | the out_dir leaf and the column prefix (`bindcraft2_traj_path`) | both tables' columns are `bindcraft2_*`, with no way to tell a backbone from a validated complex |
| `--table-label` | the table itself (`table1_bc2` instead of `table1`) | **the second run silently reserves the first run's table** |

That second one is the trap, and it is silent. `derive_new_table` names a child
`table<gen>[_<label>]` from the parent's generation, so two unlabelled runs off `table0`
both derive `table1` — and `register_table` treats identical lineage as idempotent rather
than as a conflict. No warning; the two runs just share a table. **Give every fork a
`--table-label`.** Both labels have to match between `run` and `collect`.

Note the resulting name: a child of `table0` labelled `bc2` is `table1_bc2`, not `table2`.
`table2` is gen 2 — what the atomium run off `table1` reserves.

### The lineage is a fork, joined on `hash`

```
table0  targets
  ├ table1      bindcraft2 backbones (-l traj) ──→ table2 atomium ──→ boltz / usalign / cms
  └ table1_bc2  bindcraft2 ranked    (-l bc2)      BindCraft2's own answer
```

Each branch is an ordinary generation — `create` mints a child, `lookup` walks up, so
`{expr}` and ancestor columns still resolve from the bottom. Nothing about lineage needs
bending; there is simply one more generation than a plain campaign has.

The two bindcraft2 tables are **siblings**, not ancestor and descendant, so pair them
with BindCraft2's own join key rather than with lineage: `bindcraft2_traj_hash` on a
backbone equals `bindcraft2_bc2_hash` on every accepted design that came from it. One
backbone can yield several (`_seq0`, `_seq1`, … up to `kept_sequences`), so it is a
one-to-many join.

Compare like with like: `table1_bc2` carries BC2's own AF2 numbers, the atomium branch
carries Boltz's. The honest comparison is Boltz-on-both — predict `table1_bc2`'s
sequences too (`-i bindcraft2_bc2_sequence`) rather than reading BC2's `i_pTM` against
Boltz's.

## The AlphaFold parameters are a Volume, and they are not warm

~5.3 GB, downloaded on first use into `$XDG_CACHE_HOME/bindcraft` — mounted at
`/bindcraft_cache` from the `bindcraft-cache` Volume
(`SAPIA_MODAL_VOLUME_BINDCRAFT_CACHE`). ProteinMPNN's weights ship inside the package, so
these are the only download.

**Warm it with one campaign before fanning out.** Submit N targets cold and N containers
each pull the same 5.3 GB. `BINDCRAFT_AF2_PARAMS` in `.env` overrides the location if the
weights already live somewhere shared.

Per the workspace convention the Volume is named without the `sapia-` prefix so it can be
shared. **Never delete it to tidy up** — it re-downloads, slowly.

## First run

The image builds inside the submit call: Ubuntu + jax cuda13 wheels + an editable install
of the repo. Allow **15+ minutes** before assuming a first submit is stuck; later runs are
seconds. Then:

1. One target, `--max-trajectories 10` or so, and watch it finish. This warms the weights
   Volume and proves the GPU path.
2. Check `campaigns/<name>/1_Trajectories/!_Trajectories.csv` appears. If jax fell back to
   the CPU the campaign still *runs*, ~100× slower, and says so only in a buried warning —
   that is the failure mode the image's `ldconfig` step exists to prevent.
3. Only then fan out.

Results are written incrementally, so you can read `3_Ranked/!_Ranked.csv` while a
campaign is still going.

## What it collects

Child rows keyed by BindCraft2's own design identity (campaign, modality, length, recipe
hash, `_seq<n>`), each linked to its parent target. **Not** by rank: `!_Ranked.csv` is
re-sorted by `i_pDAE` after every acceptance, so a positional key would not survive a
re-collect.

`--stage` picks what to collect, and defaults to **`auto`** — what the run recorded in
its sidecar, so you normally never pass it:

- **`trajectories`** (auto after `--trajectory-only`) — `1_Trajectories/!_Trajectories.csv`
  plus `1_Trajectories/<design>/<design>_trajectory.cif`, the hallucinated backbones.
  Adds `terminated`, `autotuned`, `length` and the derived `completed` boolean.
- **`refolded`** — `2_Refolded/!_Refolded.csv`, every scored candidate including rejects,
  with `outcome` and `failed_filters`. **This is the one to reach for when a campaign
  accepted nothing** and the question is why. Collecting it *as well as* another stage
  needs a `--reuse-campaigns` run, not just a second `-l` — see below.
- **`ranked`** (auto otherwise) — `3_Ranked/!_Ranked.csv`, only what the campaign accepted.

All three carry the same metric battery: BC2 scores a trajectory on the same filters it
runs at refold, so the columns are comparable across stages.

If `archive_trajectories` was on, the per-design folders are zipped and the collector says
so — run `bindcraft unarchive <campaign folder>` first.

Columns (leaf-prefixed `bindcraft2_`): `sequence`, `i_pDAE`, `i_pTM`, `i_pAE`, `pLDDT`,
`pTM`, `Unbound_Binder_pLDDT`, `Interface_Residues`, `Interface_BuriedArea`,
`Hotspot_Contact_Fraction`, `Surface_Hydrophobicity`, `Binder_Length`, `Binder_Net_Charge`,
`Binder_Free_Cysteines`, `rank`, `hash`, `trajectory`, `outcome`, `failed_filters`, plus
`cif_path`, `n_targets`, `binder_chain`, `binder_chain_src`, `_status` and `_path`.
`--metrics a,b,c` adds more; `--all-metrics` takes every column (~60).

### The binder's chain letter

`bindcraft2_binder_chain` is the chain **letter** the binder carries in the collected structure, and `bindcraft2_binder_chain_src` says where it was read from. Pass the letter to `atomium --chains-to-design`, `ssprofile --design-chains`, `chainsel --chains` and `cms --binder-chains`.

**The `binder_chain` *setting* is not that letter.** Its default is `null` → `"binder"`, which is the campaign's internal chain *name*. Upstream writes chains sorted on `(is_binder_chain(name), name)` (`protein.written_chains`), so **targets come first and the binder is last**: `B` behind a single-chain target, `C` behind a two-chain one, and further along when a receptor is fused from several chains. Reading the setting instead of the file is how a `--chains-to-design` ends up redesigning the target.

| `binder_chain_src` | Means |
| --- | --- |
| `stamp` | read from `_bindcraft.redesigned_residues` in the accepted mmCIF, whose spans upstream itself writes against the output letters. Exact, and correct for an oligomeric binder (several letters, e.g. `CD`). |
| `last_chain` | the last chain in the written structure. Right for one binder chain, an **under-count when `copies > 1`** — only the trajectories/refolded stages land here, since upstream stamps only accepted designs. |
| `none` (NA) | neither worked. Look at the file; do not guess. |

### A multi-target campaign: one row, one column set per target

A `--targets` campaign still collects **one row per design** — one binder sequence is one entity, however many targets it was scored against. What multiplies is the columns. Upstream collapses every per-target reading into one `;`-joined cell ordered by descending weight then name (`campaign_output.target_ordered_row`) and adds its own key columns; the collector splits each cell back out, **keyed on the row's own `targets` cell** — never by position, never by the alphabetical order of the files on disk.

| Column | Meaning |
| --- | --- |
| `bindcraft2_targets`, `bindcraft2_target_weights` | the `;` cells the split was keyed on — the audit trail |
| `bindcraft2_n_targets`, `bindcraft2_on_targets`, `bindcraft2_off_targets` | the count, and the names split into binders and counter-selections |
| `bindcraft2_<metric>_<target>` | an **on-target** reading, e.g. `bindcraft2_i_pTM_hEGFR` |
| `bindcraft2_<metric>_off_<target>` | an **off-target** reading, e.g. `bindcraft2_i_pTM_off_hERBB2` |
| `bindcraft2_path_<target>` / `bindcraft2_path_off_<target>` | that target's complex, as PDB |
| `bindcraft2_path` | the **highest-weight on-target** complex — so every downstream step that expects one structure column keeps working |
| `bindcraft2_cif_path` | the source mmCIF of `bindcraft2_path` |
| `bindcraft2_n_complexes` | how many of the targets actually had a structure |

The bare `bindcraft2_i_pTM` carries **no value** on a multi-target table — there is no single value, and inventing one would mean averaging a binding reading with a counter-selection reading. Shared (non-per-target) metrics such as `Binder_Length`, `sequence`, `hash` and `failed_filters` stay single columns.

> **Corrected 2026-09-30, measured on a real 2-target campaign.** This section previously said the bare columns "do not exist". **They DO exist — they are present and entirely `NaN`.** Observed on a trajectories-stage collect: `bindcraft2_<label>_i_pTM`, `_i_pAE`, `_pLDDT`, `_pTM`, `_Interface_Residues`, `_Interface_BuriedArea`, `_Hotspot_Contact_Fraction`, `_Surface_Hydrophobicity`, `_Binder_Net_Charge` and `_Binder_Free_Cysteines` were all present and all NaN, while the real values sat in the `_<target>` columns. Nothing is corrupted and no value is invented, but **do not write a filter that relies on a bare column being absent** — test for the per-target column you actually want, or for `notna()`. `_off_` columns are genuinely absent when there are no off-targets (`off_targets` is NaN).

**Off-targets are kept separate on purpose.** Upstream excludes detargets from its own ranking (`campaign_output.on_target_mean`), because averaging a detarget `i_pTM` into the score rewards binding the thing you are avoiding. The `off_` infix carries that separation into the table; `run bindcraft2 --targets` refuses a target actually *named* `off_*` so the infix stays unambiguous.

**Filters you can now write.** A cross-reactive binder that avoids the paralogue.

> **`-f` takes a FILE PATH, not an expression.** Verified 2026-09-30 in `core/base_run.py`: the flag is *"Path to a Python module defining an `apply_filter(df) -> df` function"*, loaded with `importlib.util.spec_from_file_location`. Passing a predicate string fails with `ImportError: Could not load module from <your expression>`. An earlier revision of this skill showed the expression form; it never worked. The predicates below are correct — they just belong *inside* a module.

```python
# filters/bc2_crossreactive.py   ->   sapia run ... -f filters/bc2_crossreactive.py
def apply_filter(df):
    return df.query(
        'bindcraft2_i_pTM_hEGFR > 0.7 and bindcraft2_i_pTM_mEGFR > 0.7 '
        'and bindcraft2_i_pTM_off_hERBB2 < 0.4'
    )
```

Other predicates worth keeping, same wrapping:

```python
'bindcraft2_status == "OK" and bindcraft2_n_complexes == bindcraft2_n_targets'
'abs(bindcraft2_i_pTM_hEGFR - bindcraft2_i_pTM_mEGFR) < 0.1'     # balanced, not lopsided
```

The path resolves against the **submit cwd** (`/runs` on Modal), so the module has to live where `sapia` runs, not only on your laptop.

### Error contracts (multi-target)

These are the point of the split, not decoration:

- **A `;` value count that disagrees with `targets` is an `error` status**, not a best-effort parse — the design's per-target columns are not written at all, and `bindcraft2_path` is empty. `target_ordered_row` always joins over *every* target name, so this should never fire; if it does, something is genuinely wrong and you want to see it.
- **`targets` and `target_weights` of different lengths** is the same error: without the weights an on-target cannot be told from a detarget.
- **A missing per-target complex is `NA`** in `bindcraft2_path_<target>`. Never a substituted sibling file. (Filter on `n_complexes == n_targets` if you need all of them.)
- **An empty position in a `;` cell is `NA`, not `0`** — upstream's own docs say do not read a blank as zero, and a `0` `i_pTM` reads as "does not bind" rather than "was not measured".
- **If the highest-weight on-target has no complex, the row is an `error`** — a detarget's complex is never promoted to `bindcraft2_path`.

### Traps found while building this

- **The filename shape changes when a second target appears.** `accepted_state_suffixes` returns an *empty* suffix below two prepared states, so a single-target campaign writes `<design>.cif` and a multi-target one writes `<design>_<target>.cif`. There is no suffix to strip and no glob that is safe: the collector matches exact names only.
- **File order is alphabetical; cell order is by descending weight.** "The first file" and "the first value in a cell" are different targets. Anything keyed on position is wrong.
- **`Timing` uses `;` too.** `timing_stamp` writes `worker=0;start=…;reprediction=…`, so on a three-target campaign a naive split produces three plausible columns of nonsense. The collector keeps an explicit never-split set (`Timing`, `failed_filters`, `hash`, `terminated`, `Binder_Sequence`, `Interface_*_Residues`, `settings_*`, …).
- **`failed_filters` is comma-joined, and its *contents* gain the target suffix** (`i_pTM.mEGFR,pLDDT.hEGFR`). It stays one cell; read it to see which target rejected a candidate.
- **A FASTA/IDR target can be cropped into several states** named `<name>_epitope_<i>`, whose written filename suffix is `<source>_<start>-<end>` instead — so the `targets` names and the filenames stop agreeing and the per-target paths come back NA. Not exercised here; if you crop a disordered target, check the paths.
- **A design that only predicted one state** in a multi-target campaign gets an *unsuffixed* file while its CSV row still names every target. The collector reports that as an error rather than attaching the unsuffixed file to an arbitrary target.

Two things about `bindcraft2_path`:

- It is the **complex** — binder *and* target — converted to PDB. `usalign` and `cms`
  want the binder alone, so `chainsel` it out first for those. Sequence redesign
  (`atomium`, `proteinmpnn`) is the exception: keep the complex and name the binder chain,
  so the target context is there.
- A design whose mmCIF will not parse still collects, with its metrics and a
  `_status` of `error: unreadable mmCIF`, pointing at the `.cif`. One bad file does not
  abort the collect.
- On a multi-target campaign it is the **highest-weight on-target** complex; the others
  are in `bindcraft2_path_<target>`. Say which target you mean before comparing poses.

`bindcraft2_sequence` is what a predictor consumes — not Boltz's default column, so a
Boltz run after bindcraft2 needs `-i bindcraft2_sequence`.

## Its own scores are not independent validation

`i_pDAE`, `i_pTM` and pLDDT here come from the **same AF2 that designed the binder**. They
are the objective the campaign optimised against, so a high value is partly a statement
that the optimiser succeeded, not that the binder is real. BindCraft2's internal filters
are a floor, not evidence.

The honest follow-ups are the ones AF2 did not see:

- **A different predictor** on the complex — `boltz` or `alphafold3` with
  `-i bindcraft2_sequence`. Agreement between two architectures means something; AF2
  agreeing with itself does not.
- **`usalign`** back to the predicted pose, to confirm the independent prediction puts the
  binder in the same place rather than merely folding it.
- **`pyrosetta`** for an energy that is not a confidence score at all.

Read `bindcraft2_rank` as position in *that* ranking and nothing more.
