---
name: bindcraft2
description: How to run the custom bindcraft2 tool on Modal — BindCraft2 binder-design campaigns, where one task is a whole campaign rather than one design. Covers the campaign cost model (--num-designs vs --max-trajectories), root vs child runs, multi-target and detargeting campaigns (--target/--extra-target/--shipped-target, weights and objectives) and how their packed metric cells are split into filterable per-target columns, hotspots with {expr}, --trajectory-only for handing backbones to atomium/proteinmpnn instead of BC2's own MPNN, --reuse-campaigns for collecting one campaign into two tables, --fetch-weights-only for warming the AlphaFold Volume, the table/dir labels that fork needs and the silent table collision when you omit them, the three collect stages, and why its own confidence scores are not independent validation. Load before composing a bindcraft2 run or reading its columns.
---

# bindcraft2

**Custom tool** (lives in `tools/bindcraft2/`, not bundled with prosapia).
**`action: create`** — mints a child table, one row per binder design a campaign produced.

Upstream: [PacesaLab/BindCraft2](https://github.com/PacesaLab/BindCraft2), pinned to tag
`v1.0.3` in `modal_image.py`. A reference checkout lives at `../BindCraft2` with a CPU venv
at `.venv-cpu`; its `docs/source/reference.md` is authoritative for settings,
`docs/source/outputs.md` for columns. `tests/test_bindcraft2.py` validates every settings
file this tool generates against **BindCraft2's own validator** in that venv — run it after
touching the tool.

> **Not yet run end to end on real hardware.** The manifest builder, the collector and the
> task script are covered by tests (including upstream validation and a real bash run), but
> no campaign has been through Modal yet. The first real submit should be one small
> campaign, not a fan-out — see *First run*.

## The thing to understand first: a task is a campaign

Every other design tool here maps one design to one unit of work. BindCraft2 does not.
One `bindcraft design` process takes a target (or several) and loops: hallucinate a
backbone with AF2, redesign its sequence with ProteinMPNN, refold it, filter it, keep or
discard — until `--num-designs` have been **accepted** or `--max-trajectories` attempts are
spent.

Consequences that shape how you drive it:

- **One manifest row = one campaign.** `-C/--max-concurrent` throttles campaigns, not designs.
- **The row count is unknown until collect.** A campaign can accept fewer designs than
  you asked for, or none at all. That is a normal outcome, not a failure.
- **`--num-designs` is a stopping condition, not a batch size.** The real cost knob is
  `--max-trajectories`: a campaign runs until one of the two fires.
- **It already contains the rfd3 → mpnn → predict → filter chain.** Do not compose those
  tools around it; compose *independent* checks after it (see the last section).

## Two shapes of run

- **Root run (no `-t`)** — starts a fresh lineage in `table0`. The target comes from
  `--target-pdb <file>` (group `<stem>_bc2`), `--shipped-target <name>` (group
  `<names>_bc2`), or one or more `--target <spec>` (group named after the binding targets).
- **Child run (`-t <table>`)** — one campaign per ready row, the **primary** target from
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

## Multi-target and detargeting

A campaign may carry **several targets**: orthologs to bind cross-reactively, and
off-targets to steer away from. This is first-class upstream — three of its shipped
examples are multi-target — and it is what the `multitarget_*` and `max_detarget_iptm_*`
settings exist for.

Each target is one `--target` (root) or `--extra-target` (added beside the table's own
target), written as **`;`-separated** `key=value` fields in BindCraft2's own `targets[]`
vocabulary. The separator is `;` and not `,` because hotspot lists already use commas.

| Field | Meaning |
| --- | --- |
| `name` | Required, distinct. Identifies the target in the tables **and in every per-target structure filename**. Defaults to the file stem. |
| `path` (= `target_path`) | PDB, mmCIF or FASTA. A FASTA target is treated as disordered and cropped (`crop_fasta_sequence`). |
| `chains` | Input chain selection, `A` or `A,B`. |
| `hotspots` / `coldspots` | Residues to contact / avoid, `A54,A56,B12-16`. |
| `weight` | Relative importance, default 1. **Negative also selects detargeting.** |
| `objective` | `target` (default) or `detarget`. |

```bash
# cross-reactive against two orthologs, steering away from a third protein
sapia run bindcraft2 <run_dir> \
    --target 'name=hPDL1;path=t/hPDL1.pdb;chains=A;hotspots=A54,A56' \
    --target 'name=mPDL1;path=t/mPDL1.pdb;chains=A;hotspots=A36,A38' \
    --target 'name=hPD1;path=t/hPD1.pdb;objective=detarget;weight=-0.5' \
    --num-designs 10

# the table's own target, plus an off-target to avoid
sapia run bindcraft2 <run_dir> -t table0 -i pdb_path \
    --extra-target 'name=hPD1;path=t/hPD1.pdb;weight=-0.5'

# two shipped targets accumulate into one cross-reactive campaign
sapia run bindcraft2 <run_dir> --shipped-target hPDL1 --shipped-target mPDL1

# or hand the whole list over in a file
sapia run bindcraft2 <run_dir> --targets-file targets.yaml
```

`--chains/--hotspots/--coldspots/--target-weight/--target-objective` apply to the
**primary** target only (the table row, or `--target-pdb`). Put per-target values inside
each `--target` spec.

### Detargeting needs a filter, not just a weight

A repulsion weight steers away from an off-target without ever establishing that the
design got away, **and interface confidence does not say it either** — upstream's own
comment records a peptide reading `0.27` i_pTM with its whole face on the off-target.
BindCraft2 therefore defaults `max_detarget_interface_residues_final` to 3 whenever an
off-target is present. Read `bindcraft2_Interface_Residues_detarget`, not just the
confidence. `max_detarget_iptm_<stage>` (unset by default) adds a per-stage ceiling.

### The one combination that is refused

Mixing `--shipped-target` with a described `--target` is **refused at submit time**.
Upstream accepts it and then silently drops the shipped ones: preset target entries are
layered *under* the campaign request, and a list replaces rather than merges
(verified against v1.0.3 — `{'target': ['hPDL1','mPDL1'], 'targets': [...]}` resolves to
the explicit list alone). A campaign that quietly designs against one target when told
three is worse than a refused submit. To combine them, write the shipped target out as its
own `--target`; its structure and hotspots are in `../BindCraft2/settings/target/<name>.json`.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--targets FILE` | none | YAML/JSON list of targets **one** binder is optimised against jointly. Negative `weight` = counter-select. Excludes `--target-pdb`/`--shipped-target` and the top-level `--chains`/`--hotspots`/`--coldspots`. See above. |
| `--trajectory-only` | off | Stop at **backbones** — no MPNN, no refold, no filter, nothing accepted. The handoff to this workspace's own designers; see below. |
| `--fetch-weights-only` | off | Run **no campaign**: one task that downloads and verifies the ~5.3 GB AlphaFold parameters into the cache Volume. Asks for no GPU. See *The AlphaFold parameters*. |
| `--reuse-campaigns` | none | `TABLE[:LABEL]`. Submit nothing; reserve a second table over an earlier run's campaigns, to collect a second stage of the same GPU hours. See below. |
| `--max-trajectories` | BC2's own | **The cost knob.** Attempts spent before giving up. Set it or a hard target can burn the whole timeout. With `--trajectory-only` it is the *only* budget (BC2 defaults it to 100). |
| `--num-designs` | BC2's own (1) | Accepted designs to stop at (`number_of_final_designs`). A target, not a guarantee. |
| `--target` / `--extra-target` | none | A described target; repeatable. See above. |
| `--shipped-target` | none | Repeatable. `hPDL1`, `hPD1`, `mPDL1`, `hIL2R`, `hIL7RA`, `dynorphin_a`; `bindcraft design --list-targets` is authoritative. |
| `--targets-file` | none | YAML/JSON holding the whole `targets` list. |
| `--hotspots` / `--coldspots` | none | Primary target's epitope, BC2 syntax (`A54,A56,A66-70`). Takes `{expr}`. Omit to let BC2 pick. |
| `--chains` | all | Primary target's chains (`A`, `A,B`). |
| `--binder-lengths` | modality's | `80`, `60-100` (range) or `60,80,100` (discrete set). Takes `{expr}`. |
| `--modality` | `binder` | `binder`, `VHH`, `peptide`, `cyclic_peptide`, `ARP`, `scFv`, `Fab`, `large_binder`, `homo_oligomer`, `multidomain`, `induced_fit`, `fold_switch`. Comma-separated to combine. |
| `--property` | none | Repeatable preset: `humanize`, `protease_stable`, `disulfide_staple`, `forced_targeting`, `initial_guess`, `mixed_topology`, `termini_together`, `termini_accessible`, `bigbang`. Validated at submit time; an unlisted one is reachable as `--set <name>=true`. |
| `--core` | none | Core profile under every preset; `benchmark` for a reproducible run. Pair with `--campaign-seed`. |
| `--metadata` | none | JSON of descriptive fields (author, project, note), recorded as `meta_*` columns. Provenance, not settings. |
| `--save-monomers` | off | Also keep the binder re-predicted **alone** (`save_binder_monomers`) → `bindcraft2_monomer_path`. Saves a later `chainsel` for self-consistency. |
| `--design-workers` / `--workers-per-gpu` | 1 / auto | Fan one campaign's trajectories across several GPUs. Only helps with `-g N` to match. |
| `--extra-settings` | none | YAML/JSON merged into every campaign's settings. Values take `{expr}`. |
| `--set KEY=VALUE` | none | Verbatim `bindcraft design --set`. **Spaces and JSON are fine** — each token is carried as its own argv entry. |
| `--no-resume` | off | See below. |

`bindcraft design --list-settings` names every setting `--extra-settings` and `--set`
accept. This tool writes the settings file itself; a key set by both a flag and
`--extra-settings` is refused at submit time rather than silently resolved — **except
`targets`**, which `--extra-settings` may own outright as long as no target flag is in play
and the run has no `-t`.

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
# 1. one FULL campaign, reserving table0. Collect its backbones.
sapia run bindcraft2 R -l traj --hotspots 'A54,A56' --num-designs 10
sapia collect bindcraft2 R -t table0 -l traj --stage trajectories

# 2. the SAME campaign's accepted designs, into table0_bc2. No GPU.
sapia run bindcraft2 R --table-label bc2 --reuse-campaigns table1:traj
sapia collect bindcraft2 R -t table1_bc2 --stage ranked

# 3. the atomium branch off the backbones, reserving table2
sapia run atomium R -t table0 -i bindcraft2_traj_path --chains-to-design <binder chain>
```

Pass the **same `-t` the original run used** (here `table0`). The collector maps each
ready design to a campaign folder by name, and those folders are named after the
*original* run's design groups — reusing with `-t table1` would hunt for campaigns named
after backbones and find none. A root campaign needs no `-t` at all; the reuse run reads
the group names off disk.

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
these are the only download. The same Volume also collects jax's compiled graphs
(`compile_cache/<card>`), so later campaigns on the same GPU model skip the compile.

**Warm it before fanning out.** Submit N campaigns cold and N containers each pull the same
5.3 GB.

```bash
sapia run bindcraft2 <run_dir> --fetch-weights-only --table-label weights
```

That submits one CPU-only task running `bindcraft fetch-weights`, which downloads **and
verifies** both checkpoint sets and exits non-zero if any is missing or unfinished — which
a dummy campaign does not. It collects nothing by design (the collector says so and
stops); give it its own `--table-label` so the empty table it reserves does not shadow a
real one. `BINDCRAFT_AF2_PARAMS` (parameters only) or `BINDCRAFT_WEIGHTS` (whole cache
root) in `.env` override the location.

Per the workspace convention the Volume is named without the `sapia-` prefix so it can be
shared. **Never delete it to tidy up** — it re-downloads, slowly.

## First run

The image builds inside the submit call: Ubuntu + jax cuda13 wheels + an editable install
of the repo. Allow **15+ minutes** before assuming a first submit is stuck; later runs are
seconds. Then:

1. `--fetch-weights-only` to warm the Volume and prove the image.
2. One target, `--max-trajectories 10` or so, and watch it finish.
3. Check `campaigns/<name>/1_Trajectories/!_Trajectories.csv` appears. If jax fell back to
   the CPU the campaign still *runs*, ~100× slower, and says so only in a buried warning —
   that is the failure mode the image's `ldconfig` step exists to prevent.
4. Only then fan out.

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
  needs a `--reuse-campaigns` run, not just a second `-l`.
- **`ranked`** (auto otherwise) — `3_Ranked/!_Ranked.csv`, only what the campaign accepted.

All three carry the same metric battery: BC2 scores a trajectory on the same filters it
runs at refold, so the columns are comparable across stages.

If `archive_trajectories` was on, the per-design folders are zipped — the collector reads
the structures straight out of the zip into `<out_dir>/.unarchived/` and leaves the
campaign folder untouched. No manual `bindcraft unarchive` needed.

Columns (leaf-prefixed `bindcraft2_`): `sequence`, `i_pDAE`, `i_pTM`, `i_pAE`, `pLDDT`,
`pTM`, `Unbound_Binder_pLDDT`, `Target_pLDDT`, `Interface_Residues`,
`Interface_BuriedArea`, `Hotspot_Contact_Fraction`, `Coldspot_Contact_Fraction`,
`Interface_Binder_Residues`, `Interface_Target_Residues`, `i_pTM_detarget`,
`i_pAE_detarget`, `Interface_Residues_detarget`, `Backbone_Clashes`, `All_Atom_Clashes`,
`Binder_Chain_Breaks`, `Surface_Hydrophobicity`, `Binder_Length`, `Binder_Net_Charge`,
`Binder_Free_Cysteines`, `targets`, `target_weights`, `rank`, `hash`, `trajectory`,
`outcome`, `failed_filters`, plus `cif_path`, `_status` and `_path`. Any `meta_*` (from
`--metadata`) and `settings_*` (which presets actually ran) column is **always** collected.
`--metrics a,b,c` adds more; `--all-metrics` takes every column (~60).

Three things about `bindcraft2_path`:

- It is the **complex** — binder *and* target — converted to PDB. `usalign` and `cms`
  want the binder alone, so `chainsel` it out first for those (or run with
  `--save-monomers` and use `bindcraft2_monomer_path`). Sequence redesign (`atomium`,
  `proteinmpnn`) is the exception: keep the complex and name the binder chain.
- In a multi-target campaign it is the **highest-weight binding target's** complex; every
  other target's is beside it as `bindcraft2_path__<target>`.
- A design whose mmCIF will not parse still collects, with its metrics and a
  `_status` of `error: unreadable mmCIF`, pointing at the `.cif`. One bad file does not
  abort the collect. A structure matched by glob rather than by target name collects with
  `_status` = `OK: structure matched by glob, not by target name` — treat that row's
  target assignment as unverified.

`bindcraft2_sequence` is what a predictor consumes — not Boltz's default column, so a
Boltz run after bindcraft2 needs `-i bindcraft2_sequence`.

### Reading a multi-target table

With **two or more** targets BindCraft2 does not write one column per target. It packs
each metric into a single semicolon-separated cell, in the order of the row's own
`targets` column — which is sorted by `(-weight, name)`, so off-targets come last. (With
fewer than two targets the row is left completely untouched, so single-target tables look
exactly as they always did.)

A packed cell is a **string**: `-f` cannot threshold it and a sort orders it as text. So
the collector keeps the packed cell verbatim *and* splits it:

```
bindcraft2_targets             hPDL1;mPDL1;hPD1
bindcraft2_target_weights      1;1;-0.5
bindcraft2_i_pTM               0.82;0.79;0.21     <- verbatim, as BC2 wrote it
bindcraft2_i_pTM__hPDL1        0.82               <- filterable
bindcraft2_i_pTM__mPDL1        0.79
bindcraft2_i_pTM__hPD1         0.21
bindcraft2_i_pTM__mean         0.805              <- binding targets only
bindcraft2_i_pTM__worst        0.79
bindcraft2_i_pTM__best         0.82
bindcraft2_i_pTM__spread       0.03
bindcraft2_i_pTM__selectivity  0.58               <- weakest binder vs best off-target
```

- `__mean`, `__worst`, `__best`, `__spread` cover the **binding targets only**, mirroring
  upstream's `on_target_mean`. Averaging an off-target in would reward binding the thing
  the campaign is avoiding.
- `__selectivity` is the one column that deliberately compares the two groups: the weakest
  binding target against the strongest off-target, in the metric's own direction, so a
  **positive margin always favours the intended distinction**. It is a computational score
  margin, not an affinity ratio.
- Summaries appear only for metrics whose direction is known (`METRIC_DIRECTION` in
  `collect_bindcraft2.py`). An unknown metric still gets its `__<target>` split.
- A cell is split only when it holds exactly as many positions as the row has targets —
  the same guard upstream applies, and what stops a residue list containing a semicolon
  from being mapped onto the wrong targets.
- `--no-split-targets` turns the extra columns off and leaves the packed cells alone.

**Do not average a binding and a detarget state together**, and do not read a blank
position as zero — upstream says so explicitly, and the split preserves blanks as empty.

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
- **`cms`** for an interface measure that is geometric rather than predicted.

Read `bindcraft2_rank` as position in *that* ranking and nothing more.
