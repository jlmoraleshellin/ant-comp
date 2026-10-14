---
name: rfdiffusion3
description: >-
  Generates protein backbones with rfdiffusion3 (rfd3): de novo, motif-scaffolded, symmetric or binder-against-a-target. Load before composing any rfdiffusion3 run or writing a contig or hotspot spec. Traps: the chain-break literal is `/0` (the colon form silently designs against half the target); `--length min-max` draws one length per batch; rfd3 renumbers every output chain from 1, so downstream residue numbers inherit an offset.
---

# rfdiffusion3

Generates protein backbones. **`action: create`** — always mints a new table.

Full reference: `docs/tools/rfdiffusion3.md` in the prosapia repo, and
`sapia run rfdiffusion3 --help` (authoritative for flags).

IMPORTANT: prosapia's bundled rfdiffusion3 is a wrapper of the original. Check the upstream repo for
model behaviour: https://github.com/RosettaCommons/foundry/tree/production/models/rfd3/docs

## Premise, and what it does not establish

rfd3 produces **backbone coordinates**. That is the whole claim. It does not produce a sequence, and
nothing it writes says whether the backbone is designable, soluble, or placed somewhere a binder
could physically sit. Every judgement of quality happens downstream: sequence design, then a
predictor, then the fold and pose gates.

So treat a batch of rfd3 output as **candidates, not results**. The columns it collects are
provenance (which batch, which model, which length), not scores — there is nothing here to rank on.
`rfd3_ca_rmsd_to_input` is the only number with any evaluative content and it is **NaN for every de
novo run**, since there is no input to deviate from.

## Three facts that change what your run means

Read these before writing the command. Each produces a successful-looking run that did something
other than what you intended.

**1. `/0` is the only chain-break literal, and the colon form fails silently.**
`elif part == "/0"` is the only break the parser recognises. The colon form
`A18-20:B18-20` does **not** error — it matches the prefix only and **drops chain B**, so you design
against half the target and every interface number afterwards describes a different problem.

**2. `--length min-max` does not vary length within a batch.** rfd3 draws **one** length per batch,
so `--num-designs 5` gives five backbones of the *same* length. If you believe you ran a length
sweep, you ran five replicates. For a spread use several batches
(`--num-designs 1 --set n_batches=5`) or several design keys, at the cost of loading the model more
often. The drawn value lands in `extra["sampled_contig"]`.

**3. rfd3 renumbers every output chain from 1.** A target numbered 18–155 comes back as 1–138. Every
downstream step that names target residues — hotspot recall, epitope contacts, a residue mask, an
RMSD on selected residues — inherits that offset, and nothing in the table announces it. Hotspots in
the *input* spec use the input file's own numbering; everything read off the *output* does not.

## Two shapes of run

- **Root run (no `-t`)** — starts a fresh lineage in `table0`. One design group, named `denovo_diff`
  (pure de novo) or `<pdb-stem>_diff` (with `--input-pdb`). Rows collect as `<group>_0`,
  `<group>_1`, …
- **Child run (`-t <table>`)** — one design group per ready row, inputs from `-i/--input-column`.
  Mints a child table.

`{expr}` placeholders in contigs/length resolve up the lineage, so they need `-t`. **A root run must
use literal values.**

## Verified de-novo invocation

```bash
sapia run rfdiffusion3 <run_dir> --length 80-120 --num-designs 5
```

No `-t`, no contig, no input. Produced 5 backbones in `table0` in ~3 min on an A100 (≈53 s of that
was inference; the rest was container start and model load).

## Flags

| Flag | Default | Note |
| --- | --- | --- |
| `--length` | none | `N` or `min-max`. **One length drawn per batch** — see fact 2. One subunit's length when symmetric. |
| `--contigs` | none | rfd3 contig syntax: `A40-60` motif, bare `30` designed, `/0` chain break. |
| `--num-designs` | 1 | Designs per input key (`diffusion_batch_size`). Model loads once for all. |
| `--shard-size` | 10 | Design keys packed per task; they run **in series** in one process. |
| `--symmetry` | none | `C5`, `D4`, or `auto`. Write **one asymmetric unit's** contig — the sampler replicates it. There is no `--replicate`. |
| `--input-pdb` | none | Root-only, for diffusing a structure not yet in a table. |
| `--num-timesteps` | 200 | |
| `--step-scale` | 1.5 | Higher = less diverse, more designable. |
| `--partial-t` | none | Partial diffusion noise, ~5–15 Å. |
| `--extra-spec` | none | YAML/JSON of extra rfd3 spec fields merged into every design. |

Default Modal resources: **A10**, 8 CPU, 32 GiB, 4 h timeout. Override with `--modal-gpu`.

## Binder design: hotspots and the two-target contig

There is no binder flag. A binder job is a **root run** that holds the target chains as motif and
appends one designed chain, with hotspots supplied through `--extra-spec`.

```bash
sapia run rfdiffusion3 <run_dir> \
    --input-pdb target/7ojg_AB.pdb \
    --contigs 'A18-155,/0,B18-155,/0,70-100' \
    --extra-spec target/hotspots.yaml \
    --num-designs 2 --set n_batches=4 --modal-gpu A100
```

Verified: 8 backbones, 3 chains each (`A` 138 / `B` 138 / **`C` = the binder**), ~4 min on an A100
for a 276-residue motif plus a ~90-residue binder.

### `select_hotspots` — the syntax

Not documented anywhere in prosapia (`grep hotspot docs/` returns nothing), and **classic
RFdiffusion's `ppi.hotspot_res=[A30,A33]` list form does not carry over.** In rfd3 the field is typed
`Optional[InputSelection]` on `DesignInputSpecification`, and `InputSelection.from_any` accepts only
**a contig-style string, a bool, or a dict** — a list raises
`ValueError: Cannot convert <class 'list'> to InputSelection`.

```yaml
# hotspots.yaml — a whole-residue selection (ALL atoms of each residue)
select_hotspots: "A47,A123,B42,B155"     # ranges work too: "A47-49"
infer_ori_strategy: hotspots             # places the origin token out along the
                                         # outward normal from the hotspot COM
```

**Correction, verified from source 2026-10-09 (`rc_foundry/rfd3/utils/inference.py`,
`infer_ori_from_hotspots`): the push-out distance is a hardcoded `10.0 Å`, not 12 Å.** 12 Å is a
*different* constant — a neighbor-search radius used only to find the "core" motif atoms near the
hotspots, which establishes the core→hotspot outward direction the origin is then pushed along.
Neither constant is a config field or a `--set`-able override; both are Python literals (the source
itself flags this as a wart: `"Both this and the distance cutoff should definitely be configs with
defaults"`). **There is a second strategy, `infer_ori_strategy: com`**
(`INFER_ORI_STRATEGIES_TO_FUNCTIONS = {"hotspots": ..., "com": infer_ori_from_com}`). **It is not a
binder-aware centroid — read this before using it.** `infer_ori_from_com` is a plain unweighted
`np.mean` over every finite-coordinate atom in the whole assembled complex, no chain/motif/hotspot
selection of any kind (contrast `infer_ori_from_hotspots`, which explicitly masks on
`is_atom_level_hotspot`). It runs *before* a design's free/designed atoms are zeroed for diffusion,
so on a typical binder contig — a large fixed target plus a small motif plus a free/designed chain —
the mean is **dominated by the fixed target's own centroid**, diluted only slightly by the motif.
It is not anchored on any interface or epitope and is independent of whether `select_hotspots` was
ever set. Exposes zero parameters (no distance, no selection, no weighting) — if anything more
rigid than the hotspots strategy, just with nothing to tune. **Verify this is actually what you
want before using it on a target that heavily outweighs the motif in residue count** — it is easy
to read the name "com" and assume it means "center of the binder." **To bypass both strategies and their hardcoded distances entirely**, set `ori_token:
[x, y, z]` directly (sibling field, same `DesignInputSpecification`) — it takes precedence, and
needs **no compensation for any hidden prior centering**: `set_com` (the function this resolves
into) performs the *only* raw-frame translation applied to the loaded input coordinates for a
non-symmetric run, so a coordinate computed from the raw input PDB can be passed as-is. One
caveat, unresolved by static reading: whether the sampler's own per-timestep stochastic
re-centering (`centre_random_augment_around_motif`, `rfd3/model/inference_sampler.py`) preserves
that placement exactly through every denoising step, or only in distribution.

```yaml
# atom-level form: the dict picks which atoms carry the annotation
select_hotspots:
  A47: TIP        # sidechain tip atoms
  B42: BKBN       # backbone
  B155: [CA, CB]  # explicit atom names
```

That string-vs-dict choice is exactly what the field's "atom-level or token-level" docstring refers
to.

**Specify few hotspots.** rfd3 was trained with hotspots present in 75% of PPI examples, showing only
a random subset of up to **20%** of the true hotspot atoms (ground truth = target atoms within 4.5 Å
of the binder). A dense patch is off-distribution. It is real conditioning, not a no-op: the shipped
checkpoint carries trained `token_initializer.…is_atom_level_hotspot.weight` tensors.

Constraints worth knowing:

- Hotspots need `input` set in the same spec, and must lie **inside the contig's motif ranges** —
  annotations on residues outside the contig never enter the built structure.
- Chain/residue ids are the **input file's own numbering**, not renumbered. (The *output* is
  renumbered — fact 3.)
- prosapia resolves `{expr}` in every `--extra-spec` string, so avoid literal braces.
- Setting `contig`/`length`/`input`/`symmetry`/`partial_t` in both a flag and `--extra-spec` is a
  `SpecConfigError`. The task script prevalidates, so a bad spec fails fast and cheaply — **let it,
  rather than building a probe container to check.**

### Low-temperature sampling — verified from source 2026-10-09

rc-foundry's own binder-design documentation recommends `is_non_loopy: true` plus
`inference_sampler.step_scale=3` and `inference_sampler.gamma_0=0.2` together, for "greatly
increased PPI designability" over the defaults. Confirmed mechanics, since two of these three are
easy to get wrong:

- **`is_non_loopy`** is a real top-level field on `DesignInputSpecification`
  (`Optional[bool]`) — goes in `--extra-spec`, same as `select_hotspots`.
- **`--step-scale`** (the CLI flag already in the Flags table, default 1.5) **is**
  `inference_sampler.step_scale` — confirmed by tracing `run_rfdiffusion3.py`'s override-building
  straight into the Hydra `compose(overrides=...)` call `rfd3 design` makes. No separate mechanism
  needed; `--step-scale 3` is the whole fix for that one.
- **`gamma_0`** has **no dedicated CLI flag**. The only way to set it is the generic
  `--set inference_sampler.gamma_0=0.2` — confirmed to resolve as a genuine dotted Hydra override,
  the same path `step_scale` goes through internally. Default confirmed `0.6` (matches
  `rc_foundry/rfd3/configs/inference_engine/rfdiffusion3.yaml`).
- So one command needs: `--step-scale 3 --set inference_sampler.gamma_0=0.2` plus
  `is_non_loopy: true` in `--extra-spec`. The CLI flag and `--set` compose on the same invocation
  without conflict (different keys).

### Burial conditioning — `select_exposed` / `select_buried` / `select_partially_buried`

Three sibling fields to `select_hotspots`, same `Optional[InputSelection]` type, same contig-style
string syntax, same **input-file numbering** (not renumbered output), same requirement that `input`
be set first. Verified exact field names from source — don't assume the names without checking;
this codebase moves (`dev` keeps moving — see the top-level project note).

```yaml
select_exposed: "D30,D31,D33"     # RASA-exposed conditioning
select_buried: "D29,D51"          # RASA-buried conditioning
select_partially_buried: "D60"    # the middle bin, if wanted
```

**These three are mutually exclusive, and it is enforced as a hard error, not a soft preference.**
A residue appearing in more than one of the three raises `ValueError` at spec-build time
(`assert_exclusivity`, the same per-row build path any other bad field goes through) — no silent
"last one wins," no warning-only. Catch an overlapping residue list before submitting; it costs
nothing to check and the error arrives at build time regardless, so a bad list fails the whole
design, not just one annotation.

### Contig traps for a two-chain target

- **Use `/0` for the chain break** — see fact 1.
- **Do not pass `--length` with a binder contig.** It is the *total*, and
  `length_min -= num_motif_residues`, so `--length 400` against a 276-residue motif demands a
  124-residue binder and fails validation. Let the contig's `70-100` govern.
- **The designed length is drawn once per spec build** — see fact 2.

## Gotchas

- **Checkpoints live on a Volume** (`SAPIA_MODAL_VOLUME_RFD3_CKPT`, mounted at `/checkpoints`) and are
  **not installed automatically**. If a run fails for a missing checkpoint, the Volume is empty — say
  so rather than retrying. Populating it is a one-off
  `foundry install rfd3 --checkpoint-dir /checkpoints` from a container of the tool's image
  (~2.5 GiB, ~3 min).
- **This repo runs a FORK** (`tools/rfdiffusion3/`, same `name`, so it shadows the built-in and
  every column is unchanged; since 2026-10-07). Two deviations, both in its module docstrings:
  (1) the Modal image pins `rc-foundry[all]==0.2.0` + `atomworks[ml]==2.2.1` — unpinned, pip
  takes atomworks 3.0.0 and every task exits 1 at import with
  `ImportError: cannot import name 'get_inferred_polymer_bonds'`; bump the two together, never
  one alone. (2) `$RFD3_CKPT` is ignored when the run's executor is modal (`-e`, which the workstation
  defaults to `modal`) — the shared `.env` holds a
  vib cluster path there, whose signature on Modal is `ckpt_path=/data/groups/csb/...` in the
  task's overrides line. On Modal the checkpoint is auto-discovered from `/checkpoints`; to force
  one, pass `--ckpt-path /checkpoints/<file>` explicitly. vib behaviour is unchanged.
- Expect harmless stderr noise: an `atomworks` warning about an unset env var, and a log line showing
  Modal's internal `/__modal/volumes/...` path. Files still land under `/runs`.

## What it collects

Columns (leaf-prefixed `rfdiffusion3_`):

| Column | Meaning | How to read it |
| --- | --- | --- |
| `_status` | `OK` is the only proof the design succeeded | check before anything else |
| `_path` | the **PDB** collect converted from rfd3's `.cif.gz`, under `<run_dir>/.cif_to_pdb/` | this is the column to feed downstream |
| `_iteration`, `_rfd3_batch`, `_rfd3_model` | provenance: which draw, which batch, which checkpoint | use to tell replicates apart from genuine variants |
| `_rfd3_ca_rmsd_to_input` | CA RMSD to the input structure | **NaN for de novo.** Only meaningful for motif scaffolding or partial diffusion |
| per-chain length | the drawn length(s) | the actual length, which `--length min-max` only bounded |

There are **no quality columns here by design.** Nothing in this table says a backbone is good.

## Filters you can write against it

```python
# successful backbones in a length band, for a sequence-design step
def apply_filter(df):
    ok = df["rfdiffusion3_status"] == "OK"
    return df[ok & df["rfdiffusion3_length"].between(80, 110)]

# one backbone per batch, to buy a cheap end-to-end probe before committing
def apply_filter(df):
    ok = df[df["rfdiffusion3_status"] == "OK"]
    return ok.groupby("rfdiffusion3_rfd3_batch", as_index=False).head(1)

# motif scaffolding only: backbones that actually kept the motif
def apply_filter(df):
    return df[df["rfdiffusion3_rfd3_ca_rmsd_to_input"] < 1.5]
```

Confirm the exact length column name from the table header before writing the first filter — it is
per-chain and the suffix depends on the run shape.

## What it is blind to

- **Designability.** Nothing here says a sequence exists that folds to this backbone. That is the
  MPNN-then-predictor round trip, and the fold gate in `usalign --mm 0`.
- **Composition and solubility.** *Measured elsewhere:* a 4.5 %-charged, 70 %-TSVIG β-sandwich folds
  at pLDDT 0.96. Check % charged, net charge, % aromatic, % β-branched directly on the designed
  sequence; the backbone cannot show it.
- **Whether the hotspots were used.** Hotspot conditioning is real, but no column reports which
  target residues the binder ends up touching — and the generated pose is not the predicted pose.
  Hotspot recall must be measured on the **prediction** (`epitope`, `ringfit`), never on the rfd3
  output.
- **Whether the pose is physically possible.** A binder placed against two trimmed protomers may
  occupy space a symmetry-related protomer fills in the real assembly. `ringfit` against the full
  assembly answers that; rfd3 will generate it happily.
- **Its own diversity.** Because one length is drawn per batch, a batch can look like a varied
  ensemble and be replicates. Use `_rfd3_batch` to tell which you have.
