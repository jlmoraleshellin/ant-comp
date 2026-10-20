---
name: rfdiffusion3
description: >-
  Which backbones can rfdiffusion3 (rfd3) generate: de novo, motif-scaffolded, symmetric or binder-against-a-target? Load before composing any rfdiffusion3 run or writing a contig or hotspot spec. Traps: the chain-break literal is `/0` (the colon form silently designs against half the target); `--length min-max` draws one length per batch; rfd3 renumbers every output chain from 1, so downstream residue numbers inherit an offset.
---

# rfdiffusion3

Generates protein backbones. **`action: create`**: always mints a new table. Full reference: `docs/tools/rfdiffusion3.md` in the prosapia repo, and `sapia run rfdiffusion3 --help` (authoritative for flags). prosapia's bundled rfdiffusion3 is a wrapper of the original: check the upstream repo for model behaviour (https://github.com/RosettaCommons/foundry/tree/production/models/rfd3/docs). A project-local `tools/rfdiffusion3/` shadows the bundled tool: its only functional difference is a pinned install in `modal_image.py` (`rc-foundry[all]==0.2.0`, `atomworks<3`, with build-time import checks); delete the fork only after building the built-in's unpinned image and confirming `import rfd3.inference.input_parsing` succeeds in it.

## Premise, and what it does not establish

rfd3 produces **backbone coordinates**. That is the whole claim: no sequence, and nothing it writes says whether the backbone is designable, soluble, or placed where a binder could physically sit. Every quality judgement happens downstream (sequence design, a predictor, the fold and pose gates), so treat a batch as **candidates, not results**. The collected columns are provenance (which batch, which model, which length), not scores; `rfd3_ca_rmsd_to_input` is the only number with evaluative content and it is **NaN for every de novo run** (no input to deviate from).

## Three facts that change what your run means

Each produces a successful-looking run that did something other than intended.

**1. `/0` is the only chain-break literal, and the colon form fails silently.** `elif part == "/0"` is the only break the parser recognises. The colon form `A18-20:B18-20` does **not** error: it matches the prefix only and **drops chain B**, so you design against half the target and every interface number afterwards describes a different problem.

**2. `--length min-max` does not vary length within a batch.** rfd3 draws **one** length per batch, so `--num-designs 5` gives five backbones of the *same* length; believing you ran a length sweep, you ran five replicates. For a spread use several batches (`--num-designs 1 --set n_batches=5`) or several design keys, at the cost of loading the model more often. The drawn value lands in `extra["sampled_contig"]`.

**3. rfd3 renumbers every output chain from 1.** A target numbered 18–155 comes back as 1–138. Every downstream step naming target residues (hotspot recall, epitope contacts, a residue mask, an RMSD on selected residues) inherits that offset and nothing in the table announces it. Hotspots in the *input* spec use the input file's own numbering; everything read off the *output* does not.

## Two shapes of run

- **Root run (no `-t`)**: starts a fresh lineage in `table0`. One design group, named `denovo_diff` (pure de novo) or `<pdb-stem>_diff` (with `--input-pdb`); rows collect as `<group>_0`, `<group>_1`, …
- **Child run (`-t <table>`)**: one design group per ready row, inputs from `-i/--input-column`; mints a child table.

`{expr}` placeholders in contigs/length resolve up the lineage, so they need `-t`: **a root run must use literal values.**

## De-novo invocation

```bash
sapia run rfdiffusion3 <run_dir> --length 80-120 --num-designs 5
```

No `-t`, no contig, no input. Typical: 5 backbones in `table0` in ~3 min on an A100 (≈53 s inference, the rest container start and model load).

## Flags

| Flag | Default | Note |
| --- | --- | --- |
| `--length` | none | `N` or `min-max`. **One length drawn per batch** (fact 2). One subunit's length when symmetric |
| `--contigs` | none | rfd3 contig syntax: `A40-60` motif, bare `30` designed, `/0` chain break |
| `--num-designs` | 1 | Designs per input key (`diffusion_batch_size`); the model loads once for all |
| `--shard-size` | 10 | Design keys packed per task; they run **in series** in one process |
| `--symmetry` | none | `C5`, `D4`, or `auto`. Write **one asymmetric unit's** contig (the sampler replicates it); there is no `--replicate` |
| `--input-pdb` | none | Root-only, for diffusing a structure not yet in a table |
| `--num-timesteps` | 200 | |
| `--step-scale` | 1.5 | Higher = less diverse, more designable |
| `--partial-t` | none | Partial diffusion noise, ~5–15 Å |
| `--extra-spec` | none | YAML/JSON of extra rfd3 spec fields merged into every design |

Default Modal resources: **A10**, 8 CPU, 32 GiB, 4 h; override with `--modal-gpu`. On vib: one GPU, 8 CPU, 32 GiB, 4 h per task.

## Binder design: hotspots and the two-target contig

There is no binder flag. A binder job is a **root run** that holds the target chains as motif and appends one designed chain, with hotspots supplied through `--extra-spec`.

```bash
sapia run rfdiffusion3 <run_dir> \
    --input-pdb target/<target>.pdb \
    --contigs 'A18-155,/0,B18-155,/0,70-100' \
    --extra-spec target/hotspots.yaml \
    --num-designs 2 --set n_batches=4 --modal-gpu A100
```

Typical: 8 backbones, 3 chains each (`A` 138 / `B` 138 / **`C` = the binder**), ~4 min on an A100 for a 276-residue motif plus a ~90-residue binder.

### `select_hotspots`: the syntax

Not documented anywhere in prosapia, and **classic RFdiffusion's `ppi.hotspot_res=[A30,A33]` list form does not carry over.** In rfd3 the field is typed `Optional[InputSelection]` on `DesignInputSpecification`, and `InputSelection.from_any` accepts only **a contig-style string, a bool, or a dict**; a list raises `ValueError: Cannot convert <class 'list'> to InputSelection`.

```yaml
# hotspots.yaml: a whole-residue selection (ALL atoms of each residue)
select_hotspots: "A47,A123,B42,B155"     # ranges work too: "A47-49"
infer_ori_strategy: hotspots             # places the origin token 12 Å out along the outward normal from the hotspot COM
```

```yaml
# atom-level form: the dict picks which atoms carry the annotation
select_hotspots:
  A47: TIP        # sidechain tip atoms
  B42: BKBN       # backbone
  B155: [CA, CB]  # explicit atom names
```

The string-vs-dict choice is what the field's "atom-level or token-level" docstring refers to.

**Specify few hotspots.** rfd3 was trained with hotspots present in 75% of PPI examples, showing only a random subset of up to **20%** of the true hotspot atoms (ground truth = target atoms within 4.5 Å of the binder), so a dense patch is off-distribution. It is real conditioning, not a no-op (the shipped checkpoint carries trained `token_initializer.…is_atom_level_hotspot.weight` tensors).

- Hotspots need `input` set in the same spec and must lie **inside the contig's motif ranges**: annotations on residues outside the contig never enter the built structure.
- Chain/residue ids are the **input file's own numbering**, not renumbered (the *output* is: fact 3).
- prosapia resolves `{expr}` in every `--extra-spec` string, so avoid literal braces.
- Setting `contig`/`length`/`input`/`symmetry`/`partial_t` in both a flag and `--extra-spec` is a `SpecConfigError`. The task script prevalidates, so a bad spec fails fast and cheaply: **let it, rather than building a probe container to check.**

### Contig traps for a two-chain target

- **Use `/0` for the chain break** (fact 1).
- **Do not pass `--length` with a binder contig.** It is the *total*, and `length_min -= num_motif_residues`, so `--length 400` against a 276-residue motif demands a 124-residue binder and fails validation; let the contig's `70-100` govern.
- **The designed length is drawn once per spec build** (fact 2).

## Gotchas

- **There is no checkpoint flag.** On vib the rfd3 activation script exports `RFD3_CKPT` and `rfdiffusion3.sh` passes it as `ckpt_path` at task time; each task's `.out` has a `checkpoint:` line naming it, and that line is the ONLY record of which model made the backbones (not a table column). `--set ckpt_path=...` overrides it for one run. Never put `RFD3_CKPT` in `.env` (Modal tasks receive it). *Signature:* a task dies after ~70 s with `Invalid checkpoint: rfd3` when no checkpoint reached rfd3.
- **Checkpoints live on a Volume (Modal)** (`SAPIA_MODAL_VOLUME_RFD3_CKPT`, mounted at `/checkpoints`) and are **not installed automatically**. If a run fails for a missing checkpoint the Volume is empty: say so rather than retrying. Populating it is a one-off `foundry install rfd3 --checkpoint-dir /checkpoints` from a container of the tool's image (~2.5 GiB, ~3 min).
- Expect harmless stderr noise: an `atomworks` warning about an unset env var, and a log line showing Modal's internal `/__modal/volumes/...` path. Files still land under `/runs`.
- **A large fixed context needs the 80 GB card (Modal).** `--modal-gpu A100` can give the 40 GB A100: a motif of about 940 residues (e.g. 3 × 152 target + 2 × 80 + 2 × 75 fixed copies + a 170-residue fused chain with one linker diffused between two fixed proteins) dies in the first forward pass with `torch.OutOfMemoryError ... GPU 0 has a total capacity of 39.52 GiB`; `--modal-gpu A100-80GB --num-designs 1 --set n_batches=4` runs clean (about 74 s per batch). A failed run leaves a phantom registered table with no `.tsv` and an empty leaf dir: delete only that leaf dir and rerun with the SAME `--table-label` to land in the same registered table.
- **Linking two fixed proteins with one diffused segment works with the contig alone**: motif chains `X1-n,/0,...`, then the fused chain last as `D1-80,<L>,S1-75` (no `/0` between the motif pieces and the designed number). rfd3 renumbers the output chains (identify them by length) and keeps the fixed segments to about 0.1 Å CA RMSD, but it does NOT avoid the other chains: diffused linkers can clip a neighbouring chain, so gate the new segment's sterics (`chainsel` cut + `epitope`, bb+CB, before any sequence design).

## What it collects

Columns (leaf-prefixed `rfdiffusion3_`):

| Column | Meaning | How to read it |
| --- | --- | --- |
| `_status` | `OK` is the only proof the design succeeded | check before anything else |
| `_path` | the **PDB** collect converted from rfd3's `.cif.gz`, under `<run_dir>/.cif_to_pdb/` | the column to feed downstream |
| `_iteration`, `_rfd3_batch`, `_rfd3_model` | provenance: which draw, which batch, which model index within the batch (NOT which checkpoint: the table does not record it) | tell replicates apart from genuine variants |
| `_rfd3_ca_rmsd_to_input` | CA RMSD to the input structure | **NaN for de novo**; meaningful only for motif scaffolding or partial diffusion |
| per-chain length | the drawn length(s) | the actual length, which `--length min-max` only bounded |

There are **no quality columns here by design**: nothing in this table says a backbone is good.

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

Confirm the exact length column name from the table header before writing the first filter: it is per-chain and the suffix depends on the run shape.

## What it is blind to

- **Designability.** Nothing here says a sequence exists that folds to this backbone: that is the MPNN-then-predictor round trip and the fold gate in `usalign --mm 0`.
- **Composition and solubility.** A 4.5%-charged, 70%-TSVIG β-sandwich can fold at pLDDT 0.96: check % charged, net charge, % aromatic, % β-branched directly on the designed sequence; the backbone cannot show it.
- **Whether the hotspots were used.** Hotspot conditioning is real, but no column reports which target residues the binder ends up touching, and the generated pose is not the predicted pose. Measure hotspot recall on the **prediction** (`epitope`, `ringfit`), never on the rfd3 output.
- **Whether the pose is physically possible.** A binder placed against two trimmed protomers may occupy space a symmetry-related protomer fills in the real assembly: `ringfit` against the full assembly answers that; rfd3 will generate it happily.
- **Its own diversity.** One length per batch means a batch can look like a varied ensemble and be replicates: use `_rfd3_batch` to tell which you have.
