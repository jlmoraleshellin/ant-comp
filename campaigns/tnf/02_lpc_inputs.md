# LPC inputs: what Proteina-Complexa needs, and the entries for this campaign

> **Part 2's hotspot list is SUPERSEDED.** It shows the 20-residue draft. The
> campaign uses the 6-residue trimmed set decided in Part 3 —
> `A32, A33, A145, C87, C91, C92`, footprint 24.6 A. The live entries are in
> `targets/targets_entry.yaml` and are already installed in the LPC checkout.
> Parts 1, 3 and 4 are still current.

**Date:** 2026-10-07
**Epitope:** A/C seam, conserved core (decided from `01_recon_human_mouse.md`).

---

## Part 1 — What LPC needs, in general

Five things. Only the last two are campaign-specific; the first three are
one-time setup.

### 1. The environment

```bash
git clone https://github.com/NVIDIA-BioNeMo/Proteina-Complexa
cd Proteina-Complexa
./env/build_uv_env.sh --minimal
```

`--minimal` skips JAX, ColabFold and tmol. Those are the AlphaFold2 / force-field
**reward** models. We are validating with Boltz instead, so we do not need them,
and skipping them avoids downloading AF2 parameters and a RoseTTAFold3 checkpoint.

**Consequence to understand:** the rewards are what drive LPC's inference-time
*search* — the "test-time compute" the paper is named after. With `--minimal` you
get the base generative model and no search. That is the cheap starting mode. It
is not the configuration that produced the published numbers.

### 2. Model checkpoints

```bash
complexa init        # writes a .env template
complexa download    # interactive wizard, pulls weights from NVIDIA NGC
```

Two files matter for protein-binder design:

| file | what it is |
| --- | --- |
| `complexa.ckpt` | the main flow-matching generative model |
| `complexa_ae.ckpt` | the autoencoder for the latent side-chain/sequence space |

Their locations are set in the pipeline config as `ckpt_path`, `ckpt_name` and
`autoencoder_ckpt_path`.

### 3. `.env`

Written by `complexa init`. Configs reference it as `${oc.env:VAR}`. With
`--minimal` most entries are unused; the ones that still matter are
`LOCAL_CODE_PATH`, `COMMUNITY_MODELS_PATH` and `DATA_PATH`.

If a variable a config needs is missing, Hydra fails at startup with
`InterpolationKeyError` naming it. That is deliberate — it is the error you want,
not a silent default.

### 4. A target structure — **a `.pdb` file, not `.cif`**

Your downloads were mmCIF. Converted, cleaned (ligands, waters, hydrogens and
alternate conformations removed) and written to:

```
campaigns/tnf/targets/1tnf_human.pdb    chains A/B/C, residues 6-157
campaigns/tnf/targets/2tnf_mouse.pdb    chains A/B/C, residues 9-157
```

Script: `campaigns/tnf/scripts/mkpdb.py`.

### 5. A target entry in `configs/targets/targets_dict.yaml`

This is the actual campaign definition. Six fields:

| field | meaning |
| --- | --- |
| `target_path` | path to the `.pdb` |
| `target_input` | **which chains and residue ranges the model sees.** Comma-separated `<chain><first>-<last>`. Everything outside this is invisible to the model. |
| `hotspot_residues` | **where to aim.** A list of `<chain><number>`. The model steers the binder onto these. May span several chains — that is how a composite epitope is specified. |
| `binder_length` | `[min, max]` residues for the designed binder |
| `target_filename`, `source`, `pdb_id` | bookkeeping only |

---

## Part 2 — The entries for this campaign

Paste both into `configs/targets/targets_dict.yaml`, under `target_dict_cfg:`,
keeping the two-space indentation.

```yaml
  90_TNFa_HUMAN_xreact:
    source: tnf_xreact
    target_filename: 1tnf_human
    target_path: /home/lpignatti/ant-comp/campaigns/tnf/targets/1tnf_human.pdb
    target_input: "A6-157,B6-157,C6-157"
    hotspot_residues: ["A32", "A33", "A65", "A67", "A113", "A144", "A145", "A147",
                       "C82", "C87", "C91", "C92", "C93", "C94", "C95", "C96",
                       "C123", "C124", "C125", "C126"]
    binder_length: [60, 110]
    pdb_id: "1tnf"

  91_TNFa_MOUSE_xreact:
    source: tnf_xreact
    target_filename: 2tnf_mouse
    target_path: /home/lpignatti/ant-comp/campaigns/tnf/targets/2tnf_mouse.pdb
    target_input: "A9-157,B9-157,C9-157"
    hotspot_residues: ["A32", "A33", "A65", "A67", "A113", "A144", "A145", "A147",
                       "C82", "C87", "C91", "C92", "C93", "C94", "C95", "C96",
                       "C123", "C124", "C125", "C126"]
    binder_length: [60, 110]
    pdb_id: "2tnf"
```

### Why these values

**The hotspot list is identical for both species.** This is the useful accident of
these two PDB entries: 2TNF is numbered to align with human, skipping 73 rather
than renumbering after the deletion. So the same residue numbers mean the same
positions in both. Verified — all 20 exist and are the same amino acid in both
structures (`scripts/checkhs.py`, 0 problems).

**`target_input` differs between them** only because mouse chain A is modelled
from residue 9 and human from residue 6.

**`binder_length: [60, 110]`** is narrowed from LPC's shipped `[50, 120]`. The
epitope is a composite one spanning two protomers; very short binders will not
reach across it. Adjust once you see real designs.

---

## Part 3 — Two things to decide before running

### The hotspot list is large

LPC's shipped targets use **1-7** hotspots. This list has **20**. Hotspots are a
steering signal, not a contact requirement, and over-specifying may either
over-constrain generation or dilute the signal across too wide an area.

The chain-A core alone already spans **32.6 A**; adding twelve chain-C residues
widens it further, while a 60-110 residue binder covers roughly 20-30 A. No single
binder will touch all twenty.

Two options:

- **(a) Run as listed.** Treat the 20 as "this whole seam is fair game" and let the
  model choose where to sit. Simple; costs nothing extra to try.
- **(b) Trim to a ~25 A core**, e.g. `A65, A67, A113, C87, C91, C94` — the mutually
  closest conserved residues spanning the seam. Tighter, more reproducible, more
  likely to give a consistent binding mode.

My suggestion: **run (a) first as the cheap control**, see where designs actually
land, then trim to (b) for the real campaign. Same logic as the staged plan.

### These are generation hotspots, not the cross-reactivity mechanism

Running the human target and the mouse target separately gives two **independent**
sets of binders. It does not give one sequence that binds both.

Dual-species binding comes from the **atomium tied-positions step** (Option B):
one binder backbone placed against both targets, with binder positions tied across
the two copies, so one sequence must satisfy both environments. These target
entries are the input to that; they are not a substitute for it.

---

## Part 4 — Commands

```bash
cd Proteina-Complexa
source .venv/bin/activate

# check the config resolves and every path/env var exists (no GPU needed)
complexa validate design configs/search_binder_local_pipeline.yaml \
    ++generation.task_name=90_TNFa_HUMAN_xreact

# same for mouse
complexa validate design configs/search_binder_local_pipeline.yaml \
    ++generation.task_name=91_TNFa_MOUSE_xreact
```

`validate` is the Stage 0 test: it resolves the whole Hydra config and fails
loudly on a missing path or variable, without touching a GPU. Run it on both
entries before any Modal spend.

Generation itself will **not** run on this laptop — the 4 GB RTX 2050 cannot hold
an all-atom model of a ~500-residue complex. That step goes to Modal.

---

## Status

- [x] Targets converted to PDB and cleaned
- [x] All 20 hotspots verified present and conserved in both species
- [x] Target entries written
- [ ] Decide hotspot list (a) or (b)
- [ ] LPC environment built
- [ ] `complexa validate` passes on both entries
- [ ] Modal configured
