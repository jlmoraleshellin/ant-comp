---
name: cms
description: >-
  GPU contact molecular surface (CMS) and shape complementarity (SC) of a binder/target interface via cms-cuda, on Modal. Load before composing a cms run or ranking binders by interface quality. Covers which structure column to score, `--binder-chains`/`--target-chains`, batching, and the GPU/CPU device rule.
---

# cms

**Custom tool** (lives in `tools/cms/`, not bundled with prosapia). Wraps
[cms-cuda](https://github.com/ullahsamee/cms-cuda), the CUDA port of Brian Coventry's
`py_contact_ms` (Longxing Cao's contact molecular surface) plus Lawrence & Colman SC.
GPU results match the original NumPy code to ~1e-13.

**`action: update`**: an interface score is a property of a design that already exists,
so it annotates the table it reads, in place. `-t` is required. **GPU by default** (L4);
`-g 0` switches to the CPU path.

## What it measures

- **CMS** = Σ over the target's interface surface dots of `area · exp(-0.5 · d²)`, where
  `d` is the gap to the binder's surface. It is a distance-weighted interface area in Å², so
  it rewards an interface that is both **large and tight**. It is the interface filter
  used in Baker-lab binder campaigns (Cao et al. 2022). By convention it is the
  **target side** (`cms_target`). The binder side (`cms_binder`) comes from the same
  calculation.
- **SC** is shape complementarity, a single 0–1 statistic for the whole interface.
  Well-packed protein interfaces are ~0.5–0.75 and antibody/antigen interfaces are lower.
  It comes off the same molecular surface at almost no extra cost.

Neither is a binding energy. They measure **geometry only**: atom positions plus a fixed
radius table, with no chemistry, no charges and no solvation. Pair them with `pyrosetta`'s
`if_dG` / `buried_unsat` and with prediction confidence.

## Invocation

```bash
sapia run cms <run_dir> -t table1 -i boltz_path \
    --binder-chains C --target-chains A,B -l pred
# poll .exit as usual, then:
sapia collect cms <run_dir> -t table1 -l pred
```

| Flag | Default | Notes |
| --- | --- | --- |
| `-i/--input-column` | **effectively required** | Sentinel default `"not applicable"` (the usalign/chainsel convention). The builder **raises** unless the table has the column. Reads PDB or mmCIF. |
| `--binder-chains` | **required** | Comma-joined, e.g. `C` or `H,L`. |
| `--target-chains` | **required** | Comma-joined, e.g. `A,B`. Must be disjoint from the binder (checked at submit). |
| `--exclude-resnames` | `""` | Resnames dropped from **both** sides, e.g. lipids or a ligand. Waters, single-atom ions and hydrogens are always dropped. |
| `--no-sc` | off | Skip SC. Rarely worth it. |
| `--max-cms` | off | Adds each side's maximum possible CMS and `frac_*` = CMS / max. Mostly useful for small-molecule targets. |
| `--device` | `auto` | `auto` = `cuda` on a GPU task, `cpu` under `-g 0`. `cuda` never falls back: without a usable GPU, every design gets an error status instead of running ~100× slower. |
| `--designs-per-task` | `100` | Designs per Modal container. See below. |
| `-l/--dir-label` | `""` | **In practice required.** See the traps. |

### Batching: why one task holds many designs

Measured on an L4: the **first design in a container takes ~60 s**, which is the one-time CUDA
kernel compilation (JIT, cached for that container only). **Every later design takes ~0.6 s**,
mostly structure parsing. The builder therefore packs
`--designs-per-task` designs into a sub-manifest (`<out_dir>/cms_tasks/task_<t>.tsv`) and
submits one task per chunk. With the default of 100, 250 designs make 3 tasks, each ~2 min. Do not set it to 1: every design would then pay the minute-long compile. Lower the
value only for more parallelism on very large tables. **`n_tasks` in `cms_modal.json` is
the number of chunks, not the number of designs.**

## Choosing the structure column

CMS and SC are only as meaningful as the pose they are computed on:

| Column | What it tells you |
| --- | --- |
| `rfdiffusion3_path` (table0) | The **designed** interface. rfd3 backbones are poly-Gly/Ala or lack real side chains, so CMS on them is a backbone-level estimate. Use it as a coarse early filter only. |
| `boltz_path` / `alphafold3_path` | The **predicted** complex, with full side chains. **This is the number to rank on**, but only for designs whose prediction actually kept the pose (check `usalign` / `ringfit` first). CMS on a binder that drifted away measures the wrong interface. |
| a pyrosetta relaxed structure | The **refined** interface. It is usually higher and smoother than the raw prediction. Compare it only with other relaxed structures. |

Check chain letters in an actual file before choosing chains. The predictor gives the
binder the **last** letter of the `mkcomplex` entity order, and rfd3 may use another letter
(see the chainsel skill). CMS depends on the binder/target split: a swapped or partial
split gives a plausible, wrong number.

## Traps

**An absent chain is an error, but a wrong chain is not.** A chain in
`--binder-chains`/`--target-chains` that is missing from a design gives
`cms_status = "error: chain X not in <file> (have: …)"`. A chain that is present but
**wrong** (the target named as binder, one protomer forgotten) scores without complaint.
To guard against it, compare `cms_n_atoms_binder` with the binder size you designed (≈ 7–8
heavy atoms per residue).

**A multi-protomer target must be named in full.** Scoring against only protomer `A` of a
two-protomer site drops half the interface and does not error. Pass `--target-chains A,B`.

**CMS = 0 is a real result, not a failure.** It means no binder atom lies within 8 Å of the
target: the binder fell off, or you picked the wrong chains. The status is still `OK`.

**`-l/--dir-label` is effectively required.** The leaf `cms` names both the output dir and
every column, so scoring `boltz_path` and then `rfdiffusion3_path` on the same table with
no labels would overwrite the first run's columns. Use `-l pred` and `-l design` to get
`cms_pred_*` and `cms_design_*`. **Pass the same `-l` to `collect`.**

**Failures still exit 0.** Errors are recorded as data, so `.exit` files that are all `0`
do not mean every design succeeded. Only the `cms_status` column is reliable. A worker
crash writes `error: worker crashed` for every design in that chunk.

**Rows can drop out of the manifest silently.** A design whose input file does not exist
is printed as `MISSING … (skipping)` and left out.

**CPU path memory is quadratic.** `-g 0` runs the original NumPy code, which builds
atom × atom distance matrices. A 4 k-atom target with a 1 k-atom binder was OOM-killed on
a 15 GB machine. If you must use the CPU, raise `--mem` well above the default. On the GPU,
memory is linear and this does not come up.

**`n_radius0_*` > 0 means some atoms are invisible.** Those atoms have no entry in the CMS
radius table (Br, I, B, most metals), so they contribute no surface. For protein/protein
interfaces this should be 0. If it isn't, check what is in the selection.

**Only compare like with like.** CMS grows with interface size, so compare designs against
the same target with similar binder types. Never compare CMS on backbones with CMS on
predictions.

## Reading the numbers

Rough guide for de novo mini-protein binders (target side, on a predicted or relaxed
complex):

| `cms_target` | Reading |
| --- | --- |
| < 200 Å² | Weak or barely touching. Usually discarded. |
| 300–400 Å² | Typical cutoff range in published binder campaigns. |
| > 450 Å² | Large, tight interface. |

These numbers are orientation, not law: the right cutoff depends on target and epitope
size. Rank within the campaign, and read CMS together with SC, `if_dG` and prediction
confidence. Reference values from PDB complexes (the cms-cuda tutorial): p53 peptide/MDM2
485 Å², SC 0.76; VHH/lysozyme 515 Å², SC 0.77.

## Columns collected (`cms_` prefix, or `cms_<label>_` with `-l`)

| Column | Meaning |
| --- | --- |
| `_target` | **CMS on the target side, Å².** The conventional CMS number. |
| `_binder` | CMS on the binder side, Å². |
| `_sc` | Shape complementarity (0–1). NA with `--no-sc`. |
| `_sc_area` | Trimmed interface area of both sides, Å². |
| `_sc_median_dist` | Median interface separation, Å. |
| `_max_target` / `_max_binder` | Max possible CMS of each side alone (`--max-cms` only). |
| `_frac_target` / `_frac_binder` | CMS / max (`--max-cms` only). |
| `_n_atoms_binder` / `_n_atoms_target` | Heavy atoms scored on each side. This is your chain-selection sanity check. |
| `_n_radius0_binder` / `_n_radius0_target` | Atoms invisible to the surface (no radius). |
| `_device` | `cuda` or `cpu`: what actually ran. |
| `_seconds` | Wall time for this design's calculation. |
| `_path` | **Per-residue CMS table** (`side, chain, resnum, resname, cms`) for both sides. These are the interface hotspots. Target-side rows show which epitope residues the binder covers. |
| `_status` | `OK`, `error: <reason>`, or `missing`. This is the only reliable success signal. |

**NA (empty) means not applicable, not zero.** On disk:
`<run_dir>/<table>/cms[_<label>]/<name>.tsv` and `<name>_per_residue.tsv`.

## Verified

- **GPU (Modal L4, this exact image), 2026-09-28:** upstream's examples scored with
  binder `B` and target `A`. linearpeptide gave CMS 471.5166 / SC 0.6415, **identical** to
  the CPU path and to upstream's `calculate_cms_sc.py`. binder.pdb gave 853.03 / 0.622,
  Nanobody 629.74 / 0.659, scFv 1057.80 / 0.465. A repeated structure gave an identical
  result. First design 61 s (kernel compile), then 0.6–5 s per design.
- Error paths (absent chain, missing file, no GPU with `--device cuda`) are recorded as
  error data. They were checked locally.
- Image traps found on the way: the upstream package **cannot be pip-installed** (so the
  pinned source is copied in), and a CUDA *runtime* base image fails with "Failed to find
  CUDA headers" (so it uses *devel*). See `tools/cms/modal_image.py`.
- **The full `sapia run cms` → `collect` loop on the `sapia-runs` Volume has not been run
  yet.** Treat the first real run as a shakedown: run 2 designs, then read `cms_status`,
  `cms_device` and `cms_n_atoms_*` before scaling.
