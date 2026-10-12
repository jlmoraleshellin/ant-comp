---
name: protonpottsmpnn
description: >-
  Proton-PottsMPNN, a PottsMPNN with an explicit protonation-state alphabet that designs pH-switchable binder sequences by pinning protonated centres (HIS-P/ASP-P/GLU-P) and redesigning around them. Load before composing a protonpottsmpnn run or reading its columns. Traps: `--explicit-centers` has a residue-numbering trap, `--placement-region interface` needs a complex, check `centers_verified`, and energies never compare across backbones.
---

# protonpottsmpnn

**Custom tool** (lives in `tools/protonpottsmpnn/`, not bundled with prosapia).
**`action: create`** — mints a child table, one row per designed sequence.
**CPU only** — `gpus_per_task` is forced to `0`, so you never pass `-g 0`.

Proton-PottsMPNN (Jacobsen et al. 2026, [repo](https://github.com/christian-creator/ProtonPottsMPNN))
is **PottsMPNN with an explicit protonation-state alphabet**. Each titratable residue is
two tokens, not one:

```yaml
HIS:  HIS-P  # protonated, +1          HIS-S  # neutral
ASP:  ASP-P  # protonated, neutral COOH  ASP-D  # deprotonated, −1
GLU:  GLU-P  # protonated, neutral COOH  GLU-D  # deprotonated, −1
          # HIS-A / ASP-A / GLU-A = ambiguous; this tool forbids them outright
```

Because the learned Potts energy is protonation-aware, the engine can **pin** protonated
centres and redesign their neighbourhood so binding **switches with pH** — the mechanism
is that a residue which is charged at pH 5 and neutral at pH 7.4 makes or breaks the
interface as the pH moves.

It sits in the same slot as `proteinmpnn` / `atomium`: backbone in, designed binder
sequences out, one child row each. What differs is the **objective**.

## The one idea you must hold: lambda

Every design minimises a weighted sum of two z-scored terms:

```yaml
O: (1 − λ)·z(H_stab)  +  λ·z(Σ selective)
#       │                      │
#       │                      └─ selective_energy: how much the centres PREFER being protonated
#       └───────────────────────── potts_energy: how stable the whole thing is
# λ = 0 → pure stability (an ordinary good binder, no switch)
# λ = 1 → pure pH selectivity (a strong switch that may not fold)
```

`--num-designs N` walks λ from `--lambda-min` to `--lambda-max` in N steps. **One run
gives you the whole trade-off curve**, and the collector marks which designs sit on the
Pareto front of that curve. This is the repo's own sweep (`inference/design_ph.py`), not
something invented here.

```yaml
# --num-designs 8 on one backbone
λ=0.00  ●  pareto ✔   stablest, least selective
λ=0.14  ●  pareto ✔
λ=0.29  ○  pareto ✘   dominated — another design beats it on BOTH axes
λ=0.43  ●  pareto ✔
…
λ=1.00  ●  pareto ✔   most selective, least stable
```

## Verified invocation

```bash
# the common case: place centres at the interface, walk the whole front
sapia run protonpottsmpnn <run_dir> -t table0 --table-label ph \
    --binder-chain A --num-designs 8 --placement-region interface
sapia collect protonpottsmpnn <run_dir> -t <the table the run reserved>
```

`default_input_column` is **`rfdiffusion3_path`** — same choice as `atomium`, because the
bundled proteinmpnn's default is the older `rfdiffusion_path` (no 3) and silently submits
nothing after an rfd3 run. From a BindCraft2 table pass `-i bindcraft2_path` (the
complex) or `-i bindcraft2_traj_path`.

## Premise — when its numbers are meaningless

**The input must be a binder + target complex.** `--binder-chain` names the binder; every
other chain is held fixed and seen by the encoder as context. That is what makes this a
*binder* redesign.

| Situation | What happens |
| --- | --- |
| Single-chain input + `--placement-region interface` | **Refused at submit.** The interface mask is built from binder-CA/target-CA contacts, so it would be empty. |
| Single-chain input, any other region | Runs, but it is a **monomer** redesign. The pH switch is then about folding, not binding. |
| `--binder-chain` naming a chain that is not in the file | **Refused at submit**, per design, with the chains it actually found. |
| Poly-glycine backbone + `--placement-region core`/`surface` | **Refused at submit** — no side chains means no RASA split, so the region is empty. |

**Energies never compare across backbones.** `potts_energy` and `selective_energy` are
model energies z-scored *within one featurised backbone*. Ranking the whole child table by
`potts_energy` ranks nothing. **Filter and rank within a parent.**

**This tool designs; it does not fold and it does not measure binding.** A low
`selective_energy` is the model's claim that the protonated microstate is preferred at the
pinned centres. It is not evidence of a pH switch. The chain afterwards is unchanged:
`mkcomplex` → `boltz` → `usalign` / `cms` / `pyrosetta`.

## HBPLUS — the build dependency the upstream README denies

**The repo's README table says design does not need HBPLUS. That is wrong**, and it is the
single most expensive thing to learn late.

`prepare_potts_input` builds its inference pipeline through
`get_protonation_state_transforms`, which runs `CalculateHbondsPlus` **unconditionally**
(`pipelines/potts_mpnn.py:269`). The v6 vocabulary's FLAML labeller consumes those bonds as
its feature pool. **Every design call shells out to the HBPLUS binary.**

HBPLUS is a C program by Ian McDonald (UCL/EBI). Not on PyPI (404), not on conda-forge
(no such package), and the EBI download sits behind a licence form. It is therefore
**vendored** into the tool:

```
tools/protonpottsmpnn/vendor/hbplus.tar.gz     # sha256 937467447bd2e429… (148 KB)
tools/protonpottsmpnn/vendor/README.md         # how to obtain it, and the citation duty
```

**The tarball is git-ignored and must stay that way.** This repo is public; HBPLUS ships
under a confidentiality agreement (`confid.txt`) requiring it be kept in confidence and
"in a reasonably secure place to prevent unauthorised access". Each person who runs this
tool supplies their own copy. Publications must cite McDonald & Thornton (1994) *JMB*
238:777-793 — clause 2 of that agreement.

`modal_image.py` untars, `make`s and exports `HBPLUS_PATH`. With the tarball absent, the
image build fails immediately with that instruction — deliberately, because the
alternative is far worse:

> If `HBPLUS_PATH` is unset, `bond_annotation.calculate_hbonds` falls back to a
> **hardcoded path on the original author's laptop** (`/Users/chrjac/…/hbplus`). The
> failure is then a `FileNotFoundError` on a stranger's OneDrive directory, hundreds of
> frames deep. The worker therefore checks `HBPLUS_PATH` by hand, first thing, and fails
> with a sentence that explains itself.

It is pure C (the bundled `accall.f` belongs to a different program), builds from gcc +
`-lm`, and needs no runtime data files — only an optional `$HOME/.hbplusrc`.

### The second missing-library trap: OpenMP

**Measured.** With HBPLUS present but no OpenMP runtime, the task dies at
`AnnotateProtonationStates` with *"XGBoost Library could not be loaded"* — i.e. **after**
HBPLUS has already succeeded, so it reads like a protonation-labelling bug rather than a
missing system package. The v6 labeller is a pickled FLAML/xgboost model, and xgboost
links `libgomp` at import. `modal_image.py` names `libgomp1` explicitly for this reason.

## Flags

| Flag | Default | Note |
| --- | --- | --- |
| `--binder-chain` | `A` | The chain to redesign. Everything else is fixed target context. Checked against the real structure per design. |
| `--num-designs` | 8 | Designs per backbone = points on the λ ladder. |
| `--samples-per-design` | 1 | Stochastic draws **at each λ**. Total rows per backbone = `--num-designs × this`. Only useful with a raised `--temperature`. |
| `--lambda-min` / `--lambda-max` | `0.0` / `1.0` | Ends of the ladder. Set equal to run every design at one fixed λ. |
| `--center-types` | `HIS-P,ASP-P,GLU-P` | The exact **composition** of centres to place (one of each by default). Repeat to place several: `'HIS-P,HIS-P'`. Only the three protonated states are valid — the deprotonated ones are what each is *scored against*, never placed. |
| `--explicit-centers` | `""` | Pin centres by hand: `'45:HIS-P,78:ASP-P'`. Overrides `--center-types`. Supports `{expr}` — see below. |
| `--placement-region` | `all` | Comma-separated subset of `all,interface,core,surface`. `interface` = binder residue with a target CA within 6 Å; `core`/`surface` split on residue RASA at 0.2. |
| `--placement-by` | `scan_potts` | How candidates are ranked: `scan_potts` (Potts gap), `scan_mpnn` (decoder log-lik gap), `random`. |
| `--neighbour-k` | 16 | Per-centre cap on designable kNN neighbours — the redesign *extent*. 0 = full neighbourhood. |
| `--max-mutations` | 20 | Hard cap on total designable positions = a direct mutation budget. 0 = uncapped. |
| `--block-size` | 3 | Block-descent block size. Cost grows as `V**block_size`: 1 greedy, 2 pairwise, 3 triples. |
| `--temperature` | 0.05 | Block readout temperature. Near-deterministic by default. |
| `--seed` | 0 | **`0` is a real seed here**, unlike atomium where 0 means "pick randomly". |
| `--seed-column` | `""` | Redesign *from* a sequence column (`proteinmpnn_sequence`, `atomium_sequence`) instead of the native one. |
| `--checkpoint` | shipped v6 | Must be an `extended_vocab="v6"` checkpoint or the 30-token weight load fails. |
| `--n-jobs` | 0 (= `--cpus-per-task`) | CPU workers the λ ladder is fanned across *within* one task. The backbone is featurised once and shared read-only via fork — close to free. |
| `--set FIELD=VALUE` | — | Escape hatch onto any other `PHDesignCriteria` field, repeatable. An unknown field fails the task loudly. |

Default Modal resources: **no GPU**, 8 CPU, 16 GiB, 2 h timeout. All weights (21 MB Potts
checkpoint + 43 MB labeller folds) ship inside the image — no cache Volume.

### `--explicit-centers` and the residue-numbering trap

Residue numbers are the **input structure's own `res_id`**, not 1..L sequence positions.
The engine matches centres through `res_id_to_pos`.

```bash
--explicit-centers '45:HIS-P,78:ASP-P'     # literal res_ids
--explicit-centers '{motif_end}:HIS-P'     # resolved up the lineage, per design
```

Generators routinely renumber a chain from 1. If rfd3 handed you a binder numbered 1..90,
then `45` means *the 45th residue of the binder*, not residue 45 of your target. The
collector reports **`resnum_offset`** (res_id of the first binder residue, minus 1) so you
can tell which world you are in — `0` means renumbered from 1.

### `--placement-region` fails soft upstream, loud here

The engine classifies regions inside `try: … except Exception: pass`, so a failed SASA or
contact computation leaves an **all-False mask** and placement quietly has nowhere to go.
The worker checks every requested region for a non-empty mask and **raises** instead,
because the alternative is a design with centres somewhere you did not ask for and no
indication that anything went wrong.

## What it collects

Child rows named `<parent>_p0`, `<parent>_p1`, … (the `_p` suffix distinguishes them from
proteinmpnn's `_f` and atomium's `_a`). Row order is the engine's `sorted_by_energy`, so
**`_p0` is the stablest design, not the most selective one** — rank by a column, never by
the suffix.

Columns, leaf-prefixed `protonpottsmpnn_`:

```yaml
science:
  sequence:          # 1-letter binder sequence — what you fold
  extended_tokens:   # the parallel 3-letter + state string (… ASP-P … HIS-P …)
  potts_energy:      # whole-system Potts H — LOWER = stabler
  selective_energy:  # Σ(e_P − e_D) over centres — LOWER = more pH-selective
  global_dh:         # binder-wide H(all protonated) − H(all deprotonated)
  combined_lambda:   # the λ this design was optimised at
  pareto:            # True = not dominated on BOTH energies, within this backbone
  centers:           # 'resnum:STATE;…' in INPUT-structure numbering
  center_seqpos:     # the same centres as 1-based positions in `sequence`
  n_centers:         #
trust:               # read these first
  centers_verified:  # fraction of pinned centres whose output token IS the pinned state
  resnum_offset:     # res_id of first binder residue − 1. 0 = renumbered from 1
  n_mut:             # mutations vs the sequence the redesign STARTED from
  seq_rec:           # identity to that same starting sequence
  seed_source:       # 'native' | 'inverse' (--seed-column)
  binder_chain:      #
  binder_len:        #
  n_designable:      # positions the optimiser was allowed to touch
```

Plus `_status` and `_path` (the per-backbone `designs.tsv`).

### `centers_verified` is the column that decides whether to trust the row

It is the fraction of pinned centres whose token in `extended_tokens` really is the
protonated microstate that was pinned. **Below 1.0, the design does not carry the centre it
was optimised for, and its `selective_energy` is meaningless.** Always gate on it:

```bash
-f "protonpottsmpnn_status == 'OK' and protonpottsmpnn_centers_verified == 1.0"
```

### `n_mut` can EXCEED `--max-mutations`, by up to `n_centers`

**Measured.** `--max-mutations 20` with 3 centres produced `n_mut` of 21 and 22.
That is not a bug: `--max-mutations` and `n_designable` cap the **designable
neighbourhood** around the centres. The pinned centres themselves are mutations too,
and they are not in that set.

```yaml
λ=1.0 run on the PD-L1 example, --max-mutations 20, 3 centres:
  n_designable: 20            # the neighbourhood the optimiser could touch
  centres:      8 T→D  ·  17 A→E  ·  28 G→H     # 3 more changes, outside that cap
  n_mut:        22            # 19 neighbourhood + 3 centres
```

So the real bound is `--max-mutations + n_centers`. If you are filtering to stay close
to a validated parent sequence, filter on `n_mut` itself — do not assume
`--max-mutations` is the ceiling.

### `n_mut == 0` means nothing was designed

The seed was echoed back. Usually `--max-mutations` or `--neighbour-k` was too tight, or
placement found no candidate position.

## Reading the numbers

Both energies are **lower = better**, and both are **arbitrary units, per backbone**.

`-f/--filter` takes a **Python module defining `apply_filter(df) -> df`**, not an inline
expression. Write one and point `-f` at it:

```python
# filters/ph_switch.py
def apply_filter(df):
    return df[
        (df["protonpottsmpnn_status"] == "OK")
        & (df["protonpottsmpnn_centers_verified"] == 1.0)   # centre really is in the sequence
        & (df["protonpottsmpnn_selective_energy"] < 0)      # protonated state is preferred
        & (df["protonpottsmpnn_pareto"])                    # on its backbone's front
        & (df["protonpottsmpnn_n_mut"] > 0)                 # something was actually designed
    ]
```

```bash
sapia run mkcomplex <run_dir> -t table1_ph -f /runs/filters/ph_switch.py
```

`selective_energy < 0` is the sign test that the protonated state is preferred at all. How
much more negative is "good" is **not calibrated** — the manuscript ranks within a
campaign, and so should you.

## Companion documents

| Document | For |
| --- | --- |
| `docs/protonpottsmpnn-modal-guide.md` | the Slack-shareable how-to: quickstart, expected outputs, reading the table |
| `docs/protonpottsmpnn-build-notes.md` | why it is built this way, what was measured, known gaps. **Read this if asked "how does it work" or "why this design".** |
| `/runs/tests/protonpottsmpnn/` (Volume) | the PD-L1 test input, `seed_table0.py`, `ph_switch.py` |

## Verified end to end

Run on the repo's own `inference/examples/pdl1_seed_binder.pdb` (PD-L1 seed binder,
chain A = 114 aa binder, chain B = target), centres `ASP-P,HIS-P`, λ ∈ {0.0, 0.5}:

| λ | potts_energy | selective_energy | n_mut | seq_rec | centers_verified | pareto |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0 | −51650.0 | −13.99 | 8 | 0.930 | 1.0 | ✔ |
| 0.5 | −51624.3 | **−21.11** | 11 | 0.904 | 1.0 | ✔ |

**The trade-off behaves as advertised**: raising λ bought 7 units of selectivity and paid
26 units of stability. Both designs put `ASP-P` at position 33 and `HIS-P` at 100, and the
two sequences differ *only* at positions 23–35 and 95–98 — the centres' neighbourhoods,
as `--neighbour-k 16` / `--max-mutations 20` require. `resnum_offset` came back **0**:
this binder is numbered 1..114, so its `centers` are sequence positions, not target
numbering — exactly the case that column exists to expose.

### Also verified on Modal

Same input, same two λ, through the real `sapia run` → `.exit` → `sapia collect` loop on
the `sapia-runs-alvaro` Volume. Task exited `0`; `table1_ph` came back with **26 columns**
and lineage stamped (`parent_name=pdl1_seed`, `gen=1`):

| row | λ | potts_energy | selective_energy | centers | verified | n_mut |
| --- | --- | --- | --- | --- | --- | --- |
| `pdl1_seed_p0` | 0.0 | −51651.0 | −14.62 | `33:ASP-P;100:HIS-P` | 1.0 | 7 |
| `pdl1_seed_p1` | 0.5 | −51631.7 | **−20.62** | `33:ASP-P;100:HIS-P` | 1.0 | 9 |

**Results are not bitwise reproducible across machines.** The Modal numbers differ from the
local macOS run in the third significant figure (−14.62 vs −13.99, 7 vs 8 mutations) at the
same seed. The centres chosen, the direction of the trade-off and `centers_verified` were
identical. The cause is the protonation labeller: a different BLAS/thread count and a newer
xgboost reading the pickles change a few borderline labels, which shifts the energies.
**Compare designs within one run, not across runs on different backends.**

Expect a **benign `.err`** on every task: an xgboost "model generated by an older version"
warning (the FLAML pickles predate xgboost 3.x), a biotite `chararray` DeprecationWarning,
and a `fork() may lead to deadlocks` warning from the CPU pool. Exit code 0 with ~2 KB of
`.err` is the normal, healthy signature here.

### The interface run — the premise this tool is actually for

Same PD-L1 example, `--placement-region interface --num-designs 8` (default 3-centre
composition), on Modal. Task exit `0`, 8 rows, every one `status OK` and
`centers_verified 1.0`, 8 distinct sequences:

| λ | potts_energy | selective_energy | global_dh | n_mut | pareto |
| --- | --- | --- | --- | --- | --- |
| 0.000 | −51650.7 | −3.75 | 24.2 | 15 | ✔ |
| 0.143 | −51645.4 | −5.62 | 24.1 | 14 | **✘ dominated** |
| 0.286 | −51647.3 | −7.31 | 23.9 | 15 | ✔ |
| 0.429 | −51616.0 | −13.90 | 21.9 | 16 | ✔ |
| 0.571 | −51587.2 | −19.98 | 20.2 | 16 | ✔ |
| 0.714 | −51514.9 | −24.25 | 18.5 | 18 | ✔ |
| 0.857 | −51457.5 | −27.05 | 15.5 | 21 | ✔ |
| 1.000 | −51406.0 | −27.78 | 12.4 | 22 | ✔ |

**Monotonic in both directions**, which is the strongest evidence the objective is wired
correctly: selectivity improves 7-fold from λ=0 to λ=1 while stability degrades by ~245
units. `global_dh` falls alongside (24.2 → 12.4) without being optimised at all — the
whole binder drifts toward preferring protonation, not just the three centres.

The one dominated design, λ=0.143, is beaten by λ=0.286 on **both** axes — so the Pareto
column is doing real work rather than marking everything true.

**Interface placement put the centres at 8/17/28**, versus 33/100 for the same structure
under `--placement-region all` — and all mutations clustered at 7–19, 27–31, 38–48,
106–107, i.e. the centres' neighbourhoods. Expect 7 of 8 on the front for a smooth
trade-off like this; a front with far fewer points means the ladder is hitting a wall
(usually `--max-mutations` too tight).

**Timing: ~2 min of CPU per λ** for a 114-residue binder (4 min wall for 2 λ on 2 forked
workers, Apple M-series). The λ ladder parallelises across `--n-jobs`, so `--num-designs 8`
on 8 CPUs costs about the same wall time as 1. Budget per *backbone*, not per design.

## What it is blind to

- **Whether the sequence folds.** No structure is produced. Fold with `boltz` and check
  self-consistency with `usalign`, exactly as for any other designer.
- **Whether the switch is real.** The model scores one backbone in one protonation
  assignment. It does not simulate pH, does not predict pKa shifts in the folded complex,
  and does not model the conformational change a switch implies.
- **Binding affinity at either pH.** Use `cms` / `pyrosetta` on the folded complex.
- **Anything about the target.** Only the binder chain is designed and only its residues
  appear in `sequence`.

## What this tool deliberately does not wrap

The repo has three other halves. None of them are this tool:

| Repo part | Why not here |
| --- | --- |
| `labeller/` | Annotates protonation states of a structure. Different question (annotate, not design) — would be an `update` tool. |
| `scoring/` + `inference/fold_rf3.py` | RF3 folding and pH-bond / charge-clash read-outs. RF3's ~3 GB weights are not shipped, and folding is `boltz`'s job here. |
| `benchmarks/`, `training/` | Not pipeline steps. |

## Traps found while building

- **HBPLUS is required at design time** and the README says otherwise. See above. This is
  the one thing that will stop a first run.
- **`--seed-column` must give the binder chain ALONE.** A complex-sequence column
  (`mkcomplex_sequence`) carries the target too, and the length check rejects it with that
  explanation. A `--seed-column` that is empty on a row is an **error**, not a silent
  fallback to the native sequence — the two give different designs and the table has to
  record which one ran.
- **`_p0` is the stablest, not the best.** Sorted by `potts_energy`, which at λ=1 is the
  term being traded away.
- **The deprotonated states are not placeable.** `--center-types HIS-S` is rejected;
  `HIS-S` is the contrast `HIS-P` is scored against.
- **`--block-size` is exponential.** 4 is not "a bit slower" than 3.
- **The whole λ ladder runs in one task**, fanned over `--n-jobs` forked CPU workers
  sharing one featurisation. Raising `--num-designs` is much cheaper than raising the
  number of backbones.
