---
name: atomiumph
description: >-
  AtomiUM's `design_ph.py` for pH-selective binder design: centre modes, required chain flags, the lambda sweep, trust columns and protonation energies. Load before composing an atomiumph run or interpreting its numbers.
---

# atomiumph

**Custom tool** (`tools/atomiumph/`), **Modal only** — no `SAPIA_ACTIVATE_ATOMIUMPH` exists, so a
slurm task stops at activation. **`action: create`** — mints a child table, **one row per designed
sequence** (`<parent>_p0`, `<parent>_p1`, …).

Wraps `design_ph.py` from the private AtomiUM repo (branch `pure_wo_jit`). It is the **sibling** of
`tools/atomium`, not a mode of it: same repo, different entry point, **its own commit pin and image**.
`tools/atomium` designs sequences for a backbone and reports sequence recovery — no protonation
states, no apo reference, no `dE`/`ddE`. `tools/protonpottsmpnn` is a **different program** reporting
Potts energies with **no complex-vs-apo contrast**, which is the whole quantity `ddE` exists to isolate.

## Premise, and when its numbers are meaningless

`design_ph.py` pins titratable binder residues (His/Asp/Glu) into a chosen protonation state and
redesigns their environment so that state is favoured **in the complex more than in the free binder**:

```
dE  = log p(target_state) - log p(comparison_state)      one centre, one context
ddE = dE(complex) - dE(apo binder)                       the selection metric
```

**The input must be a complex.** The apo reference is the same file with the target chains deleted.

| Situation | What happens |
| --- | --- |
| No target chain in the input | Refused at submit (`--target-chains` is required). |
| `--target-chains` names the wrong chain | **It still runs.** `ddE` is computed against the wrong apo reference and reads like a result. **Nothing in the output detects this.** |
| `--binder-chains` names the target | **It still runs, and redesigns the target.** `atomiumph_seq_len` against the binder length you expect is the only cheap tell. |

**There is no fixed binder chain letter — establish it per structure.** The house convention is that
**the binder is the LAST chain**, so which letter that is depends on how many chains the target
carried (two-chain target → `C`; one-chain → `B`; a trimer → `D`). `design_ph.py`'s own default
`--binder_chains A` is wrong here. Establish it from an upstream `chainsel` `chain_n_res`, from an
`epitope` chain list, or from the **pre-flight** (`--lam 0.0 --num-seq 2`), which prints the
chains it found before generating anything. Both chain flags are required and must be disjoint; the
tool asserts nothing about which is which.

**Energies are model log-probabilities on one backbone.** Comparable between draws of the same parent
and across the lambda sweep. They are **not** a calibrated pKa, and `ddE_mean` across unrelated
backbones is weakly meaningful at best — rank within a parent first.

## `--center-mode` — which pH behaviour you get

Getting this backwards produces a healthy-looking table of the opposite molecule. **Required, no
default.**

| `--center-mode` | Pins | Binds | Releases | Charge prior seeks |
| --- | --- | --- | --- | --- |
| `his_s` | neutral His (0) | **neutral** (plasma 7.4) | acid (endosome 6.0–6.5) | **positive** (Arg/Lys) |
| `his_p` | protonated His (+1) | **acid** (tumour 6.0) | neutral (7.4) | **negative** (Asp/Glu) |
| `acid_p` | protonated Asp/Glu | acid | neutral | negative |
| `all_p` | every protonated type | acid | neutral | negative |

- **Recycling / "sweeping" binder** (bind in plasma, dump the antigen in the endosome, anti-TNF) is
  **`his_s`**.
- **Tumour-microenvironment binder** (bind in the acidic tumour, release at neutral, anti-EGFR) is
  **`his_p`**.
- `acid_p`/`all_p` pin Asp/Glu, which sit near pKa 4 and only switch below ~5.8 — **diagnostic only,
  unusable for either His challenge**. That is why upstream's `acid_p` default is not inherited.

**The upstream README's section-1 use-case table has Binds/Releases SWAPPED.** Its own section-4
`--center_mode` table, its section-5 examples and the client README all say the opposite, and swapping
section 1 is what gives biologically correct designs. **The table above is right — do not re-derive
the direction from the README.**

**The `his_s` charge prior is sign-reversed.** The prior adds `-dq * S`, where `S` is the signed sum of
nearby charged **target** residues. `his_p`/`acid_p` have `dq = +1` and seek a **negative**
environment, matching intuition. **`his_s` has `dq = -1` and seeks a POSITIVE (Arg/Lys) environment** —
so "put the histidine next to an aspartate" is exactly backwards for the recycling direction. Target
histidines are skipped by the prior.

**Pre-flight tells you whether the epitope can work at all.** `--lam 0.0 --num-seq 2` prints into the
task `.out`, before generating anything: the chains found, the banned tokens, **the charge-field
range** (`Charge prior on (weight=0.5); field range [-1.51, +1.94]`) and the centres with their scores
(`pos 10 HIS-S placement_score +0.708 model +0.657 charge +0.10`). **The charge-field range and the
model/charge decomposition reach no column** — the `.out` is the only place to read them. If the field
has the wrong sign for the mode, that epitope lacks the complementary charges and no amount of
sampling fixes it.

## Fork with `--table-label`, never with `-l`

`-l/--dir-label` becomes part of the **column leaf**, so it renames every column this tool writes:

```
no -l        →  atomiumph_status, atomiumph_ddE_mean, …
-l nv2_n03   →  atomiumph_nv2_n03_status, atomiumph_nv2_n03_ddE_mean, …
```

So every filter written as `df["atomiumph_status"]` raises `KeyError` against a labelled leaf, and
**`collect` then needs the SAME `-l`** or it looks in `…/atomiumph/`, finds nothing, and reports every
parent as having produced no output. `--table-label` separates both the table and the leaf directory
while keeping plain `atomiumph_` columns. Reach for `-l` only in its one real case: **two leaves
feeding ONE table**, where the rename is what keeps their columns apart.

## Invocation

```bash
# 1. pre-flight ONE structure: prints chains, banned tokens, charge-field range and
#    centres before generating anything.
sapia run atomiumph <run_dir> -t table0 --table-label preflight \
    --binder-chains <binder> --target-chains <targets> \
    --center-mode his_s --lam 0.0 --num-seq 2 -f pick_one_design.py

# 2. the sweep. Only the three flags below are required; the rest is the recipe.
sapia run atomiumph <run_dir> -t table0 --table-label ph \
    --binder-chains <binder> --target-chains <targets> \
    --center-mode his_s          # his_p for a tumour/acid-binding campaign

sapia collect atomiumph <run_dir> -t <the table the run reserved>
#   ...and the SAME -l as the run, if the run used one.
```

`default_input_column` is **`rfdiffusion3_path`** (rfd3 binder output is already a complex). From any
other parent pass `-i`, and check the file still holds **both** chain sets: `-i boltz_path` on a
binder-only prediction, or a `chainsel` output that extracted the binder, has no target and the run is
meaningless even though it will submit.

Resources: **L4, 8 CPU, 16 GiB, 2 h**. The image builds from the private repo at build time using the
`github-token`/`github-username` Modal Secrets; the first submit pays a long build (torch + cu128 +
weights), so allow 15+ minutes.

## Flags

Defaults are README v2's TNF-α production recipe, not `design_ph.py`'s own argparse defaults. The
`client`/`script` columns show what would apply if a flag were simply omitted upstream.

| Flag | ours | client | script | |
| --- | --- | --- | --- | --- |
| `--binder-chains` | **required** | `A` | required | the binder is the LAST chain here |
| `--target-chains` | **required** | `""` | `""` | empty target ⇒ every ddE structurally zero |
| `--center-mode` | **required** | `acid_p` | `acid_p` | a default would answer a question nobody asked |
| `--centers-per-type` | **`2`** | `1` | `1` | README v2 |
| `--redesign-scope` | **`coupling`** | `neighbourhood` | `neighbourhood` | README v2: "the default choice for production" |
| `--interface-cutoff` | **`8.0`** | `12.5` | `12.5` | README v2's strict first shell |
| `--coupling-radius` | **`12.0`** | `15.0` | `15.0` | README v2 |
| `--model-noise` (alias `--noise`) | **`n03_v2`** | `n05` | `n05` | user decision |
| `--apo-weight` | `0.5` | `0.5` | `0.0` | client ≠ script, so it is always passed |
| `--protect-titratable` | `1` | `1` | `0` | client ≠ script, so it is always passed |
| `--lam` | `0.0 0.25 0.5 0.75` | same | same | **must contain `0.0`** |
| `--num-seq` | `8` | `8` | `8` | **PER LAMBDA** |
| `--temperature` | `0.1` | `0.1` | `0.1` | README's example uses 0.15 |
| `--max-mutations` | `20` | `20` | `20` | README's example uses 12 |
| `--seed-mode` / `--seed` | `native` / `42` | same | same | |
| `--charge-prior` / `--charge-cutoff` | `0.5` / `8.0` | same | same | sign follows `--center-mode` |
| `--placement-region` | `interface` | same | same | |
| `--neighbour-k` / `--burial-cutoff` | `12` / `16` | same | same | |
| `--coupling-probe` / `--ban-cys` | `ala` / `1` | same | same | |
| `--output-scope` | `binder` | same | same | |
| `--no-control`, `--set` | — | — | — | ours only |
| `--designs-per-task` | **`5`** | — | — | ours only; structures per task (atomium's is 10) |

Where the meaning is not obvious:

| Flag | Note |
| --- | --- |
| `--centers-per-type` | Centres per **TOKEN TYPE**. `his_p`/`his_s` have one type (so `2` → **2 centres**), `acid_p` two (Asp, Glu → **4**), `all_p` three (→ **6**). Why 2: one histidine gives a transition too shallow to separate pH 6.5 from 7.4 when the window is ≤1.5 units. Centres scale the `coupling` scan linearly and each adds a value to the `;`-joined per-centre cells. |
| `--interface-cutoff` | `interface` is a pure Cα-distance mask with **no burial filter**, and nothing excludes the N-terminal residue — a centre at index 0 is normal. Upstream's 12.5 admits a second shell, so an "interface" centre need not touch the target at all. |
| `--model-noise` | `n03`, `n05`, `n03_v2`, `n05_v2` → `ft_weights_for_pH/ft_model_epoch10_<noise>.pt`, all fine-tuned **30-class** pH checkpoints. `tools/atomium`'s `n00…n07` are the **21-class base** weights and are refused at submit. **`_v2` are DIFFERENT weights, not newer builds** — nothing upstream says what changed in training. |
| `--lam` | The selective-term weight sweep, **per value**. Duplicates and values colliding at 6 dp are refused. λ ≥ 1.0 is degenerate (sampling goes uniform); the upstream client refuses it outright. |
| `--no-control` | Deliberate opt-out of the `0.0` requirement. Only honest when that control already exists for the same parents in another table. |
| `--apo-weight` | Weight of the apo term **in the sampling objective**, not only in scoring. Above 0 the apo graph enters sampling — what pushes a design toward *selectivity* rather than toward merely being good in the complex — and roughly doubles the infill forward passes. |
| `--protect-titratable` | `1` removes every **natively** His/Asp/Glu binder residue from the redesign candidate pool in all scopes, so it stays as in the input. Centres are chosen *before* this, so a protected residue is still eligible as a centre. |
| `--seed-mode` | `native` redesigns the sequence already in the PDB. **`generate` is needed for a raw poly-Gly backbone**, where `native` would seed the model from nothing. |
| `--redesign-scope` | `coupling` ranks candidates by **measured** influence on the centres' `dE` and is the only scope writing `couplings.csv`. `neighbourhood` is the cheap pass (graph neighbours, capped by `--neighbour-k`/`--max-mutations`). `binder` is every designable position, for a raw backbone or an untrusted sequence. |
| `--seed` | Every structure seeded identically, and the pin **reseeds before every lambda arm** — so draw *i* sees the same decoding order at every lambda and differences across lambda are attributable to lambda alone. `0` is not special. |
| `--output-scope` | `all` joins every chain with `/`, which makes `seq_len` meaningless as a binder length and disables the centre sequence check. Leave it at `binder`. |
| `--set` | Raw flag forwarded verbatim, repeatable, **for unmodelled flags only**: `--contact_cutoff` (8.0), `--bias_tokens`, `--fixed_positions` (`"B:12,15 C:3"`, **1-based per chain**), `--placement_metric` (`ddE`), `--coupling_seeds` (4), `--coupling_metric` (`dE_complex`). Tokens split on whitespace. **A token repeating a flag the tool already sets is refused** (argparse takes the last occurrence, which would override a validated value invisibly). `design_ph.py` has no `--checkpoint` — that belongs to the upstream client. |

### Cost

Per structure: a **coupling scan** of `S × C × (N+1)` forward passes (`S` distinct seeds — 1 under
`native`; `C` centres; `N` candidate residues, binder *and* target, within `--coupling-radius` of the
nearest centre); **infill** of about `4F` forwards per lambda at `--apo-weight 0.5` on a graph
replicated `--num-seq` times, with `F ≤ --max-mutations`; and **scoring** of `2 × C` per design. That
is hundreds of forward passes, not a heavy job — roughly **2–2.5 min per task** at the defaults
(8 draws/lambda, coupling, 2 centres), against 70–107 s for a `neighbourhood`/1-centre run.

**Time one structure before batching.** The task wall is 2 h and a chunk that hits it **discards every
structure in the chunk**, including ones already written, because outputs only return when the Modal
call completes. Use `--designs-per-task 1` first, then size with a factor of two in hand.

## `lambda = 0.0` is the control and it is unrecoverable

Sampling is from `R(a) = log p_target(a) - lambda * log p_comparison(a)`. At `lambda = 0` the selective
term is off, so those designs are **what this backbone produces with no pH selection pressure at all**.
Without them, a run full of good `ddE_mean` cannot be told from *"this backbone was simply good"* — and
the control cannot be added later without re-running the sweep on GPU. Upstream only warns; **this tool
makes it a submit-time error.**

## Trust columns — read these before any science column

The sequence of a row comes from `<stem>.fa`; its numbers come from `<stem>.metrics.csv`. Two files —
if their orders diverge, every design silently inherits another's `ddE_mean`. The collector reports how
it paired them.

| Column | Expected | How to read it |
| --- | --- | --- |
| `join_mode` | **`key`** | **Read first.** `key` = matched on an identifier in both files (the design id). `positional` = matched by order, meaning the header format changed — check the next two rows. `failed` never appears on an OK row. |
| `join_key` | `name` | Which metrics column keyed the join, matched **value by value, not by order**. `name,lambda` is the composite fallback. Neither unique ⇒ the parent FAILS; it is never paired by position. |
| `join_checked_fields` | `name,lambda,n_mutations` | What was compared on every pair. The rounded measurement echoes (`dE_complex`, `ddE`) are **deliberately not** compared — the header prints 3 dp and the CSV 4, and the real data contains a tie that would fail a correct pairing. Non-empty with `positional` = order was cross-checked; empty with `positional` = genuinely unverified. |
| `centers_seq_check` | **`True`** | **The one check that catches a wrong index base.** For every centre, `sequence` at that 0-based index must carry the residue the token names. It is **not** corroboration of which residue the centre faces — the model was told to pin that token there. It catches a changed index base upstream, and a designed chain that is not the one you think. `False` never appears on an OK row (it fails the parent). `NA` = could not run. |
| `centers_resnum_key` | **`position_index0`** | **A centre is a 0-BASED INDEX into the FASTA string / graph node order, never a PDB residue number.** To name the residue, count residues along the chain as the structure lists them and take the (index+1)-th. **Do not add a constant offset** unless the chain is gapless and unpadded — a numbering gap becomes an `X` node that still consumes an index. `centers_resnum_key == "residue"` means distrust the middle field entirely (in several schemas that field is the residue *name*). |
| `placement_metric` | `ddE` | **Which quantity ranked the candidate centres.** The same structure gets DIFFERENT centres under `ddE` and `dE_complex`, so split a merged table on this as on `model_noise`. |
| `placement_score_min` | **no fixed bar** | `min(placement_score)` — how the **weakest** centre scored. **Read `placement_metric` first and never compare across metrics.** Nats, **unbounded, may be negative**, and already includes the charge prior, so it is not purely "how much the model wanted it". Take the bar from THIS run's distribution (e.g. a percentile of the `lambda == 0.0` arm) and record it in the campaign log. The MINIMUM, not the mean: one weak centre makes the whole placement doubtful. NA, never 0. |
| `centers_match_metrics` | `True` | **A FILE-PAIRING check, not an independent one** — both sides are written from the same in-memory object, so `True` is near-tautological. It catches a stale `centers.json` left in a reused leaf beside a fresh `metrics.csv`. NA = unreadable, which is blindness, not failure. |
| `model_noise` | `n03_v2` | **Which checkpoint produced these centres.** The noise level decides WHICH CENTRES are picked, so rows from two levels are not comparable and must be split before ranking. |
| `image_commit` | the pinned commit | Which build of a moving private branch produced the row. **Never mix two values in one table** — they are two different programs. |
| `apo_weight` | `0.5` | What the SUBMIT asked for. **Cross-check against `metrics_apo_weight`**, design_ph's own echo: a disagreement means the flag did not reach the program. |
| `lam_sweep` | `0.0 0.25 0.5 0.75` | Whether the `0.0` control arm was in the run at all — a per-row `lambda` cannot tell you that. |
| `binder_chains` / `target_chains` / `center_mode` / `output_scope` | as invoked | Which molecule, and which pH behaviour, these energies describe. |
| `seq_len` | the binder length you expect | The cheapest detector of a swapped chain assignment. |
| `ext_join_mode` / `ext_checked_fields` / `ext_status` | `key` / … / `OK` | How `sequence_ext` was paired. On failure the column is emptied and `ext_status` says why — **the row's own `_status` stays `OK`**, because that is one column's problem. |
| `n_metrics_rows` / `n_fa_entries` / `n_expected` | equal; `len(--lam) × --num-seq` | `rows < expected` ⇒ upstream dropped draws. `rows > expected` ⇒ a preamble entry was collected as a design. |
| `n_preamble_dropped` | **`1`** native, **`1 + --num-seq`** generate | Leading FASTA entries that are not designs, dropped by **exact** match on the first comma-token (never a substring test). A lower number means an entry was not recognised, and the parent will have failed. |
| `output_predates_meta` | `False` | **Advisory.** `True` = some files are older than the leaf's `.meta.json`, so either the leaf was submitted twice and these are the earlier run's results (making the recorded parameters wrong for these rows), or the run_dir was copied without preserving mtimes — a plain `copytree` flips it without a byte changing. Treat as "check this leaf's history", not proof. |
| `centers_source` | `centers` | `unparsed` ⇒ the centre columns are NA and say nothing about the metrics columns, which come from a different file. |
| `entry_index` / `header` / `metrics_path` | — | The `_p<i>` suffix, the raw FASTA header, and the file — for checking a join by hand. |

**An empty string comes back as `NaN`, so never filter one with `!= ""`.** prosapia reads tables
without `keep_default_na=False`, so every column written as `""` is `NaN` by the time a `-f` module
sees it, and `NaN != ""` is **True** — a `!= ""` gate admits every row it meant to exclude. Use
`.notna()`.

**Contradictions fail; they never become a positional guess.** The parent gets a single `<parent>_p0`
row with `status = "error: …"`, the counts, and **no numbers at all** — not a `ddE_mean` of 0 — when:
counts differ; design ids are not unique; the two files name different draws; a cross-checked field
disagrees; a centre does not match the designed sequence; a metrics row has more fields than its
header; a required file is missing; or the stem's files appear in two places.

## Gotchas, each with its signature

- **TWO DIFFERENT UNITS OF "A DISTINCT DESIGN".** `--num-seq` draws are frequently not distinct, and
  which column tells you depends on the question:

  | Question | Dedupe on |
  | --- | --- |
  | how many designs to FOLD | `<leaf>_sequence` |
  | how many designs an ENERGY describes | **`<leaf>_sequence_ext`** |

  Two rows with an identical canonical `sequence` differed in ONE protonation token (`GLU-P` vs
  `GLU-D`) and their energies differed accordingly (0.0198 vs 0.0826 nats). Those are two protonation
  hypotheses on one amino-acid sequence; averaging them, or dropping one as a duplicate, mixes or
  discards a real result. Rows with identical `sequence_ext` always had identical energies (≤ 2.4e-6).
  *Signature:* `nunique()` on `sequence_ext` exceeding `nunique()` on `sequence`.
- **Draws within a lambda arm are often the SAME design, and arms are PAIRED by sample index.**
  Upstream reseeds before each arm, so `_s3` at λ=0 and `_s3` at λ=0.5 saw the same decoding order —
  **compare arms pairwise**, not as independent samples. Unique sequences out of 8 draws:

  | run | λ=0 | λ=0.25 | λ=0.5 | λ=0.75 |
  | --- | --- | --- | --- | --- |
  | one centre | 2–4 | 4–6 | 5–6 | 8 |
  | two centres | 6–7 | 7–8 | 7–8 | 8 |

  Diversity rises with λ and with centre count. **Never report a row count as a design count** — at
  λ=0 with one centre you may be paying for 8 predictions of 2 designs. Raise `--temperature` or widen
  `--lam` to buy real diversity; raising `--num-seq` alone may buy none.
- **Row count is a product.** `len(--lam) × --num-seq` = **32 rows per structure** at the defaults, not
  8. *Signature:* the child table is 4× what you budgeted and `lambda` has four distinct values.
- **`ddE_mean` hides whether both centres switch.** With 2 centres a mean of 0.45 is "two at 0.45" and
  "one at 0.80, one at 0.10" alike, and only the second has a histidine doing nothing. **Rank on
  `ddE_mean`, then gate on `ddE_min`.** *Signature:* `ddE_min` far below `ddE_mean`.
- **`dE_complex` and `dE_apo` go NA the moment there are two centres** — they are `;`-joined
  per-centre strings upstream. *Signature:* a filter on `dE_apo` that worked on a one-centre run
  silently selects nothing. **Filter on `dE_apo_mean`/`dE_complex_mean`**, which are always numeric,
  and read `dE_apo_centres` for the breakdown.
- **`--centers-per-type` is per TOKEN TYPE**, so `2` is 4 centres under `acid_p` and 6 under `all_p`.
  *Signature:* `n_centers` is 4 where you budgeted 2, and `dE_complex`/`dE_apo` go NA.
- **A wrong-but-present chain letter is not caught.** The disjointness and non-empty checks catch
  typos, not a wrong assignment. *Signature:* `seq_len` is the target's length while `binder_chains`
  shows exactly the letter you asked for. Nothing else in the table looks wrong.
- **A structure that fails inside a batch is silent at task level.** `design_ph.py --pdb_dir` does not
  abort the batch on a per-structure failure (usually **a requested chain not in that file**): the task
  exits 0 and that stem has no output. *Signature:* collect prints `N of M parent(s) produced no
  design_ph output` and those parents have **no child rows at all**. Read design_ph's trailing failure
  list at the end of that task's `.out`. **After a filtered run (`-f`) that line is EXPECTED** — collect's
  ready set is the whole parent table, so every excluded design is counted there.
- **Re-running the SAME dir label with different settings is refused at submit.** `grp_<g>` numbering
  restarts at 0, nothing prunes old output, and `.meta.json` is rewritten — so collect would stamp the
  first run's outputs with the second run's parameters. The submit refuses, names the stems and deletes
  nothing. **Use a new `--table-label` for a changed run.**
- **A chunk timeout loses structures that already finished.** *Signature:* a non-zero `.exit` for one
  task and several parents with no rows while their neighbours collected fine.
- **An atomium-habit noise level is refused by design.** *Signature:* `--model-noise: 'n07' is not one
  of n03, n05, n03_v2, n05_v2`. `n00…n07` are the 21-class base weights.
- **`--set` cannot "just override" a flag** in any position — every word is inspected. *Signature:* a
  submit-time `ValueError: --set: --binder_chains is set by this tool…`. A `*`, `?` or `[` in a value
  is safe (`set -f` is applied around the call).
- **`sequence_ext` is not a sequence a predictor can read** — it carries protonation tokens. Feeding it
  to Boltz/AF3 gives garbage or a refusal.
- **`summary.csv` is ignored on purpose** — it ranks designs across structures, which the table does
  natively with lineage intact.
- **The `_p<i>` row-name suffix is shared with `protonpottsmpnn`.** The leaf prefix is the only thing
  telling two child tables apart.

## What it collects

All columns leaf-prefixed `atomiumph_`. **Read the trust columns above first.**

| Column | Meaning | How to read it |
| --- | --- | --- |
| `status` | `OK` or `error: …` | **The only proof of success.** An error row is the single `<parent>_p0` row, with counts and **no science columns**. |
| `sequence` | Canonical 1-letter binder sequence | **What folding tools consume** (`-i atomiumph_sequence`, via `mkcomplex`). |
| `sequence_ext` | The same design with explicit protonation tokens | For reading and an orthogonal pKa check. **Never a predictor input.** Validated per row; emptied on failure with `ext_status` saying why. |
| `lambda` | The sweep value this draw was sampled at | `0.0` rows are the control arm. |
| `ddE_mean` | `mean(dE_complex) − mean(dE_apo)` over centres | **The selection metric. Higher is better.** Equal improvement in complex and apo is a pKa shift, not a switch — which is why `dE` alone must not be ranked on. **Read it WITH `ddE_min`.** |
| `ddE_min` | **Minimum** per-centre ddE | **Whether BOTH centres carry the switch or one carries all of it** — the question `ddE_mean` cannot answer. NA when an apo value is `nan` (no target chains). Equals `ddE_mean` with one centre. |
| `ddE_centres` | `;`-joined per-centre ddE | In the same order as `centers` — **by token type then descending placement score**, not by position. |
| `dE_complex` | `log p(target) − log p(comparison)` in the complex | Higher = the pinned state is more favoured with the target present. **Not the selection metric alone. Numeric only with ONE centre; NA with several.** |
| `dE_apo` | The same on the apo binder | Strongly **positive** means the free binder's centre will not change state in the accessible pH range — nothing switches. **NA with several centres; filter on `dE_apo_mean`.** |
| `dE_complex_mean` / `dE_apo_mean` | Full-precision means | Always numeric. Upstream's selection recipe names `dE_apo_mean`. |
| `dE_complex_centres` / `dE_apo_centres` | The verbatim `;`-joined cells | `"nan;nan"` on the apo side means there were no target chains, and nothing in the row is a switch measurement. |
| `interface_fraction` | Share of the REDESIGN SET within `--contact_cutoff` of the target | **READ IT BEFORE BELIEVING A FLAT ddE CURVE.** With the redesign set inside the binder rather than at the interface, `dE` improves equally in complex and apo and **`ddE` cannot move by construction** — λ appears to do nothing and the null result is an artefact of the objective. **Below 0.5, a flat ddE-vs-λ curve means nothing.** `NaN`, not 0, without target chains. |
| `n_mutations` | Mutations vs the starting sequence, **as that commit defined it** | Bounded by `--max-mutations`; `0` means the seed came back unchanged. **The definition changed between commits under an unchanged name — never compare across `image_commit`.** Use `n_diff_native`. |
| `n_diff_native` | Hamming distance `sequence` vs `native_sequence` | **What makes `n_mutations` checkable**, computed from the sequences so it cannot inherit a header's bookkeeping. NA when the echo is absent or lengths differ. |
| `native_sequence` | The input structure's own binder sequence | Identical on every child row of one parent. **Under `--seed-mode native` it is also the seed.** |
| `n_his` / `n_his_native` | Histidines in the designed and input binder | Under `his_p`/`his_s` with 1 centre the design is normally `n_his_native + 1`; `n_his <= n_his_native` means it traded a His rather than adding a site. |
| `centers` | `<index0>:<token>;…` | The audit of what was pinned; repeats across all children of one parent. **A 0-BASED INDEX into the designed sequence, not a residue number** — read `centers_resnum_key` and `centers_seq_check` first. No chain field exists in the file. |
| `n_centers` / `n_redesign` | Counts from `centers.json` | `n_redesign` = `len(free_positions)`, matching `Stage 2: redesign positions = N` in the `.out`. NA, never 0, when absent. |
| `charge_field_min` / `charge_field_max` | **ALWAYS NA** | The range is printed in the task `.out` only. Kept so another build that records it needs no code change. |
| `metrics_name` / `metrics_centers` / `metrics_apo_weight` | Upstream's own cells | Renamed so they cannot collide with the table's `name` index or a trust column. **Compare `metrics_apo_weight` against `apo_weight`.** |
| any other metrics column | Passed through verbatim | Upstream can add a column and it appears with no code change. |
| `couplings_path` | Path to `<stem>.couplings.csv`, else `""` | **A file, not columns** — a per-position influence map is not one value per design. Columns: `position, chain, is_binder, is_designable, residue, dist_to_center, dE_shift_<pos>_<TOKEN>, max_abs_shift{,_complex,_apo}, in_redesign_set`. Written only under `--redesign-scope coupling`. |
| `atomiumph_path` | The parent's canonical FASTA | Provenance only. **NEVER USE IT AS AN INPUT COLUMN — see below.** |

On disk per stem (globbed; the observed layout is `<out_dir>/grp_<g>/out/`): `<stem>.fa`,
`<stem>.ext.fa`, `<stem>.metrics.csv`, `<stem>.centers.json`, `<stem>.couplings.csv` (coupling scope
only), plus upstream's ignored `summary.csv`.

### Never pass `-i atomiumph_path` downstream

It is the **parent's whole `<stem>.fa`** — every design on that parent *plus the `_native` echo* — and
it is the **same path on every child row**. A downstream tool would read N sequences where it expected
one, silently include the undesigned input, and repeat the work once per child row.

**The per-design sequence column is `atomiumph_sequence`**: `sapia run mkcomplex <run_dir> -t <child>
-i atomiumph_sequence …`. The same applies to `atomiumph_metrics_path`.

## Filters

### Two helpers every `-f` module here should start with

Derive the prefix (never hard-code `atomiumph_`), and never read a bool-with-NA column with
`.astype(bool)`:

```python
def _leaf(df):
    # "atomiumph", or "atomiumph_<dir_label>" if the run used -l. Requiring the
    # _status/_path pair excludes the tool's own _ext_status (which has no _ext_path).
    hits = [
        c[: -len("_status")]
        for c in df.columns
        if c.startswith("atomiumph")
        and c.endswith("_status")
        and c[: -len("_status")] + "_path" in df.columns
    ]
    if len(hits) != 1:
        raise KeyError(f"expected exactly one atomiumph leaf, found {hits or 'none'}")
    return hits[0]


def _flag(df, col, *, na):
    # A bool-with-NA column as a real mask. After the TSV round trip it can be object
    # dtype holding the STRINGS "True"/"False", and bool("False") is True.
    raw = df[col]
    text = raw.astype(str).str.strip().str.lower()
    return text.eq("true").where(~raw.isna() & text.ne("nan"), na)
```

`pd.Series(["True", "False"]).astype(bool)` is `[True, True]`, so `.fillna(False).astype(bool)` can
admit **every row it was written to exclude**. The `na` argument differs **on purpose**:

| Column | NA means | `na=` |
| --- | --- | --- |
| `centers_seq_check` | the index mapping was never checked | **`False`** — unverified is not verified |
| `centers_match_metrics` | a separator this tool cannot parse | **`True`** — our blindness is not the row's fault |
| `output_predates_meta` | no sidecar to compare against | **`False`** — no evidence of staleness |

### Upstream's selection recipe

```python
def apply_filter(df):
    L = _leaf(df)
    ok = df[f"{L}_status"] == "OK"
    # Status OK only says the files parsed, not that the sequence and the numbers
    # belong to each other. .notna(), NOT != "": an empty cell round-trips as NaN.
    mapped = (df[f"{L}_join_mode"] == "key") | df[f"{L}_join_checked_fields"].notna()
    centres_ok = _flag(df, f"{L}_centers_seq_check", na=False)
    switch = df[f"{L}_ddE_mean"] > 0.5            # campaign-specific; see below
    apo = df[f"{L}_dE_apo_mean"] < 1.0            # _mean: dE_apo is NA with 2 centres
    # With >1 centre, require the WEAKEST centre to switch too, not just the average.
    every_centre = df[f"{L}_ddE_min"].fillna(df[f"{L}_ddE_mean"]) > 0.2
    return df[ok & mapped & centres_ok & switch & apo & every_centre].sort_values(
        f"{L}_ddE_mean", ascending=False)
```

**The `ddE_mean` bar is a property of the campaign, not of the table.** There is no published
threshold. Set it from the `lambda == 0.0` control of *these* backbones — e.g. the 90th percentile of
the control's `ddE_mean` — and record the number in the campaign log.

### The full trust gate

```python
def apply_filter(df):
    L = _leaf(df)
    ok = df[f"{L}_status"] == "OK"
    mapped = (df[f"{L}_join_mode"] == "key") | df[f"{L}_join_checked_fields"].notna()
    centres_ok = _flag(df, f"{L}_centers_seq_check", na=False)
    # PLACEHOLDER: no portable number exists. Read <leaf>_placement_metric, then take
    # the bar from THIS run's distribution.
    MIN_PLACEMENT = ...        # e.g. df[f"{L}_placement_score_min"].quantile(0.25)
    wanted = df[f"{L}_placement_score_min"].fillna(float("-inf")) >= MIN_PLACEMENT
    centres_agree = _flag(df, f"{L}_centers_match_metrics", na=True)
    fresh = ~_flag(df, f"{L}_output_predates_meta", na=False)
    one_model = df[f"{L}_model_noise"] == "n03_v2"        # never mix noise levels
    one_build = df[f"{L}_image_commit"] == df[f"{L}_image_commit"].mode()[0]
    one_placement = df[f"{L}_placement_metric"].fillna("") == "ddE"
    right_molecule = df[f"{L}_seq_len"].between(60, 120)  # your binder's length
    return df[ok & mapped & centres_ok & centres_agree & fresh & wanted
              & one_model & one_build & one_placement & right_molecule]
```

### De-duplicate on the right key

Before FOLDING, dedupe on `sequence` — a predictor cannot tell two protonation hypotheses apart:

```python
def apply_filter(df):
    L = _leaf(df)
    ok = df[df[f"{L}_status"] == "OK"]
    return (ok.sort_values(f"{L}_ddE_mean", ascending=False)
              .drop_duplicates(subset=[f"{L}_sequence"], keep="first"))
```

Before counting or averaging ENERGIES, dedupe on `sequence_ext` instead — same `sequence`, different
tokens is **two results, not a duplicate**:

```python
def apply_filter(df):
    L = _leaf(df)
    ok = df[df[f"{L}_status"] == "OK"]
    return ok.drop_duplicates(subset=[f"{L}_sequence_ext"], keep="first")
```

## Invariants to check on the first real output

Cheap, loud on failure, and none is checked by the collector — they span the table, the task `.out` and
the input structure, which no single row can see.

1. **The binder is the molecule you named.** These must be equal: `len(native_sequence)`; the `N` in
   `Native binder: N residues` in the `.out`; the residue count of the `--binder-chains` chain(s) in
   the input; and `seq_len`. **If they disagree the run redesigned something other than the binder.**
2. **Under `--seed-mode native`, `n_diff_native == n_mutations` on every row.** A difference means a
   position outside the redesign set moved, or `native_sequence` is not the seed. (Under `generate`
   they are expected to differ.)
3. **Every design id distinct**, and `n_fa_entries == n_metrics_rows == n_expected`. Check
   `n_preamble_dropped`: **1** in `native`, **1 + --num-seq** in `generate`.
4. **`apo_weight == metrics_apo_weight`** on every row — the submit record against design_ph's echo.
5. **`len(sequence_ext.split()) == seq_len`**, i.e. `ext_status` is `OK` everywhere.
6. **The dtypes of the three bool-with-NA columns after the round trip.** If they came back as object
   dtype holding `"True"`/`"False"`, any `.astype(bool)` filter admits everything — use `_flag`.

Read once from the first task's `.out`: the `design_ph.py --help` block, the resolved checkpoint path
and size, the `Charge prior … field range […]` line and the per-centre `placement_score` lines — none
of which reach a column.

## What it is blind to

- **Whether the sequence folds, or folds to the designed backbone.** Nothing here is a structure
  prediction. Fold `sequence` (`mkcomplex` → `boltz`), test self-consistency with `usalign`, read
  confidence with `boltzconf`.
- **Whether the binder still binds at all.** `ddE` is a *difference* — a design can improve its switch
  while destroying the interface. Re-score after folding with `cms`, `pyrosetta` (`if_dG`, `if_dSASA`),
  `epitope` (did it stay on the epitope) and `dockq` (did it keep the pose).
- **The actual pKa, and so the actual pH window.** Model log-probabilities, not a titration. Confirm
  with an orthogonal pKa calculation before believing a pH number.
- **Whether the chains you named are the chains you meant.** `chainsel`/`epitope` establish chain
  identity.
- **Whether a FLAT ddE-vs-λ curve is a result or an artefact.** `interface_fraction` answers it; below
  0.5 a flat curve means nothing. Raise `--apo-weight` or tighten `--interface-cutoff` and re-run.
- **Whether the centre sits in the RIGHT-SIGN charge environment.** `placement_score_min` is the
  **sum** of the model and charge terms, so a strong model term hides a wrong-sign prior. Only the
  `.out` decomposes them.
- **WHERE the centre is in the target's frame.** `centers` is an index into the binder and says nothing
  about which target residue it faces — that is `epitope`'s job.
- **Nothing upstream records the run parameters.** Noise, seed, `--num-seq`, the lambda list, chains,
  `--centers-per-type`, `--charge-prior` and `--protect-titratable` appear in **no** output file; they
  are in the table only because this tool writes them to `.meta.json`. A leaf whose sidecar and
  `atomiumph_tasks/` are both gone cannot be reconstructed.
- **Expression, aggregation and developability.** `--ban-cys 1` removes one failure mode; it is not a
  developability screen.
