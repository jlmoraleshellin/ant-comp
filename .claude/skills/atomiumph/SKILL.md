---
name: atomiumph
description: >-
  Which sequences make a pinned His/Asp/Glu protonation state favoured in the complex more than in the free binder (a pH-selective binder)? Wraps AtomiUM's `design_ph.py`: centre modes, required chain flags, the lambda sweep, `--binder-alone`, trust columns, `dE`/`ddE` energies and the numeric `centre<k>_idx1` positions a downstream flag pins through `{expr}`. Load before composing an atomiumph run, interpreting its numbers, or handing its centres to another tool. Traps: `--center-mode` has no default and getting it backwards yields the opposite molecule; a wrong binder letter still runs; `--output-scope all` suppresses `centre<k>_idx1` rather than publishing a position on the wrong molecule.
---

# atomiumph

**Custom tool** (`tools/atomiumph/`), **`action: create`**: mints a child table, **one row per designed sequence** (`<parent>_p0`, `<parent>_p1`, …). Runs on both executors: on vib through `SAPIA_ACTIVATE_ATOMIUMPH`, which points at the lab's shared `atomium.sh` activation script (**H100-only** there: the shared environment is PyTorch 2.1.2 / CUDA 12.3), and on Modal through its own image.

Wraps `design_ph.py` from the private AtomiUM repo (branch `pure_wo_jit`), the **sibling** of `tools/atomium`, not a mode of it: same repo, different entry point, **its own commit pin and image**. `tools/atomium` designs sequences for a backbone and reports sequence recovery (no protonation states, no apo reference, no `dE`/`ddE`). `tools/protonpottsmpnn` is a different program reporting Potts energies with **no complex-vs-apo contrast**, which is the quantity `ddE` isolates.

## Premise, and when its numbers are meaningless

`design_ph.py` pins titratable binder residues (His/Asp/Glu) into a chosen protonation state and redesigns their environment so that state is favoured **in the complex more than in the free binder**:

```
dE  = log p(target_state) - log p(comparison_state)      one centre, one context
ddE = dE(complex) - dE(apo binder)                       the selection metric
```

**The input must be a complex.** The apo reference is the same file with the target chains deleted.

| Situation | What happens |
| --- | --- |
| No target chain in the input | Refused at submit (`--target-chains` is required) **unless `--binder-alone`** |
| `--target-chains` names the wrong chain | **It still runs.** `ddE` is computed against the wrong apo reference and reads like a result; nothing in the output detects it |
| `--binder-chains` names the target | **It still runs and redesigns the target.** `seq_len` against the expected binder length is the only cheap tell |

**There is no fixed binder chain letter.** The house convention is that the binder is the LAST chain, so the letter depends on the target (two-chain target → `C`, one-chain → `B`, trimer → `D`); `design_ph.py`'s own default `--binder_chains A` is wrong here. Establish it from `chainsel`'s `chain_n_res`, an `epitope` chain list, or the **pre-flight** (below). Both chain flags are required and must be disjoint; the tool asserts nothing about which is which.

**Energies are model log-probabilities on one backbone.** Comparable between draws of the same parent and across the lambda sweep; **not** a calibrated pKa, and `ddE_mean` across unrelated backbones is weakly meaningful at best: rank within a parent first.

## `--center-mode`: which pH behaviour you get

Getting this backwards produces a healthy-looking table of the opposite molecule. **Required, no default.**

| `--center-mode` | Pins | Binds | Releases | Charge prior seeks |
| --- | --- | --- | --- | --- |
| `his_s` | neutral His (0) | **neutral** (plasma 7.4) | acid (endosome 6.0–6.5) | **positive** (Arg/Lys) |
| `his_p` | protonated His (+1) | **acid** (tumour 6.0) | neutral (7.4) | **negative** (Asp/Glu) |
| `acid_p` | protonated Asp/Glu | acid | neutral | negative |
| `all_p` | every protonated type | acid | neutral | negative |

- **Recycling / "sweeping" binder** (bind in plasma, release the antigen in the endosome): **`his_s`**. **Tumour-microenvironment binder** (bind acidic, release at neutral): **`his_p`**.
- `acid_p`/`all_p` pin Asp/Glu (pKa ~4, switch only below ~5.8): **diagnostic only, unusable for either His challenge**, which is why upstream's `acid_p` default is not inherited.
- **The upstream README's section-1 use-case table has Binds/Releases SWAPPED**; its section-4 `--center_mode` table, section-5 examples and the client README say the opposite. The table above is right: do not re-derive the direction from the README.
- **The `his_s` charge prior is sign-reversed.** The prior adds `-dq * S` with `S` the signed sum of nearby charged **target** residues; `his_p`/`acid_p` have `dq = +1` and seek a negative environment, but **`his_s` has `dq = -1` and seeks a POSITIVE (Arg/Lys) environment**, so "put the histidine next to an aspartate" is exactly backwards for the recycling direction. Target histidines are skipped by the prior.

**Pre-flight tells you whether the epitope can work at all.** `--lam 0.0 --num-seq 2` prints into the task `.out`, before generating anything: the chains found, the banned tokens, **the charge-field range** (`Charge prior on (weight=0.5); field range [-1.51, +1.94]`) and the centres with their scores (`pos 10 HIS-S placement_score +0.708 model +0.657 charge +0.10`). **The charge-field range and the model/charge decomposition reach no column**; the `.out` is the only place. If the field has the wrong sign for the mode, that epitope lacks the complementary charges and no sampling fixes it.

## `--binder-alone`: the binder's OWN core, and no ddE at all

`design_ph.py` defines `core` as a **Cb-neighbour count within 10 Å over every chain in the graph**, masked to the binder only afterwards (`neighbours = (cdist(cb, cb) < 10).sum(dim=1) - 1`, `core = (neighbours > burial_cutoff) & is_design`). So **`--placement-region core` on a complex means "buried in the complex"**: a centre buried under the *target* qualifies while fully exposed in the free binder. *Signature:* free-binder burial near 0% for the chosen centres while the run reported them as core.

`--binder-alone` (opt-in, default off) puts **no target chain in the graph**, so `core` is the binder's own core. It is a **different question, not a cheaper run**: where a centre sits inside the binder and what it costs in `dE_complex`; it says **nothing** about a pH-dependent switch.

| In `--binder-alone` mode | |
| --- | --- |
| `--target-chains` | must be **empty**; passing both is refused |
| `--lam` | only `0.0` accepted (lambda scales the apo-vs-complex term; other arms would be duplicates) |
| `atomiumph_placement_metric` | reads **`dE_complex`**, the metric that actually ranked the centres (without an apo graph upstream falls back whatever was requested; the collector records the effective value) |
| `--charge-prior`, `--apo-weight` | **forced to `0.0`** and recorded as `0.0` in `atomiumph_charge_prior` / `atomiumph_apo_weight` (the prior sums TARGET charges only; `apo_weight` is zeroed rather than trusted to be ignored) |
| `--placement-region interface` | **refused** (the mask is all False without a target) |
| `--set '--placement_metric ddE'` / `'--coupling_metric ddE'` | **refused** (both fall back to `dE_complex`; `--coupling_metric` reaches no column) |
| input structure | **checked per file at submit**: any chain other than `--binder-chains` (ATOM or HETATM) refuses the submit, naming file and chain (a stray chain could add padding rows at the origin that the unmasked `cdist` counts as neighbours, inflating burial invisibly). Cut the binder out with `chainsel` and point `-i` at that column |
| `--table-label` | **its own, always**: collect stamps a whole leaf from one `.meta.json`, so a shared leaf would mislabel every row of the other mode; flipping the flag in a leaf that holds output is refused |
| `atomiumph_binder_alone` | **`True` on every row** |
| NA on every row, never 0 | `ddE_mean`, `ddE_min`, `ddE_centres`, `dE_apo`, `dE_apo_mean`, `dE_apo_centres`, `interface_fraction`, `charge_field_min`, `charge_field_max` |
| still real | `dE_complex`, `dE_complex_mean`, `placement_score_min`, `n_centers`, `centers`, `centers_seq_check`, `n_his`, `n_mutations`, `n_diff_native`, `seq_len` |
| not blanked, not a column | `couplings_path`: the CSV's apo columns are **structural zeros** here |
| row status | `OK`: a valid run of a different question |

**`dE_complex` here is the FREE-BINDER energy** (the graph built is the binder alone): **not comparable** with a `dE_complex` from a run that had a target, nor is `placement_score_min`. **The NA list is an allowlist of column names**: a target-dependent column a future design_ph adds passes through under its own name, possibly as a structural zero; re-read the metrics header on a new pin. **At `lambda 0.0` with one centre the draws repeat** (see the draws gotcha): set `--num-seq` high.

## Fork with `--table-label`, never with `-l`

`-l/--dir-label` becomes part of the **column leaf**, renaming every column the tool writes (`atomiumph_status` → `atomiumph_nv2_n03_status`), so every filter written as `df["atomiumph_status"]` raises `KeyError`, and **`collect` then needs the SAME `-l`** or it looks in `…/atomiumph/`, finds nothing and reports every parent as having produced no output. `--table-label` separates both the table and the leaf directory while keeping plain `atomiumph_` columns. Use `-l` only for **two leaves feeding ONE table**, where the rename keeps their columns apart.

## Invocation

```bash
# 1. pre-flight ONE structure: prints chains, banned tokens, charge-field range and centres before generating
sapia run atomiumph <run_dir> -t table0 --table-label preflight \
    --binder-chains <binder> --target-chains <targets> \
    --center-mode his_s --lam 0.0 --num-seq 2 -f pick_one_design.py

# 2. the sweep: only the three flags below are required; the rest is the recipe
sapia run atomiumph <run_dir> -t table0 --table-label ph \
    --binder-chains <binder> --target-chains <targets> \
    --center-mode his_s          # his_p for a tumour/acid-binding design

# 3. --binder-alone: its own --table-label, an input column holding the binder chain ALONE (chainsel output), no target
sapia run atomiumph <run_dir> -t table0 --table-label corealone -i chainsel_path \
    --binder-alone --binder-chains D --center-mode his_s --placement-region core \
    --lam 0.0 --centers-per-type 1 --max-mutations 5

sapia collect atomiumph <run_dir> -t <the table the run reserved>   # ...and the SAME -l as the run, if it used one
```

`default_input_column` is **`rfdiffusion3_path`** (rfd3 binder output is already a complex). From any other parent pass `-i` and check the file still holds **both** chain sets: `-i boltz_path` on a binder-only prediction, or a `chainsel` output that extracted the binder, has no target and the run is meaningless though it submits.

Resources: **L4, 8 CPU, 16 GiB, 2 h** (Modal). **(Modal)** the image builds from the private repo using the `github-token`/`github-username` Secrets; the first submit pays a long build (torch + cu128 + weights), allow 15+ minutes.

## Flags

Defaults are the upstream README v2 production recipe, not `design_ph.py`'s own argparse defaults; the last column shows what applies if a flag were omitted upstream (client / script, only where they differ from ours).

| Flag | ours | upstream (client / script) | |
| --- | --- | --- | --- |
| `--binder-chains` | **required** | `A` / required | the binder is the LAST chain here |
| `--target-chains` | **required** | `""` / `""` | empty target ⇒ every ddE structurally zero; empty allowed ONLY with `--binder-alone` |
| `--binder-alone` | **off** | n/a | no target in the graph; `core` becomes the binder's own core, every target-dependent column NA |
| `--center-mode` | **required** | `acid_p` / `acid_p` | a default would answer a question nobody asked |
| `--centers-per-type` | **`2`** | `1` / `1` | README v2 |
| `--redesign-scope` | **`coupling`** | `neighbourhood` / `neighbourhood` | README v2 "the default choice for production" |
| `--interface-cutoff` / `--coupling-radius` | **`8.0`** / **`12.0`** | `12.5` / `15.0` for both | README v2's strict first shell |
| `--model-noise` (alias `--noise`) | **`n03_v2`** | `n05` / `n05` | chosen default, not upstream's |
| `--apo-weight` / `--protect-titratable` | `0.5` / `1` | client `0.5` / `1`, script `0.0` / `0` | client ≠ script, so both are always passed |
| `--lam` | `0.0 0.25 0.5 0.75` | same | **must contain `0.0`** |
| `--num-seq` | `8` | same | **PER LAMBDA** |
| `--temperature` / `--max-mutations` | `0.1` / `20` | same | README's examples use 0.15 / 12 |
| `--seed-mode` / `--seed-temperature` / `--seed` | `native` / `0.1` / `42` | same | |
| `--charge-prior` / `--charge-cutoff` | `0.5` / `8.0` | same | sign follows `--center-mode` |
| `--placement-region` / `--output-scope` | `interface` / `binder` | same | |
| `--neighbour-k` / `--burial-cutoff` | `12` / `16` | same | |
| `--coupling-probe` / `--ban-cys` | `ala` / `1` | same | |
| `--no-control`, `--set`, `--designs-per-task` | none, none, **`5`** | none | ours only; `--designs-per-task` = structures per task (atomium's is 10) |

| Flag | Note |
| --- | --- |
| `--centers-per-type` | Centres per **TOKEN TYPE**: `his_p`/`his_s` have one type (`2` → **2 centres**), `acid_p` two (Asp, Glu → **4**), `all_p` three (→ **6**). Why 2: one histidine gives a transition too shallow to separate pH 6.5 from 7.4 when the window is ≤ 1.5 units. Centres scale the `coupling` scan linearly and each adds a value to the `;`-joined per-centre cells. **Only the first 12 centres get a `centre<k>_idx1` column**; past that `n_centres_idx1 < n_centers` and `centres_idx1_basis` reads `truncated`, so the truncation is never silent |
| `--interface-cutoff` | `interface` is a pure Cα-distance mask with **no burial filter**, and nothing excludes the N-terminal residue (a centre at index 0 is normal). Upstream's 12.5 admits a second shell, so an "interface" centre need not touch the target |
| `--model-noise` | `n03`, `n05`, `n03_v2`, `n05_v2` → `ft_weights_for_pH/ft_model_epoch10_<noise>.pt`, fine-tuned **30-class** pH checkpoints. `tools/atomium`'s `n00…n07` are the **21-class base** weights and are refused at submit. **`_v2` are DIFFERENT weights, not newer builds** |
| `--lam` | Selective-term weight sweep, **per value**; duplicates and values colliding at 6 dp are refused; λ ≥ 1.0 is degenerate (uniform sampling) |
| `--no-control` | Opt-out of the `0.0` requirement; honest only when that control exists for the same parents in another table |
| `--apo-weight` | Weight of the apo term **in the sampling objective**, not only scoring: above 0 it pushes toward *selectivity* rather than merely good in the complex, and roughly doubles the infill passes |
| `--protect-titratable` | `1` removes every **natively** His/Asp/Glu binder residue from the redesign pool in all scopes; centres are chosen *before* this, so a protected residue can still be a centre |
| `--seed-mode` | `native` redesigns the sequence already in the PDB; **`generate` is needed for a raw poly-Gly backbone** (`native` would seed from nothing) |
| `--redesign-scope` | `coupling` ranks candidates by **measured** influence on the centres' `dE` and is the only scope writing `couplings.csv`; `neighbourhood` is the cheap pass (graph neighbours, capped by `--neighbour-k`/`--max-mutations`); `binder` is every designable position, for a raw backbone or untrusted sequence |
| `--seed` | Every structure is seeded identically and the pin **reseeds before every lambda arm**, so draw *i* sees the same decoding order at every lambda. `0` is not special |
| `--output-scope` | `all` joins every chain with `/`, making `seq_len` meaningless, disabling the centre sequence check **and suppressing `centre<k>_idx1` entirely** (`centres_idx1_basis == "output_scope_all"`): under `all` the designed sequence is the whole complex, so an index into it is not a binder position and the tool refuses to publish one. Leave at `binder` |
| `--set` | Raw flag forwarded verbatim, repeatable, **for unmodelled flags only**: `--contact_cutoff` (8.0), `--bias_tokens`, `--fixed_positions` (`"B:12,15 C:3"`, **1-based per chain**), `--placement_metric` (`ddE`), `--coupling_seeds` (4), `--coupling_metric` (`dE_complex`). Tokens split on whitespace. **A token repeating a flag the tool already sets is refused** (argparse takes the last occurrence and would override a validated value invisibly). `design_ph.py` has no `--checkpoint` |

Cost per structure: a **coupling scan** of `S × C × (N+1)` forward passes (`S` distinct seeds, 1 under `native`; `C` centres; `N` candidate residues, binder *and* target, within `--coupling-radius` of the nearest centre), **infill** of about `4F` forwards per lambda at `--apo-weight 0.5` (`F ≤ --max-mutations`) and **scoring** of `2 × C` per design: typically 2–2.5 min per task at the defaults, against 70–107 s for a `neighbourhood`/1-centre run. **(Modal)** a chunk that hits the 2 h task wall **discards every structure in the chunk**, finished ones included (outputs return only when the call completes; *signature:* a non-zero `.exit` for one task and several parents with no rows): time one structure (`--designs-per-task 1`) before batching.

## `lambda = 0.0` is the control and it is unrecoverable

Sampling is from `R(a) = log p_target(a) - lambda * log p_comparison(a)`. At `lambda = 0` the selective term is off, so those designs are **what this backbone produces with no pH selection pressure**. Without them a run full of good `ddE_mean` cannot be told from "this backbone was simply good", and the control cannot be added later without re-running the sweep on GPU. Upstream only warns; **this tool makes it a submit-time error.**

## Trust columns: read these before any science column

The sequence of a row comes from `<stem>.fa`; its numbers from `<stem>.metrics.csv`. If their orders diverge, every design silently inherits another's `ddE_mean`; the collector reports how it paired them.

| Column | Expected | How to read it |
| --- | --- | --- |
| `join_mode` | **`key`** | **Read first.** `key` = matched on an identifier in both files (the design id); `positional` = matched by order (the header format changed: check the next two rows); `failed` never appears on an OK row |
| `join_key` | `name` | Which metrics column keyed the join, matched **value by value, not by order**; `name,lambda` is the composite fallback; neither unique ⇒ the parent FAILS, never paired by position |
| `join_checked_fields` | `name,lambda,n_mutations` | What was compared on every pair; the rounded echoes (`dE_complex`, `ddE`) are deliberately not (header 3 dp, CSV 4 dp). Non-empty with `positional` = order cross-checked; empty with `positional` = unverified |
| `centers_seq_check` | **`True`** | **The one check that catches a wrong index base**: `sequence` at each centre's 0-based index must carry the residue the token names. Not corroboration of which residue the centre faces (the model was told to pin it); it catches a changed upstream index base and a wrong designed chain. `False` never appears on an OK row; `NA` = could not run |
| `centers_resnum_key` | **`position_index0`** | **A centre is a 0-BASED INDEX into the FASTA string / graph node order, never a PDB residue number.** To name the residue, count along the chain as the structure lists it and take the (index+1)-th; **do not add a constant offset** unless the chain is gapless and unpadded (a numbering gap becomes an `X` node that still consumes an index). `"residue"` means distrust the middle field entirely |
| `centres_idx1_basis` | **`index0+1`** | **Read before any `centre<k>_idx1`.** Says whether the numeric centre positions were written at all, and if not, which condition refused them (full list in the column table below). The ones to expect in the wild besides `index0+1` are `output_scope_all` and `multi_binder_chain`. |
| `n_centres_idx1` | **= `n_centers`** | How many centres got a position column; anything smaller means some centre has none (cap, or a refused basis) |
| `placement_metric` | `ddE` | **Which quantity ranked the candidate centres.** The same structure gets DIFFERENT centres under `ddE` and `dE_complex`: split a merged table on this as on `model_noise` |
| `placement_score_min` | **no fixed bar** | `min(placement_score)`, the **weakest** centre. **Read `placement_metric` first; never compare across metrics.** Nats, unbounded, may be negative, includes the charge prior. Take the bar from THIS run's distribution (e.g. a percentile of the `lambda == 0.0` arm). Minimum, not mean: one weak centre makes the placement doubtful. NA, never 0 |
| `centers_match_metrics` | `True` | **A FILE-PAIRING check, near-tautological** (both sides written from one in-memory object); it catches a stale `centers.json` left in a reused leaf beside a fresh `metrics.csv`. NA = unreadable (blindness, not failure) |
| `model_noise` | `n03_v2` | **Which checkpoint produced these centres.** The noise level decides WHICH CENTRES are picked: rows from two levels are not comparable, split before ranking |
| `image_commit` | the pinned commit | Which build of a moving private branch produced the row. **Never mix two values in one table** (two different programs) |
| `apo_weight` | `0.5` | What the SUBMIT asked for; **cross-check against `metrics_apo_weight`**, design_ph's own echo (disagreement = the flag did not reach the program) |
| `lam_sweep` | `0.0 0.25 0.5 0.75` | Whether the `0.0` control arm was in the run at all |
| `binder_chains` / `target_chains` / `center_mode` / `output_scope` | as invoked | Which molecule and which pH behaviour these energies describe |
| `binder_alone` | **`False`** | **Read before any ddE column.** `True` = no target in the graph: `ddE_mean`/`ddE_min`/`ddE_centres`/`dE_apo*`/`interface_fraction`/`charge_field_*` are **NA by construction, not missing data**, and `placement_metric` reads `dE_complex` |
| `seq_len` | the binder length you expect | The cheapest detector of a swapped chain assignment |
| `ext_join_mode` / `ext_checked_fields` / `ext_status` | `key` / … / `OK` | How `sequence_ext` was paired; on failure the column is emptied and `ext_status` says why while the row's own `_status` stays `OK` |
| `n_metrics_rows` / `n_fa_entries` / `n_expected` | equal; `len(--lam) × --num-seq` | `rows < expected` ⇒ upstream dropped draws; `rows > expected` ⇒ a preamble entry was collected as a design |
| `n_preamble_dropped` | **`1`** native, **`1 + --num-seq`** generate | Leading FASTA entries that are not designs, dropped by exact match on the first comma-token; a lower number means an entry was not recognised and the parent will have failed |
| `output_predates_meta` | `False` | **Advisory.** `True` = some files are older than the leaf's `.meta.json`: the leaf was submitted twice (recorded parameters wrong for these rows) or the run_dir was copied without preserving mtimes; check the leaf's history |
| `centers_source` | `centers` | `unparsed` ⇒ centre columns are NA and say nothing about the metrics columns (a different file) |
| `entry_index` / `header` / `metrics_path` | | The `_p<i>` suffix, the raw FASTA header and the file, for checking a join by hand |

**An empty string comes back as `NaN`, so never filter one with `!= ""`**: prosapia reads tables without `keep_default_na=False`, and `NaN != ""` is **True**, so a `!= ""` gate admits every row it meant to exclude. Use `.notna()`.

**Contradictions fail; they never become a positional guess.** The parent gets a single `<parent>_p0` row with `status = "error: …"`, the counts and **no numbers at all** (not a `ddE_mean` of 0) when: counts differ; design ids are not unique; the two files name different draws; a cross-checked field disagrees; a centre does not match the designed sequence; a metrics row has more fields than its header; a required file is missing; or the stem's files appear in two places.

## Gotchas, each with its signature

- **TWO DIFFERENT UNITS OF "A DISTINCT DESIGN".** Dedupe on `<leaf>_sequence` to count designs to FOLD, on **`<leaf>_sequence_ext`** to count designs an ENERGY describes. Two rows with an identical canonical `sequence` can differ in one protonation token (`GLU-P` vs `GLU-D`) with different energies: two protonation hypotheses on one sequence, which averaging or dropping as a duplicate mixes or discards. Rows with identical `sequence_ext` always had identical energies. *Signature:* `nunique()` of `sequence_ext` exceeding that of `sequence`.
- **Draws within a lambda arm are often the SAME design, and arms are PAIRED by sample index** (upstream reseeds before each arm, so `_s3` at λ=0 and `_s3` at λ=0.5 saw the same decoding order): compare arms pairwise, not as independent samples. Typical unique sequences of 8 draws: one centre 2–4 at λ=0 rising to 8 at λ=0.75; two centres 6–7 rising to 8. **Never report a row count as a design count**; raise `--temperature` or widen `--lam` for real diversity (raising `--num-seq` alone may buy none).
- **Row count is a product**: `len(--lam) × --num-seq` = **32 rows per structure** at the defaults, not 8. *Signature:* the child table is 4× the budget and `lambda` has four distinct values.
- **`ddE_mean` hides whether both centres switch.** With 2 centres a mean of 0.45 is "two at 0.45" and "one at 0.80, one at 0.10" alike, and only the second has a histidine doing nothing. **Rank on `ddE_mean`, then gate on `ddE_min`** (*signature:* `ddE_min` far below `ddE_mean`).
- **`dE_complex` and `dE_apo` go NA the moment there are two centres** (`;`-joined per-centre strings upstream). *Signature:* a filter on `dE_apo` that worked on a one-centre run silently selects nothing. Filter on `dE_apo_mean`/`dE_complex_mean` (always numeric) and read `dE_apo_centres` for the breakdown.
- **`--centers-per-type` is per TOKEN TYPE**: `2` is 4 centres under `acid_p`, 6 under `all_p` (*signature:* `n_centers` 4 where 2 was budgeted, `dE_*` NA).
- **`atomiumph_centers` has no consumer and `centre<k>_idx1` does.** A packed `;`-joined cell cannot ride through an `{expr}` island (one island resolves to exactly one integer, and the commas between islands are literal text), so a downstream `--centers` / `--fixed-positions` must name one numeric column per centre. *Signature of the mistake:* `design 'x': column 'atomiumph_centers'='5:HIS-S;10:HIS-S' is not an integer`.
- **`centre<k>_idx1` is 1-based and `centers` is 0-based**, deliberately: the audit cell keeps upstream's own numbering and the position column is in the numbering every prosapia position flag uses. *Signature of the mistake:* an off-by-one that puts a pin on the residue before the His, which `hisbond` catches as `not a histidine` and a sequence-design tool does not catch at all.
- **The positions are NA rather than wrong when they cannot be verified.** *Signature:* every `centre<k>_idx1` NA on an `OK` row; read `centres_idx1_basis` (the commonest value is `output_scope_all`).
- **A wrong-but-present chain letter is not caught.** The disjointness and non-empty checks catch typos only. *Signature:* `seq_len` is the target's length while `binder_chains` shows exactly the letter asked for.
- **A structure that fails inside a batch is silent at task level.** `design_ph.py --pdb_dir` does not abort on a per-structure failure (usually **a requested chain not in that file**): the task exits 0 and that stem has no output. *Signature:* collect prints `N of M parent(s) produced no design_ph output` and those parents have **no child rows**; read design_ph's trailing failure list in that task's `.out`. After a filtered run (`-f`) that line is expected (collect's ready set is the whole parent table).
- **Re-running the SAME dir label with different settings is refused at submit** (`grp_<g>` numbering restarts at 0, nothing prunes old output and `.meta.json` is rewritten, so collect would stamp the first run's outputs with the second run's parameters). The submit names the stems and deletes nothing: **use a new `--table-label` for a changed run.** **A SMALLER same-label re-submit is refused too, whatever the parameters** (the check is about which stems would survive untouched): **retry failed parents by re-submitting the WHOLE same filter under the same label.** *Signature:* `… already holds design_ph output for N stem(s) this submit will NOT regenerate`, N roughly the parents that succeeded.
- **A refused or abandoned submit still moves `.meta.json`'s mtime** (the driver writes its sidecar before the manifest builder runs), so collect reads `output_predates_meta == True` **leaf-wide**; advisory only. *Signature:* `True` on every row of a leaf, including rows you know are fine.
- **2–4 array tasks starting in the same second on one node die 4–8 s in** with `AttributeError: module 'torch_geometric.nn.conv.transformer_conv_TransformerConv_propagate' has no attribute 'propagate'`. Cause: PyG renders a jinja template to a shared `<PYG_HOME or ~/.cache/pyg>/message_passing/<module>.py`, so one task truncates the file while another imports it. Fixed by a per-task `PYG_HOME` exported by the **vib activation script** (`$SAPIA_VIB_SHARED_ACTIVATION/atomium.sh`, shared by `SAPIA_ACTIVATE_ATOMIUM` and `SAPIA_ACTIVATE_ATOMIUMPH`), not by the task script and not needed under Modal. *Check:* the `PYG_HOME=` line echoed into each task `.out` must be **distinct for every task**; if two tasks print the same path or none prints one, look at the shared activation script.
- **An atomium-habit noise level is refused by design.** *Signature:* `--model-noise: 'n07' is not one of n03, n05, n03_v2, n05_v2`.
- **`--set` cannot "just override" a flag** in any position (every word is inspected). *Signature:* submit-time `ValueError: --set: --binder_chains is set by this tool…`. A `*`, `?` or `[` in a value is safe (`set -f` around the call).
- **`sequence_ext` is not a sequence a predictor can read** (it carries protonation tokens): garbage or a refusal in Boltz/AF3. **`summary.csv` is ignored on purpose** (it ranks across structures, which the table does natively). **The `_p<i>` suffix is shared with `protonpottsmpnn`**: the leaf prefix is the only thing telling two child tables apart.

## What it collects

All columns leaf-prefixed `atomiumph_`; read the trust columns above first.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `status` | `OK` or `error: …` | **The only proof of success.** An error row is the single `<parent>_p0` row with counts and **no science columns** |
| `sequence` | Canonical 1-letter binder sequence | **What folding tools consume** (`-i atomiumph_sequence`, via `mkcomplex`) |
| `sequence_ext` | The same design with explicit protonation tokens | For reading and an orthogonal pKa check; **never a predictor input**; validated per row, emptied on failure with `ext_status` saying why |
| `lambda` | The sweep value this draw was sampled at | `0.0` rows are the control arm |
| `ddE_mean` | `mean(dE_complex) − mean(dE_apo)` over centres | **NA on a `--binder-alone` run.** **The selection metric, higher is better.** Equal improvement in complex and apo is a pKa shift, not a switch, which is why `dE` alone must not be ranked on. **Read WITH `ddE_min`** |
| `ddE_min` | **Minimum** per-centre ddE | Whether BOTH centres carry the switch; NA when an apo value is `nan` (no target); equals `ddE_mean` with one centre |
| `ddE_centres` | `;`-joined per-centre ddE | Same order as `centers`: **by token type then descending placement score**, not by position |
| `dE_complex` | `log p(target) − log p(comparison)` in the complex | Higher = the pinned state is more favoured with the target present. **Numeric only with ONE centre; NA with several** |
| `dE_apo` | The same on the apo binder | Strongly **positive** means the free binder's centre will not change state in the accessible pH range: nothing switches. **NA with several centres; filter on `dE_apo_mean`** |
| `dE_complex_mean` / `dE_apo_mean` | Full-precision means | Always numeric; upstream's selection recipe names `dE_apo_mean` |
| `dE_complex_centres` / `dE_apo_centres` | The verbatim `;`-joined cells | `"nan;nan"` on the apo side = no target chains: nothing in the row is a switch measurement |
| `interface_fraction` | Share of the REDESIGN SET within `--contact_cutoff` of the target | **READ IT BEFORE BELIEVING A FLAT ddE CURVE.** With the redesign set inside the binder rather than at the interface, `dE` improves equally in complex and apo and **`ddE` cannot move by construction**: λ appears to do nothing and the null result is an artefact. **Below 0.5, a flat ddE-vs-λ curve means nothing.** NaN, not 0, without target chains |
| `n_mutations` | Mutations vs the starting sequence, **as that commit defined it** | Bounded by `--max-mutations`; `0` = the seed came back unchanged. **The definition changed between commits under an unchanged name: never compare across `image_commit`**; use `n_diff_native` |
| `n_diff_native` | Hamming distance `sequence` vs `native_sequence` | Makes `n_mutations` checkable; computed from the sequences; NA when the echo is absent or lengths differ |
| `native_sequence` | The input structure's own binder sequence | Identical on every child row of a parent; under `--seed-mode native` also the seed |
| `n_his` / `n_his_native` | Histidines in the designed and input binder | With 1 centre under `his_p`/`his_s` normally `n_his_native + 1`; `n_his <= n_his_native` means it traded a His rather than adding a site |
| `centers` | `<index0>:<token>;…` | **A HUMAN-READABLE AUDIT CELL THAT NO TOOL PARSES.** The record of what was pinned, repeated across a parent's children, and the only record of `HIS-S` vs `HIS-P`. **A 0-BASED INDEX into the designed sequence, not a residue number**; no chain field exists; read `centers_resnum_key` and `centers_seq_check` first. To hand a position to another tool use `centre<k>_idx1` below, never this string |
| `centre<k>_idx1`, k = 1..12 | **1-BASED** position of the k-th centre along the designed sequence (`sequence`) **of a SINGLE binder chain**, i.e. `index0 + 1` | **THE COLUMN A DOWNSTREAM `{expr}` ISLAND CONSUMES**, one island per centre: `hisbond --centers '{atomiumph_centre1_idx1},{atomiumph_centre2_idx1}'`. Same k-order as the `;`-entries of `centers`, so the two can never disagree about which centre is `centre<k>` (one pass over `centers.json`). It is a position: no direction, neither high nor low is better. NA past the row's own centre count, past k = 12, and on every error row -- **never 0**. Written only when `centres_idx1_basis == "index0+1"`; read that first |
| `n_centres_idx1` | **must equal `n_centers`** | How many centres got a position column. Anything smaller, 0 included, means some centre has none: read `centres_idx1_basis` for which reason |
| `centres_idx1_basis` | **`index0+1`** | Why the positions are, or are not, written. **Eleven values, each a different cause — the basis is the diagnosis, so read it rather than guessing from the NAs.** Written: `index0+1` = all of them; `truncated` = only the first 12 (`--centers-per-type` raised very high). Not written, from `centers.json` itself: `unavailable` = absent, unreadable or empty; `not_a_mapping` = a centre entry was not a JSON object; `no_resnum_key` = an entry carried none of the known number keys; `not_index0` = a centre's number came from a `resnum`/`resi`-style key, so it is a residue number and `+1` would be a guess; `bad_index` = an index0 key whose value is not a non-negative integer (**`10.0` included** — one float in a future upstream JSON must not read as a numbering-convention problem). Not written, from the run: `binder_chains_unknown` = the sidecar recorded no `--binder-chains` at all (an older leaf), so the one-chain premise is unchecked — **not** `multi_binder_chain`, which would be a false diagnosis in the one column whose job is the diagnosis; `output_scope_all` / `output_scope_unknown` = **the run did not use `--output-scope binder`, so the designed sequence is every chain joined by `/` and an index into it is not a binder position**; `multi_binder_chain` = **`--binder-chains` named more than one chain, so the designed sequence is their CONCATENATION WITH NO `/`** — `centers_seq_check` passes, but an index into it is a position along no single chain, while every consumer (`hisbond --binder-chain`, `atomium --fixed-positions`) counts along one. Not written, from the row: `seq_unverified` = `centers_seq_check` is not `True` here; `seq_check_failed` / `error_row` = the row is an error row. **In every case but the first two the position columns are NA, which is the point: a plausible number on the wrong molecule is the expensive failure here** |
| `n_centers` / `n_redesign` | Counts from `centers.json` | `n_redesign` = `len(free_positions)`, matching `Stage 2: redesign positions = N` in the `.out`; NA, never 0, when absent |
| `charge_field_min` / `charge_field_max` | **ALWAYS NA** | The range is printed in the task `.out` only |
| `charge_prior` | What design_ph was given | `0.0` under `--binder-alone` |
| `metrics_name` / `metrics_centers` / `metrics_apo_weight` | Upstream's own cells | Renamed so they cannot collide with the table's `name` index or a trust column; compare `metrics_apo_weight` against `apo_weight` |
| any other metrics column | Passed through verbatim | Upstream can add a column and it appears with no code change |
| `couplings_path` | Path to `<stem>.couplings.csv`, else `""` | **A file, not columns.** Columns: `position, chain, is_binder, is_designable, residue, dist_to_center, dE_shift_<pos>_<TOKEN>, max_abs_shift{,_complex,_apo}, in_redesign_set`. Written only under `--redesign-scope coupling` |
| `path` | The parent's canonical FASTA | Provenance only. **NEVER an input column** (below) |

On disk per stem (globbed; layout `<out_dir>/grp_<g>/out/`): `<stem>.fa`, `<stem>.ext.fa`, `<stem>.metrics.csv`, `<stem>.centers.json`, `<stem>.couplings.csv` (coupling scope only), plus upstream's ignored `summary.csv`.

**Never pass `-i atomiumph_path` downstream.** It is the parent's whole `<stem>.fa` (every design on that parent plus the `_native` echo) and the **same path on every child row**, so a downstream tool would read N sequences where it expected one, include the undesigned input, and repeat the work per child. The per-design column is `atomiumph_sequence` (`sapia run mkcomplex <run_dir> -t <child> -i atomiumph_sequence …`); likewise never `atomiumph_metrics_path`.

## Filters

Two helpers every `-f` module here should start with: derive the prefix (never hard-code `atomiumph_`), and never read a bool-with-NA column with `.astype(bool)`.

```python
def _leaf(df):
    # "atomiumph" or "atomiumph_<dir_label>"; the _status/_path pair excludes the tool's own _ext_status
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
    # A bool-with-NA column as a real mask: after the TSV round trip it can hold the STRINGS "True"/"False"
    raw = df[col]
    text = raw.astype(str).str.strip().str.lower()
    return text.eq("true").where(~raw.isna() & text.ne("nan"), na)
```

`pd.Series(["True", "False"]).astype(bool)` is `[True, True]`, so `.fillna(False).astype(bool)` can admit every row it was written to exclude. The `na` argument differs **on purpose**:

| Column | NA means | `na=` |
| --- | --- | --- |
| `centers_seq_check` | the index mapping was never checked | **`False`**: unverified is not verified |
| `centers_match_metrics` | a separator this tool cannot parse | **`True`**: our blindness is not the row's fault |
| `output_predates_meta` | no sidecar to compare against | **`False`**: no evidence of staleness |
| `binder_alone` | the leaf predates the flag, i.e. a complex run | **`False`**: absent means the target WAS in the graph |

```python
# Upstream's selection recipe, with the trust gates
def apply_filter(df):
    L = _leaf(df)
    ok = df[f"{L}_status"] == "OK"
    # Status OK only says the files parsed. .notna(), NOT != "": an empty cell round-trips as NaN.
    mapped = (df[f"{L}_join_mode"] == "key") | df[f"{L}_join_checked_fields"].notna()
    centres_ok = _flag(df, f"{L}_centers_seq_check", na=False)
    # Every ddE column is NA on a --binder-alone row: exclude them explicitly, not by accident.
    not_alone = ~_flag(df, f"{L}_binder_alone", na=False)
    switch = df[f"{L}_ddE_mean"] > 0.5            # a property of the batch; see below
    apo = df[f"{L}_dE_apo_mean"] < 1.0            # _mean: dE_apo is NA with 2 centres
    every_centre = df[f"{L}_ddE_min"].fillna(df[f"{L}_ddE_mean"]) > 0.2   # the WEAKEST centre must switch too
    return df[ok & mapped & centres_ok & not_alone & switch & apo & every_centre].sort_values(
        f"{L}_ddE_mean", ascending=False)
```

**The `ddE_mean` bar is a property of the batch, not the table**: there is no published threshold; set it from the `lambda == 0.0` control of these backbones (e.g. the 90th percentile of the control's `ddE_mean`).

For a **full trust gate** add to the recipe: `_flag(..._centers_match_metrics, na=True)`, `~_flag(..._output_predates_meta, na=False)`, `placement_score_min` above a bar from this run's distribution (after reading `placement_metric`), `model_noise == "n03_v2"` (never mix noise levels), `image_commit` equal to the mode, `placement_metric == "ddE"` (a `--binder-alone` row cannot pass) and `seq_len` inside your binder's length band.

**De-duplicate on the right key.** Before FOLDING, `ok.sort_values(f"{L}_ddE_mean", ascending=False).drop_duplicates(subset=[f"{L}_sequence"], keep="first")` (a predictor cannot tell two protonation hypotheses apart). Before counting or averaging ENERGIES, `ok.drop_duplicates(subset=[f"{L}_sequence_ext"], keep="first")` instead (same `sequence`, different tokens is **two results**).

## Invariants a worker can check

None is checked by the collector (they span the table, the task `.out` and the input structure):

1. **The binder is the molecule you named**: `len(native_sequence)`, the `N` in `Native binder: N residues` in the `.out`, the residue count of the `--binder-chains` chain(s) in the input, and `seq_len` must be equal. **If they disagree the run redesigned something other than the binder.**
2. Under `--seed-mode native`, `n_diff_native == n_mutations` on every row (a difference means a position outside the redesign set moved, or `native_sequence` is not the seed); under `generate` they are expected to differ.
3. Every design id distinct, and `n_fa_entries == n_metrics_rows == n_expected`; `n_preamble_dropped` is **1** in `native`, **1 + --num-seq** in `generate`.
4. `apo_weight == metrics_apo_weight` on every row; `len(sequence_ext.split()) == seq_len`, i.e. `ext_status` is `OK` everywhere.
5. Read once from the first task's `.out`: the `design_ph.py --help` block, the resolved checkpoint path and size, the `Charge prior … field range […]` line and the per-centre `placement_score` lines, none of which reach a column.
6. **The positions against the sequence.** `sequence[centre<k>_idx1 - 1]` is the one-letter residue the k-th `;`-entry of `centers` names (`HIS-S` → `H`), and `n_centres_idx1 == n_centers`. **Only the sequence half is evidence**: it reads a third artefact (the FASTA) and an off-by-one introduced between `centers_seq_check` and the published column shows up there. The other half — that the k-th `;`-entry has index `centre<k>_idx1 - 1` — is **a tautology**, since the label and the position come off the same loop variable in one pass; it proves the ORDER is shared, which is all it is for. Report the two separately.
7. **The binder chain count**, from `atomiumph_binder_chains`: it must name exactly one chain, or `centres_idx1_basis` must read `multi_binder_chain` (two or more) / `binder_chains_unknown` (none recorded) with every position NA. A multi-chain binder with positions written would mean the gate is broken.

## What it is blind to

- **Whether the sequence folds to the designed backbone**: fold `sequence` (`mkcomplex` → `boltz`), check self-consistency with `usalign`, confidence with `boltzconf`.
- **Whether the binder still binds**: `ddE` is a *difference*, so a design can improve its switch while destroying the interface. Re-score after folding with `cms`, `pyrosetta` (`if_dG`, `if_dSASA`), `epitope` and `dockq`.
- **The actual pKa and pH window** (model log-probabilities, not a titration): confirm with `propka`.
- **Whether the chains you named are the chains you meant**: `chainsel`/`epitope` establish identity.
- **Whether a FLAT ddE-vs-λ curve is a result**: `interface_fraction` answers it (raise `--apo-weight` or tighten `--interface-cutoff` and re-run).
- **Whether the centre sits in the RIGHT-SIGN charge environment.** `placement_score_min` is the **sum** of the model and charge terms, so a strong model term hides a wrong-sign prior; only the `.out` decomposes them.
- **WHERE the centre is in the target's frame**: `centers` indexes the binder and says nothing about which target residue it faces (`epitope`).
- **The run parameters**: noise, seed, `--num-seq`, the lambda list, chains, `--centers-per-type`, `--charge-prior` and `--protect-titratable` appear in **no** upstream output file, only in this tool's `.meta.json`; a leaf whose sidecar and `atomiumph_tasks/` are gone cannot be reconstructed.
- **Expression, aggregation and developability**: `--ban-cys 1` removes one failure mode and is not a developability screen.
