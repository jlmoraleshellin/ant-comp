---
name: dockrmsd
description: >-
  Did a predicted complex put the binder where the design put it, allowing for the target's symmetry? Fits the MODEL's target chains onto the REFERENCE's, applies that one transform to the model binder, and reports the binder RMSD with no further fitting — minimised over every permutation of the target chains. Takes two complex columns (`--col-a` model, `--col-b` reference) and writes `rmsd` plus a permutation audit. Load before composing a dockrmsd run, and before reading any number it wrote. Three traps: **pass `--target-match sequence`** whenever model and reference can differ in modelled extent (the default `number` errors every row when a prediction models a loop the design lacks — measured on this campaign); the binders must be EQUAL LENGTH and are never sequence-aligned (same backbone, different sequences); and this metric is UNCORRELATED with TM-score, motif RMSD and every Boltz confidence column, so no existing column substitutes for it.
---

# dockrmsd

**Docking RMSD** of a predicted complex against the design it came from.
`action: update` — it annotates the table it runs against. `-t` is **required**.
CPU only; the manifest builder forces `gpus_per_task = 0`, so callers never need `-g 0`.

The measurement, in four steps and no more:

1. Superpose the **model's target chains** onto the **reference's target chains** (Kabsch,
   reflections excluded). **Target atoms only** — the binder contributes nothing to the
   transform. That is what makes this a *docking* metric and not a fold metric.
2. Apply that one transform **T** to the **model's binder chain**. Once.
3. RMSD between the transformed model binder and the reference binder, with **no further
   fitting**.
4. Repeat 1–3 for **every permutation of the target chains** and report the **minimum**,
   naming the permutation that won.

## Premise, and when its numbers are meaningless

dockrmsd assumes:

- both files hold **the same target** (same sequence) plus **exactly one binder chain** each;
- the two binders are **the same backbone, residue for residue, in the same file order** —
  typically an rfd3/BoltzGen backbone and a prediction of a sequence designed onto it.
  Their **sequences may differ** and normally do (measured 0.467–0.604 identity on this
  campaign's designs). Sequence identity is **not** the relationship and is never required;
- the target residues can be put in correspondence by **one of two modes** (`--target-match`,
  below) — author residue numbers, or a gapped sequence alignment;
- **the target's protomers are physically interchangeable** (a true homo-oligomer). The
  permutation minimum only means something if every mapping names an *equivalent* site.

**What it does when the premise does not hold.** Unlike most tools here, dockrmsd mostly
refuses rather than returning a number:

| Violation | What you get |
| --- | --- |
| binder lengths differ | `error: binder length mismatch: ...` — never a trimmed comparison |
| a chain flag names a chain the file does not have | `error: chain(s) Z (--binder-chain) are not protein chains of the file (protein chains present: ...)` |
| **any single chain pair** corresponds over fewer than `--min-target-res` residues, on every permutation | `error: no chain permutation gave a usable target superposition: <audit>`, with `skipped_few_shared_residues(min_pair=...)` and the `per_chain=` split |
| the two targets are not the same sequence (identity < `--min-target-identity` on every permutation) | the same error, with `skipped_identity(...)` in the audit |
| the model models residues the reference lacks, under `--target-match number` | the same error. **This is the common case** — switch to `sequence` |
| the model did not reproduce the target (fit > `--max-target-fit-rmsd`) | `warn: target_fit_rmsd ... A exceeds ...` — **numbers are kept** and excluded by `== "OK"` |
| the two files correspond over less than `--min-target-res-frac` of the shorter target | `warn: only ... of the shorter target entered the fit` |
| more than one model in a file, or a duplicated author residue number | `error:` naming the count / the number — never a silent `st[0]` or first-wins |
| the target is **pseudo**-symmetric (homologous but inequivalent protomers) | **a number, and it lies.** The minimum reports the best of several inequivalent sites. Pass `--pin-mapping` |
| the binder is not rigid between the two files (a different fold) | **a number, and it lies.** The fold error is folded into the placement error. dockrmsd cannot tell the two apart; `usalign`/`posrmsd` answer the fold |

The one case where it returns a plausible number for a question nobody asked is the last
one: a binder that folded differently. Gate on `usalign` TM first if that is in doubt.

## Correspondence is per role — read this before choosing flags

This is the whole design of the tool, and the reason `dockq` could not be reused.

| Role | How residues correspond | Negotiable? |
| --- | --- | --- |
| **Target** | `--target-match number` (author residue numbers, intersected per chain pair) or `--target-match sequence` (gapped global alignment) | **yes, and the default is often wrong** |
| **Binder** | **File-order index**, equal length enforced | **no, ever** |

**Why the binder is never sequence-aligned.** The two binders are *the same backbone carrying
different sequences* — an rfd3 backbone and the sequence designed onto it, measured at
0.467–0.604 identity on this campaign. Aligning them would compare only the coincidentally
identical half. That is exactly what disqualified DockQ; do not ask for it here.

**Why `--target-match sequence` is usually what you want.** Measured on this campaign's pilot
pair:

| | model (Boltz, mmCIF) | reference (rfd3, PDB) |
| --- | --- | --- |
| target chains A/B/C | 149 residues, numbered **1–149** | 139 residues, numbered **1–139** |
| binder chain D | 118, numbered 1–118 | 118, numbered 1–118 |

Both files number every chain independently **from 1**, so this is *not* a numbering-convention
problem and **no renumbering fixes it**: the model simply *models residues the reference lacks*
— the loop `RETPEGAEAK` from position 95. Under `number` the author intersection is 139
residues of which only 97 are identical (**identity 0.698**, positions 1–94 correct and the
remaining 45 shifted by 10), so every permutation is skipped and **every row errors** — loudly,
correctly, and uselessly. Under `sequence` the loop lands in a gap and the measurement works.

There is **no fallback between the modes**, by design: a caller who asked for `number` and
silently got an alignment would have no way to know which question was answered.
**`target_match` records the mode on every row.**

Under `sequence`, `n_target_res_used` is the number of aligned pairs and
`target_seq_match_frac` is the alignment's identity — so a bad alignment is visible rather
than inferred. The aligner is Needleman–Wunsch with affine gaps and identity scoring: it is
meant for *the same protein modelled to different extents*, **not** for aligning homologues.

## Trust columns — read these before `rmsd`

| Column | Read it as |
| --- | --- |
| `target_match` | `number` or `sequence` — **which question was answered.** Read it before any target-correspondence number below |
| `pin_mapping` | `1` = no permutation search was done |
| `binder_n_res` / `ref_binder_n_res` | **Must be equal** — the row errors otherwise. Index correspondence is undefined if they are not |
| `target_fit_rmsd` | **The first number to read.** How well the target superposition itself converged, in Å. Measured 0.59–0.79 Å on the validated C3 TNFα case. A large value means the binder was moved into a frame that is not the reference frame, and the docking RMSD is meaningless **at any magnitude** |
| `best_target_fit_rmsd` | The best fit any permutation achieved. Equal to `target_fit_rmsd` unless a slightly worse-fitting permutation won the binder minimum |
| `n_target_res_used` | Residue pairs the fit used (aligned pairs under `sequence`, shared numbers under `number`). Compare with `n_target_res_ref` |
| `n_target_res_used_per_chain` | How they split, e.g. `139/139/139`. **The floor is applied per pair**, so `139/0/0` can never be accepted — if it could, the "whole-target" fit would be a single-protomer fit and `--fit-tolerance` would stop working |
| `target_resnum_range_model` / `target_resnum_range_ref` | The author numbering each file spans per target chain, e.g. `1-149;1-149;1-149` vs `1-139;1-139;1-139`. **Numbering provenance as an observation**, not something to deduce from a shrunken intersection |
| `n_target_res_model` / `n_target_res_ref` | Total target polymer residues on each side (e.g. 3×149 = 447 model, 3×139 = 417 reference) |
| `target_res_used_frac` | `n_target_res_used / min(model, ref)`. **~1.0 expected.** Well under 1 means the two files disagree on numbering and the "shared" residues may be mispaired |
| `target_seq_match_frac` | Identity over the paired target residues of the winning permutation. **~1.0 expected** — it is the same target. Below `--min-target-identity` the permutation is excluded. Under `number` this doubles as a mis-registration check, and it is an imperfect one: see the gotchas |
| `binder_seq_match_frac` | Identity between the two binders. **NOT a gate and not expected to be high**: ~0.5 is normal for a designed sequence vs its parent backbone. 1.0 means you compared a structure with itself |
| `margin` | Best minus second-best candidate permutation, in Å. Large ⇒ the symmetry choice was decisive; small ⇒ genuinely ambiguous. **`NA`, never 0, when there is only one candidate** |
| `n_permutations` / `n_candidates` | Permutations enumerated (`k!`) and how many survived the identity floor and `--fit-tolerance`. On a true C3 trimer: **6 and 3** |
| `permutation` | The winning mapping, `model>reference`, **in reference order** (`B>A,C>B,A>C` means model chain B landed on reference chain A) |
| `rmsd_all` | Every permutation, `;`-joined: `[*]<mapping>:fit=<Å>:binder=<Å>:n=<n>:per_chain=<a/b/c>:id=<identity>`, or `<mapping>:skipped_identity(...)` / `:skipped_few_shared_residues(min_pair=...)`. **A leading `*` marks a candidate** — the set that competed for the minimum — so it is read off rather than re-derived. **The minimum is auditable here rather than asserted** |
| `binder_chain` / `ref_binder_chain` / `target_chains` / `ref_target_chains` / `fit_atoms` | Echo of the flags — the cheapest way to catch a wrong chain letter |
| `model_chain_ids` / `ref_chain_ids` | Every protein chain in each file, in file order. Read this when a chain flag errors |
| `n_fit_atoms` / `n_binder_atoms` | Atoms in the fit and in the binder RMSD. With `--fit-atoms ca` these equal the residue counts; with `backbone`, 4× them |

**Why the trust block matters here specifically.** rfd3 renumbers every output chain from 1,
and predictors renumber too. A reference cut from a crystal structure (target numbered
18–155) and a prediction of it (numbered 1–138) share *no* author residue numbers — dockrmsd
errors on that rather than fitting, but the first time you see the error, the fix is
`chainsel`/renumbering upstream, not a flag here.

## The question it answers, and the threshold

**"Is this prediction's binder in the designed pose?"** — a pose gate, applied *after* a fold
gate and *before* an interface-quality gate.

No bar is conventional for this metric. On the ten pilot designs measured by hand (C3 TNFα
trimer, 139-residue reference protomers vs 149-residue predicted ones, VHH binders) the
values were:

| design (prefix `5M2J_native_trimer_vhh_`) | docking RMSD (Å) |
| --- | --- |
| `scaleup4_diff_76_a3` | **6.37** |
| `scaleup0_diff_1_a7` | 9.11 (144-residue binder — the role-flip case) |
| `scaleup3_diff_3_a7` | 16.86 |
| `scaleup6_diff_43_a8` | 17.28 |
| `scaleup5_diff_50_a1` | 19.39 |
| `scaleup2_diff_19_a5` | 22.65 |
| `scaleup1_diff_77_a2` | 24.19 |
| `scaleup8_diff_44_a4` | 26.27 |
| `scaleup7_diff_58_a8` | 26.34 |
| `scaleup9_diff_76_a7` | 30.68 |

Record the campaign's own distribution before picking a cutoff; the spread above is the
whole table's, not a validated threshold. **These values were produced by the hand-run
script this tool was built from, not by this tool** — reproducing them is the real-data
verification this tool still needs.

**This metric is uncorrelated with fold quality.** The same ten designs score TM 0.88–0.97
and motif RMSD 0.8–3.9 Å while docking 6–40 Å away, and **no Boltz confidence column
(`iptm`, `pair_chains_iptm`, `ptm`, `plddt`) identifies the good one.** Contact-overlap
counting is also not a substitute: one design reproduces 8 of 14 reference contacts while
rotated 151°. A pose gate has to be a pose measurement.

## Invocation and flags

```bash
sapia run dockrmsd <run_dir> -t table2 \
    --col-a boltz_path --col-b rfdiffusion3_path \
    --binder-chain D --target-chains A,B,C \
    --ref-binder-chain D --ref-target-chains A,B,C \
    --target-match sequence

sapia collect dockrmsd <run_dir> -t table2 \
    --col-a boltz_path --col-b rfdiffusion3_path
```

`--col-a` and `--col-b` follow the **usalign/posrmsd prefix contract**, so the three compose:
columns are `dockrmsd[_<dir-label>]_<prefix>_<metric>` with `<prefix>` from `--output-prefix`
or `<col_a minus _path>_vs_<col_b minus _path>`, and per-design files live in
`<out_dir>/<prefix>/`. **collect must be given the same `--col-a`/`--col-b` (or
`--output-prefix`) and the same `-l`**, or it raises `Results dir not found`.

| Flag | Default | Meaning |
| --- | --- | --- |
| `--col-a` | **required** | The **MODEL** complex column (the prediction). Selects which rows run; read from the `-t` table itself, **not** up the lineage |
| `--col-b` | **required** | The **REFERENCE** complex column — **the design the model came from** (e.g. the rfd3 output complex), **not** the campaign's input/native PDB; see the gotchas, the substitution does not error. Resolved per design **up the lineage**, like usalign's `--col-b` |
| `--binder-chain` | **required, no default** | Binder chain in the MODEL, exactly one id |
| `--target-chains` | **required, no default** | Target chains in the MODEL, chain mini-language (`A,B,C` or `A:C`) |
| `--ref-binder-chain` | **required, no default** | Binder chain in the REFERENCE. Often differs from `--binder-chain` |
| `--ref-target-chains` | **required, no default** | Target chains in the REFERENCE; same **count** as `--target-chains`. The pairing is **searched**, not taken from the order given |
| `--fit-atoms` | `ca` | `ca` or `backbone` (N, CA, C, O), applied to **both** the target fit and the binder RMSD. Numbers from the two settings are not comparable. `backbone` on a CA-only file is an error, not a fallback |
| `--target-match` | **`number` — often the wrong default; see above** | `number` (author residue numbers) or `sequence` (gapped alignment). **Pass `sequence` whenever model and reference can differ in modelled extent.** Affects the TARGET only; the binder is always file-order index |
| `--min-target-res` | `10` | A permutation in which **any single chain pair** corresponds over fewer residues than this is excluded. **Per pair, never pooled** |
| `--min-target-res-frac` | `0.75` | Below this fraction of the shorter target entering the fit, the row is `warn:` |
| `--pin-mapping` | off | Evaluate only the chain order as given, with no permutation search. For a **pseudo**-symmetric target |
| `--min-target-identity` | `0.8` | A permutation whose shared target residues agree below this fraction is excluded. On a homo-oligomer every permutation passes; on a hetero target only the right pairings do |
| `--fit-tolerance` | `1.0` | A permutation competes for the minimum only if its target fit is within this many Å of the best permutation's. **This is the load-bearing default** — see the gotchas |
| `--max-target-fit-rmsd` | `5.0` | Above this, `warn:` instead of `OK`. Numbers kept |
| `--max-permutations` | `720` | Refuses the submit if `k!` exceeds it (so `k ≤ 6` target chains) |
| `--designs-per-task` | `25` | Designs per task. The measurement is well under a second each, so container start dominates |
| `--output-prefix` | derived | Prefix for the results subdir and every column |

`-i/--input-column` is **refused**: rows are selected by `--col-a`.

Resources: 2 CPU, 4 G, 30 min per task. Modal image (**Modal only**) is
`debian_slim + gemmi + numpy` — a few seconds to build the first time. On vib the tool needs
a `SAPIA_ACTIVATE_DOCKRMSD` entry supplying a Python with gemmi and numpy; **no such entry
exists yet**, so on vib every task would die at activation.

## Gotchas, each with its signature

- **`--fit-tolerance` is what makes the minimum honest.** On a C3 trimer, 3 of the 6
  permutations are genuine symmetry operations (fit ~1 Å) and 3 are impostors (fit ~22 Å).
  Taking the minimum over *all 6* lets an impostor win on an accidentally small binder RMSD
  — a wrong number with a `target_fit_rmsd` that screams, if anyone reads it. The default
  excludes them. *Signature of raising it:* `n_candidates` jumps to `k!` and
  `target_fit_rmsd` becomes large while `rmsd` becomes implausibly small.
- **`--target-match number` is the default and it errors every row on this campaign's
  files.** *Signature:* every row `error: no chain permutation gave a usable target
  superposition`, with `skipped_identity(0.698)` in `rmsd_all` and
  `target_resnum_range_model` `1-149;…` against `target_resnum_range_ref` `1-139;…`. The fix
  is `--target-match sequence`, not renumbering.
- **Under `number`, the identity floor is an imperfect mis-registration check.** It is
  `matches / n_used` over the intersection, so a numbering break at position *p* of *N* scores
  `(p-1)/N`. A break at 95 of 139 gives 0.698 and is caught; the same break at ~120 of 139
  gives ~0.86, **clears a 0.8 floor**, and fits the last ~19 residues against the wrong
  partners. It degrades rather than inverting, but silently. *Signature:* nothing — which is
  why `target_resnum_range_model` / `target_resnum_range_ref` exist, and why `sequence` is
  the safer mode.
- **A pseudo-symmetric target gets a number with no warning.** `--min-target-identity` admits
  homologous protomers, and the minimum then reports the best of several *inequivalent* sites:
  a binder docked at site B of a reference designed for site A scores as though it hit the
  designed site. The tool cannot detect this. *Signature:* `n_candidates == n_permutations`
  on a target you know is not a true homo-oligomer. Use `--pin-mapping`.
- **`permutation` is in REFERENCE order.** `B>A,C>B,A>C` is one mapping written left-to-right
  over reference chains A, B, C. Reading it as model order gives the inverse mapping.
- **A wrong chain letter is an error, not a number** — but only because the letter is absent.
  If you name a *real but wrong* chain of the right length, nothing errors: the fit will be
  bad (`target_fit_rmsd` large, `warn:`) or the identity floor will refuse it. *Signature:*
  `model_chain_ids` does not contain what you expected, or `target_seq_match_frac` is low.
- **gemmi chain names are `auth_asym_id` for mmCIF.** The boltzgen skill records that author
  and label chain ids differ (on 1TUP, author E,F,A,B are label A,B,C,D). Check the real file
  before choosing the flags; `model_chain_ids` / `ref_chain_ids` are the echo.
- **A design whose `--col-b` resolves nowhere is not dropped silently.** It gets a
  `missing_input:` row written at submit time. *Signature:* `Submitting N designs` below the
  table's row count, plus a `WARNING: N design(s) not submitted for a missing --col-b`.
- **`complex_rmsd` is NOT "mostly the target", and filtering on it keeps far worse designs
  than it looks.** It is *exactly*
  `sqrt((n_fit_atoms·target_fit_rmsd² + n_binder_atoms·rmsd²)/(n_fit_atoms + n_binder_atoms))`
  — no independent information, and the terms are quadratic, so past ~2 Å the binder term
  leads. On the pilots it runs ~3.1 Å (at `rmsd` 6.37) to ~14.2 Å (at `rmsd` 30.68): a
  monotone, roughly 2×-compressed copy of `rmsd`. *Signature:* a filter `complex_rmsd < 8`
  keeps designs with `rmsd` up to ~17 Å. Kept only because its exactness is a free arithmetic
  self-check. **Never gate on it.**
- **Chain ids must be single letters**, and `expand_chain_spec` **upper-cases** them: `--target-chains a,b,c` becomes `A,B,C`. A lower-case id in the file will not be found.
- **Not the small-molecule *DockRMSD*** (Bell & Zhang 2019), which scores ligand-pose atom
  correspondence. Same word, unrelated measurement.
- **A `missing` status is not a failure.** It means no result file: the row was not selected
  by `--col-a`, or the task has not finished.
- **Both `--col-a` and `--col-b` must be full complexes.** Handing it a binder-only structure
  from `chainsel` makes the target chains absent and the row errors.
- **`--col-b` is the DESIGN the model came from, never the campaign's input PDB.** These are
  two different measurements and only one of them is docking RMSD: against the rfd3 design
  output you are asking "did the prediction reproduce *the pose we designed*"; against the
  input crystal structure you are asking "did it reproduce *the native pose*", which is a
  template-recovery question. The reference implementation used the input PDB
  (`5M2J_native_trimer_vhh_scaleup0.pdb`, numbered 9–157) for its template analysis and the
  **rfd3 design output** for the docking RMSD — match that. *Signature:* nothing errors.
  `target_resnum_range_ref` reads `9-157;9-157;9-157` instead of `1-139;1-139;1-139`, and the
  numbers are plausible but answer a different question.

## What it collects

Per-comparison status: **`dockrmsd[_<label>]_<prefix>_status`** — `OK` is the only proof of
success; `warn:` keeps the numbers but marks the fit untrustworthy. The bare `<leaf>_status`
is deliberately **not** stamped (the usalign/posrmsd convention), so filter on the prefixed
one.

| Column | Meaning | How to read it |
| --- | --- | --- |
| *trust* `binder_chain`, `ref_binder_chain`, `target_chains`, `ref_target_chains`, `fit_atoms` | echo of the flags | wrong letter here ⇒ wrong molecule measured |
| *trust* `target_match` | `number` or `sequence` | which question was answered; no fallback exists |
| *trust* `pin_mapping` | `1` = no permutation search | `1` on a true homo-oligomer means the symmetry minimum was skipped |
| *trust* `model_chain_ids`, `ref_chain_ids` | protein chains in each file, file order | what the files actually contain |
| *trust* `binder_n_res`, `ref_binder_n_res` | binder residues on each side | equal, or the row errored |
| *trust* `binder_seq_match_frac` | binder identity | ~0.5 normal; 1.0 means a self-comparison |
| *trust* `target_fit_rmsd` | winning permutation's target fit (Å) | **read first**; 0.59–0.79 Å measured |
| *trust* `best_target_fit_rmsd` | best fit over all permutations (Å) | equals the above unless a slightly worse fit won |
| *trust* `n_target_res_used`, `n_target_res_model`, `n_target_res_ref`, `target_res_used_frac` | residue pairs used, and the totals | `target_res_used_frac` ~1.0 expected; below `--min-target-res-frac` ⇒ `warn:` |
| *trust* `n_target_res_used_per_chain` | the per-chain-pair split, `139/139/139` | a lopsided split cannot be accepted — the floor is per pair |
| *trust* `target_resnum_range_model`, `target_resnum_range_ref` | the numbering each file spans | `1-149;…` vs `1-139;…` says outright that the two disagree |
| *trust* `target_seq_match_frac` | target identity over shared residues | ~1.0 expected |
| *trust* `n_fit_atoms`, `n_binder_atoms` | atoms entering each stage | `ca`: = residue counts |
| *trust* `n_permutations`, `n_candidates` | `k!`, and how many survived | C3: 6 and 3 |
| *trust* `permutation` | winning mapping, model>reference, reference order | |
| *trust* `margin` | best minus second-best candidate (Å) | NA with one candidate; small = ambiguous |
| *trust* `rmsd_all` | every permutation's fit, binder RMSD, n and identity | the audit for the minimum |
| **`rmsd`** | **the docking RMSD (Å)** | the headline; meaningless unless `target_fit_rmsd` is small |
| `max_dev` | largest single binder-atom deviation (Å) | a hinge or a terminus pulling the RMSD up |
| `complex_rmsd` | target + binder under the winner (Å) | an atom-count-weighted quadrature mean of `target_fit_rmsd` and `rmsd`; always **smaller** than `rmsd`, carries nothing new, never gate on it |
| `sup_path` | the whole model moved into the reference frame (mmCIF) | open it beside the reference to see the pose |

Every non-`OK`/non-`warn` row carries **`NA`, never `0`**, in `rmsd`, `max_dev`,
`complex_rmsd` and `sup_path`. An **error row keeps the trust columns it managed to
establish** (chain ids, binder lengths, the permutation audit) — those are the diagnosis.

On disk: `<run_dir>/<table>/dockrmsd[_<label>]/<prefix>/<name>.tsv` and `<name>_sup.cif`,
plus the sub-manifests under `<prefix>_tasks/`.

## A filter a caller can write

```python
# Keep designs whose binder is in the designed pose, on a target the model
# actually reproduced. 8.0 A is a placeholder: set it from this campaign's own
# distribution and record the number in the campaign log.
PREFIX = "dockrmsd_boltz_vs_rfdiffusion3"

def apply_filter(df):
    return df[
        (df[f"{PREFIX}_status"] == "OK")
        & (df[f"{PREFIX}_target_match"] == "sequence")
        & (df[f"{PREFIX}_target_fit_rmsd"] < 1.5)
        & (df[f"{PREFIX}_target_res_used_frac"] > 0.9)
        & (df[f"{PREFIX}_rmsd"] < 8.0)
    ]
```

Note the ambiguity check a careful caller adds on a highly symmetric target:
`df[f"{PREFIX}_margin"].fillna(999) > 2.0` — a small margin means two symmetry copies fit
almost equally, and the reported minimum is a coin toss between them.

## What it is blind to

- **Whether the binder folded at all.** It measures placement, and a misfolded binder's fold
  error shows up as placement error. → `usalign` (TM-score), `posrmsd` (specific positions).
- **Whether the pose is *good*, only whether it is the *designed* one.** A prediction may
  find a better interface than the design. → `cms`, `pyrosetta` (`if_dG`, `if_dSASA`),
  `polyv`.
- **Which target residues are actually contacted.** → `epitope`, `ifacegeom`.
- **Clashes, chemistry, side chains.** CA (or backbone) atoms only. → `epitope` (`n_clash`),
  `pyrosetta`.
- **Whether the reference itself is a sensible design.** It is taken as ground truth.
- **Anything per-batch** — correlations, distributions, bimodality. One value per design is
  all it writes; batch questions belong to `table-analyst`.

## Status before use

Written 2026-10-15. **218 synthetic checks pass**, and they prove the arithmetic and the
contracts only. `tool-reviewer` has passed it for real-data verification and its findings are
landed. It has **not** been run on real data. Until it reproduces the ten values above under a
`verification` label, do not gate a campaign on its numbers.

### The reference implementation, read and compared (2026-10-15)

The hand-run script that produced the ten values above
(`scratchpad/worker_clean/arms_analysis_script.py`) was read directly. **Its fit basis is
identical to this tool's**, so **exact agreement is the acceptance bar, not rank agreement**:

| | the script | this tool |
| --- | --- | --- |
| atoms | CA only (`ca_by_num()`, `r['CA'][0]`) | `--fit-atoms ca` (default) |
| fit | one Kabsch over all three protomers pooled (`np.vstack([... for c in 'ABC'])`) | one Kabsch over all mapped chain pairs pooled |
| permutation choice | **top 3 by target fit** (`perms.sort()`, `perms[:3]`), then the minimum binder RMSD | permutations within `--fit-tolerance` of the best fit, then the minimum |

The two permutation rules differ in form and **cannot disagree on a true C3**: three mappings
fit at ~1 Å and three at ~22 Å, so "best three" and "within 1 Å of the best" are the same set.
They *would* diverge on a target where four or more mappings fit comparably — if that ever
happens, this tool's rule is the correct one and the script's `[:3]` is the hard-coded
assumption.

Its residue correspondence is also exactly what `--target-match sequence` should derive: the
script normalised both files to original 5M2J author numbering, **uniformly** for the
prediction (`num + 8`, 1–149 → 9–157, because the prediction *has* the loop) and **piecewise**
for the design (`num + 8` up to 94, `num + 18` beyond, because the design *lacks* it). That
piecewise map is the gapped correspondence, derived there from construct knowledge and here
from the sequences — prefix 1:1, suffix shifted by exactly 10, which is what
`test_inserted_loop` asserts. If the aligner ever disagreed with that shape, that test is
where it would show.

### The verification run, when the Volume is free

Ten pilot designs under a `verification` label, `--target-match sequence`, acceptance at
**exact agreement (<0.05 Å) on all ten**. Plus two deliberate failures: a chain letter absent
from the file, and — the sharper one — `--ref-target-chains A,B,D --ref-binder-chain C`,
feeding the binder in as a target, which is the only way to watch the identity floor fire on
real sequences rather than synthetic ones.

### Still open after that

`--min-target-identity 0.8` and `--fit-tolerance 1.0` have never been measured under
`--target-match sequence` on real files. Neither has `--pin-mapping`, `--min-target-res-frac`
or the per-chain-pair floor.

### Two settled decisions, recorded so they are not quietly reopened

- **`--target-match` defaults to `number` on purpose, even though `sequence` is right for the
  only known use.** The failure modes are asymmetric: forgetting the flag makes *every row
  error* with an explicit identity message — loud, instant, cheap. A `sequence` default
  meeting files that genuinely should have matched by number would let the aligner quietly
  find *some* mapping. Strict default, deliberate opt-in to the looser mode.
- **There is no `n_models` column.** A file with more than one model is refused, with the
  count in the error message. A column would read `1` on every surviving row — a vacuously
  passing trust column, which is exactly what the house contract forbids.
