---
name: boltzconf
description: >-
  How confident is Boltz about THE BINDER and its interface, rather than the whole complex? Writes binder pLDDT, `chains_ptm`, symmetrised `pair_chains_iptm`, contact-restricted interface PAE and ipSAE, with trust columns and comparability keys. Load before composing a boltzconf run or interpreting its numbers. Traps: the chain-index→letter map of the JSON is unchecked in-tool, `ipae_*` is not a BindCraft `i_pae`, and an all-NA ipSAE is a result, not a fault.
---

# boltzconf

Per design, the **binder-specific** confidence of one Boltz prediction: binder pLDDT, the binder's own `chains_ptm`, the symmetrised binder–target `pair_chains_iptm` (overall and per target chain), the PAE over the residue pairs actually in contact, and **ipSAE** (Dunbrack's interface-restricted pTM, per target chain and in both directions). **`action: update`**: annotates the table it reads, in place; `-t` is required. CPU only (the manifest builder forces `gpus_per_task = 0`, so no `-g 0`); it reads files Boltz already wrote.

## Why it exists

`boltz` collects whole-complex numbers. With N chains `iptm` is the mean over N(N−1)/2 chain pairs: on three target copies (A, B, C) plus a binder (D) that is 6 pairs, of which 3 are target–target pairs a forced template put there and exactly 1 is the binder–target pair anyone cares about. So `boltz_iptm` mostly reports how faithfully Boltz reproduced the template. The binder numbers exist in per-prediction files nothing else collects.

## The `ipsae_*` columns are a port, not a reimplementation

Source: <https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py>; paper: <https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2>. Anyone changing them must re-read that file: the d0 clamp and the normalisation are where a plausible wrong number comes from, and nothing downstream could detect it. **Agreement with the reference pins its arithmetic only**, never the residue and chain mapping it runs over.

Taken verbatim: `ptm(x, d0) = 1 / (1 + (x/d0)^2)`; the `calc_d0_array` clamp (**L raised to 26 first**, then `d0 = 1.24*(L-15)^(1/3) - 1.8`, floored at 1.0); pair selection `pae[i][j] < cutoff`, **strictly, on the raw matrix**; the **`d0res`** variant as the headline; the chain-pair value as the **max over residues** of a per-residue mean, which is **asymmetric** (A→B ≠ B→A), with the reference's "max" row the larger. The reference's `dist_cutoff` never enters the ipSAE value (it only labels its `dist1`/`dist2` report columns), so there is no `--ipsae-dist-cutoff`; `--contact-cutoff` drives `ipae_*`/`n_if_*` only.

Two deliberate differences: contributing residues are counted by **token index**, not residue number (identical whenever numbers are unique within a chain, as in Boltz output); and a chain pair with no sub-cutoff pair is **NA**, not the reference's 0.0 (`ipsae_nres_binder_<x>` is the real zero that says what happened; the binder-level aggregates keep the reference's 0.0 for such a chain). `ipSAE_d0chn`, `ipSAE_d0dom`, `ipTM_d0chn`, pDockQ, pDockQ2 and LIS are not collected.

## Premise, and when its numbers are meaningless

Each assumption is checked and **errors rather than guessing**, except the unchecked ones.

| Assumption | Signature when it fails |
| --- | --- |
| Boltz made one token per residue | `error: N pae tokens but M residues ...`. A ligand or modified residue is tokenised per atom. **Does not catch a nucleic acid** (next row) |
| The named chains are proteins | `error: chain(s) B (has DA) ... look like a nucleic acid`. A nucleotide is one token, so the row above never fires on DNA/RNA, and without this check ipSAE would use the protein `d0` floor of 1.0 instead of 2.0 and come back **inflated with every trust column clean**. An unnamed nucleic-acid chain is fine |
| Exactly one binder chain | refused at submit (`chains_ptm` and `pair_chains_iptm` are per chain) |
| The three sibling files sit beside the structure | `error: <file> is missing beside the structure`; never a partial row or another model's number |
| The structure's chain order is the order Boltz was given | **Unchecked**: the only bridge between the JSON's positional chain indices and letters. The only guard is the chain-COUNT check; `plddt_bfac_maxdiff` validates the token→residue map and nothing wider |
| The JSON describes the same complex | `error: 'pair_chains_iptm' covers K chain(s) but the structure has N` |
| The rows of `-t` are the molecules that were predicted | the submit **raises** if the structure column is absent from `-t`, unless `--allow-lineage` (a child table would score its parent's prediction under its own name) |
| Every value is finite | `error: … holds N non-finite value(s)` (a NaN would silently defeat the mapping check, since `nan > tol` is False) |

It does not know whether the chain you named is the binder, whether the pose is right, or whether the binder binds: every column is the predictor's belief about its own output. A binder confidently docked to the wrong face gets a low (good) `ipae_mean` **and a high ipSAE**; these columns rank designs that already passed the where gate (`epitope`) and the pose gate (`usalign`), they never establish either. A wrong-but-present binder letter gives a confident ipSAE about the wrong pair, and **`binder_n_res` is the only column that shows it within one row**.

## Trust columns: read these before the science

| Column | Read it as |
| --- | --- |
| `plddt_bfac_maxdiff` | **The token→residue check, and only that.** Max \|npz pLDDT ×100 − B-factor\| over every token; a right order gives **~1e-5**, > **0.05** is an error row. NA when the file has no B-factors (nothing independent checked it). Says nothing about the JSON's chain-index → letter map |
| `stem_matches_name` | `1` when the scored file is the design's own `<name>_model_0`; `0` is not an error (Boltz sanitises some names) but means the file was matched by the table cell alone |
| `input_source` | `own table` or `lineage`; `lineage` only with `--allow-lineage` and means an **ancestor's** prediction was scored |
| `n_tokens`, `n_res_struct` | Equal on every OK row by construction (a mismatch is an error); present so the equality is visible |
| `binder_n_res` | Residues of the `--binder-chains` chain. **Check against the designed length**: a wrong letter shows here and nowhere else |
| `n_target_res` | Boltz writes every residue it was given, so this must equal the expected total exactly (three 152-residue protomers: 456); any other value is a wrong input, not an unresolved loop |
| `chains` | Echo of what was asked, `binder\|targets` (`D\|A,B,C`); kept on error rows |
| `plddt_scale` | `0-1` (Boltz's own, ×100 here) or `0-100`; if it reads `0-100`, Boltz changed its output and every pLDDT threshold needs rechecking |

**The JSON side has no in-tool check**, so close it from the data:

- the `chains_ptm` of identical templated target copies should agree to about **0.05**; a much wider spread means the positional indices are not the chains you think;
- across designs `chains_ptm_binder` should correlate **positively** with `binder_plddt` and **negatively** with `pae_binder_self`; no correlation means the binder index picks another chain;
- **the strongest:** `iptm_<x>` and `chains_ptm_binder` come entirely from the JSON's positional map, while `ipsae_*` never opens it (raw PAE plus the structure's own chain order). They differ numerically but must **rank the target chains the same way across the batch**; if chain *A*'s `iptm_a` tracks chain *B*'s `ipsae_max_b`, the indices are permuted and every `iptm_*`/`chains_ptm_*` column is about the wrong chain. This cannot be closed within one row: do it once per new Boltz leaf.

## Questions, and the bar for each

| Question | Column | Conventional bar |
| --- | --- | --- |
| Is the binder a confidently folded domain? | `binder_plddt`, `binder_plddt_min`, `pae_binder_self` | Mean pLDDT 80+ for a believable fold; independent of placement |
| Does Boltz believe the interface? | `iptm_binder_max` (best target chain), `iptm_binder_top2` (two best), `iptm_binder_mean` (all) | **`_max` or `_top2` for a binder meant to span protomers.** `_mean` is diluted by every target chain the binder never touches: on three copies, 0.60 / 0.60 / 0.15 averages to 0.45, so a flat 0.5 cut on `_mean` discards exactly the spanning designs. 0.5–0.6 is a loose reference point, not a calibrated threshold. With one target chain `_top2` = `_max` |
| How sure is it where the binder sits, at the contact? | `ipae_mean`, `ipae_median` | Lower is better. The "< 10 Å" figure quoted in the `boltz` skill has no stated definition: folklore, not a threshold |
| Does the binder straddle two protomers? | `ipae_a_mean`, `ipae_b_mean`, … and `iptm_a`, `iptm_b`, … | Two decent values rather than one; a single `ipae_mean` cannot tell |
| Did it touch anything? | `n_if_pairs`, `n_if_res_binder` | `0` with status OK is a real measurement: the binder floated off |
| Same as ipTM, normalised for interface size (ipSAE)? | `ipsae_binder_max`, `ipsae_binder_top2`, `ipsae_max_<x>` | `_max`/`_top2`, never `_mean`; always read `ipsae_n0res_*` beside it |
| Which side is Boltz confident about? | `ipsae_binder_to_<x>` vs `ipsae_<x>_to_binder` | ipSAE is asymmetric (frame on the binder vs on chain *x*); a de-novo binder on a large templated target often scores higher one way |

Take every bar from your own batch's distribution, not from this table.

### There is no published threshold on ipSAE

**Do not use the ~0.3 in circulation as a pass mark.** In the paper 0.3 is the bottom of the overlap band, not a cut: *"there is overlap in the density between values of ipTM or ipSAE from 0.3 to 0.7 for the true dimers … and false dimers"*. Any gate comes from the batch's own distribution, read beside `ipsae_n0res_*` (a high score on `n0res = 2` is one residue, not an interface). The paper's numbers are anchors: non-interaction is **0.0–0.1** (worked example 0.044), true well-modelled complexes at a 10 Å cutoff are around 0.55–0.73.

**A near-zero or all-NA ipSAE across a de-novo binder batch is the paper's signature of non-interaction, not a mis-set cutoff or a tool fault.** Typical: most de-novo designs come back all-NA at a 10 Å cutoff, with per-design PAE minima of 13–26 Å (zero residue pairs under the cutoff). Raising the cutoff moves such a design from NA into the low band without making it a binder, and the two are not comparable. Diagnose a tool fault only if `n_if_pairs > 0` **and** `plddt_bfac_maxdiff` is clean **and** `ipsae_nres_binder_<x>` is 0 for every chain at 15 Å on a batch you have other reasons to believe binds.

### Comparability keys

- **This ipAE is contact-restricted and is NOT the BindCraft-style `i_pae`.** It averages only residue pairs within `--contact-cutoff` (5 Å); the mean of the whole binder×target PAE block is a different, systematically larger quantity. Compare `ipae_mean` only at the same `contact_cutoff`.
- `ipsae_pae_cutoff`: ipSAEs at different PAE cutoffs are different quantities; never pool them, and give a rerun at another cutoff its own `-l`. `--contact-cutoff` is a heavy-atom **distance** used only by `ipae_*`/`n_if_*`; `--ipsae-pae-cutoff` is a **PAE** value used only by `ipsae_*`.
- **`--target-chains` applies to the binder-level aggregates only.** `ipsae_binder_mean`'s denominator and `ipsae_binder_top2`'s width are the number of named target chains, an unreached chain contributing a hard 0.0, so `A:C` and `A` runs do not share a distribution (the same design can score 0.20 in one and 0.60 in the other); `iptm_binder_mean`/`_top2` depend on the chain set likewise. Read it off `chains`. Per-chain columns depend only on `ipsae_pae_cutoff`.
- **`ipae_*`, `pae_binder_self` and `iptm_*` use the symmetrised** `(pae[i][j] + pae[j][i]) / 2`; **`ipsae_*` use the raw asymmetric matrix** (the asymmetry is the point). A low `ipae_mean` does not imply a high ipSAE.

## Invocation and flags

```bash
sapia run boltzconf <run_dir> -t <table> -i boltz_path --binder-chains D --target-chains A:C
sapia collect boltzconf <run_dir> -t <table>

# both cutoffs stated rather than inherited
sapia run boltzconf <run_dir> -t <table> -i boltz_path --binder-chains D --target-chains A:C \
    --contact-cutoff 5.0 --ipsae-pae-cutoff 15.0
```

| Flag | Default | Note |
| --- | --- | --- |
| `-t/--table` | none | **Required** (`update`) |
| `-i/--input-column` | **`not applicable` sentinel, effectively required** | The predicted **structure** column (`boltz_path` or `boltz_<label>_path`); the sibling JSON/npz are found from that path, so not a sequence column |
| `--allow-lineage` | **off** | Reads the structure column from an ancestor table when `-t` lacks it. Off, unlike `epitope`, because a Boltz prediction belongs to one sequence; with it on, read `input_source` and `stem_matches_name` |
| `--binder-chains` | none | **Required, exactly one chain** (two are refused at submit) |
| `--target-chains` | none | **Required.** Chain mini-language (`A,B,C`, `A:C`), disjoint from the binder. **The order is the order of the per-chain columns** |
| `--contact-cutoff` | `5.0` | Heavy-atom minimum distance (Å) defining an interface residue pair. Drives `ipae_*` and `n_if_*` only; echoed into `contact_cutoff` |
| `--ipsae-pae-cutoff` | `15.0` | PAE cutoff (Å) for `ipsae_*`; a pair counts when `pae[i][j] <` this, strictly. The paper's only recommendation is "10 or 15 Å", not predictor-specific (it benchmarks no AF3 or Boltz cutoff); this default is the upper end. Echoed into `ipsae_pae_cutoff`; affects nothing else |
| `--designs-per-task` | `20` | |
| `-l/--dir-label` | empty | Own output dir and column prefix: a second boltz leaf or a second cutoff on one table |

Resources: 2 CPU, 4 GiB, 30 min per task; Modal image is `debian_slim` + `gemmi==0.7.5` + `numpy==2.5.3`. **On vib, `.env` has no `SAPIA_ACTIVATE_BOLTZCONF`** (the shared `prosapia_python.sh` would serve), so every task dies at activation until one is added; signature: `.err` holds `sapia: set SAPIA_ACTIVATE_BOLTZCONF in your .env to a tool activation script`.

## Gotchas, with their signatures

- **The input column must be the predicted structure.** An rfd3 backbone gives every row `error: confidence_<stem>.json is missing beside the structure`; the fix is the flag, not the file.
- **Siblings are found from the structure file's own stem**, not the design name: `<dir>/confidence_<stem>.json`, `pae_<stem>.npz`, `plddt_<stem>.npz` with `<stem>` = `<name>_model_0`, so a sanitised design name cannot read another design's confidence. Only `model_0`'s files are read.
- **PAE is not symmetric**: `pae[i][j]` is the error in residue *j* when the frame is aligned on *i*. `ipae_*`, `pae_binder_self` and `pair_chains_iptm` are symmetrised (typical asymmetry up to 0.08); numbers from a tool reading one direction are not comparable. `ipsae_*` are the exception; signature of getting it wrong: `ipsae_binder_to_<x>` and `ipsae_<x>_to_binder` identical on every row.
- **`ipsae_binder_top2` and `_mean` contain the reference's 0.0 for a target chain with no sub-cutoff pair, while `ipsae_max_<x>` for that chain is NA.** Deliberate: skipping a missing chain would give a one-protomer binder `_top2 == _max` and make it look spanning. All three aggregates are NA only when no target chain qualifies. Signature: `_top2` almost exactly half of `_max` on a three-copy target means one protomer reached.
- **A low `ipae_mean` with an NA ipSAE is not a contradiction.** `ipae_*` averages pairs in spatial contact regardless of PAE; ipSAE selects pairs by PAE regardless of distance. A binder touching a target Boltz is unsure about has `n_if_pairs > 0`, `ipsae_binder_to_<x>` NA and `ipsae_nres_binder_<x> = 0`, which constrains the binder→*x* direction only (`ipsae_<x>_to_binder`, hence `ipsae_max_<x>`, can still be real).
- **ipSAE is a max over residues**, so one very confident residue sets the pair's value; a 0.6 backed by `n0res = 3` is one residue being sure about three.
- **`d0res` saturates at 1.0 for every `n0res <= 26`** (the clamp raises L to 26 and the formula gives 0.9577, under the floor of 1.0): `ipsae_d0res_max_<x> == 1.0` means "small interface"; every row at 1.0 means every interface is small, not a broken column.
- **Rerun semantics bite both ways.** `sapia run … --force` deletes the per-design TSVs first, so a killed task collects as `missing`, not the old `OK`. `sapia collect` without `--force` skips rows already `OK`, so a rerun at a new `--contact-cutoff` or `--ipsae-pae-cutoff` keeps the old numbers (the echoed cutoff columns reveal it), and a later Boltz rerun that blanks a row's structure cell leaves its old `OK` and numbers in place (`filter_ready` drops the row). **Give each Boltz leaf or rerun its own `-l`** rather than re-collecting over the top.
- **Do not `sapia collect --force` over a leaf collected before the `ipsae_*` columns existed.** `ipsae_pae_cutoff` is required on an OK row and the old per-design TSVs lack it, so every row flips to `error: OK row missing ipsae_pae_cutoff` and loses its numbers; rerun the tool instead.
- **A blank own cell in the input column leaves the status EMPTY, not `error:`** (`filter_ready` drops those rows before collect) and **never borrows a parent's prediction** (the lineage is walked only when the column is absent from this table). Count non-empty `boltzconf_status` against the row count; the submit prints how many could not be resolved.
- **Per-chain columns are reset on every collect** (every `<leaf>_iptm_*`/`_ipae_*` column the table carries is blanked, including on `missing`/`error` rows and on a rerun with fewer `--target-chains`), so a stale `iptm_c` cannot sit beside numbers computed without chain C. The table is wide: 8 per-chain `ipsae_*` columns per target chain plus 4 binder-level ones.

## What it collects

Prefix `boltzconf_` (or `boltzconf_<label>_`); trust columns are in the table above. `_status` is `OK` / `error: …` / `missing`: **`OK` is the only proof of success**, and every other row is NA in every numeric column. `_path` is the scored structure (empty on non-OK rows). Per-chain columns are dynamic: the collector reads their names from each design's TSV header, and a `missing` design has none (NA).

| Column | Meaning | How to read it |
| --- | --- | --- |
| `n_res_struct` | residues in the structure, counted from atom records (runs of consecutive atoms, a repeated atom name starting a new residue), never residue numbers | TRUST; equals `n_tokens` on OK rows |
| `n_if_pairs` | interface residue pairs | `0` with OK is a real zero: the binder touches nothing at this cutoff |
| `n_if_res_binder` | binder residues in ≥ 1 pair | Small with a large `n_if_pairs` means one or two residues do all the touching |
| `binder_plddt` | mean pLDDT over binder tokens, 0–100 | Fold confidence; blind to placement |
| `binder_plddt_min` | worst binder residue | A frayed terminus shows here while the mean stays high |
| `chains_ptm_binder` | Boltz's pTM for the binder chain, 0–1 | Its own fold confidence |
| `iptm_binder_mean` | mean over **all** target chains of the symmetrised `pair_chains_iptm` | Diluted by untouched copies: do not gate a spanning binder on it |
| `iptm_binder_max` | highest per-chain value | The honest gate for "binds one protomer well" |
| `iptm_binder_top2` | mean of the two highest | The two-protomer reading; equals `_max` with one target chain |
| `ipae_mean` | mean symmetrised PAE over the **contacting** pairs, Å | **The filter column.** Lower is better; **NA when `n_if_pairs` is 0**, never 0.0 |
| `ipae_median` | same pairs, median | Less moved by a few far pairs |
| `pae_binder_self` | mean PAE inside the binder, diagonal excluded, Å | Fold confidence independent of placement; NA for a 1-residue binder |
| `contact_cutoff` | the Å cutoff used | Comparability key |
| `iptm_<x>` | symmetrised `pair_chains_iptm` against target chain *x* | One per `--target-chains` entry, letter lower-cased (`iptm_a`) |
| `ipae_<x>_mean` | `ipae_mean` restricted to pairs with chain *x* | NA when that chain has no interface pairs: the normal reading for "did not reach that protomer" |
| `ipsae_pae_cutoff` | the PAE cutoff used for every `ipsae_*` column | Comparability key, unrelated to `contact_cutoff` |
| `ipsae_binder_max` | highest `ipsae_max_<x>`, 0–1 | **The headline filter column**; NA only when no target chain has a qualifying pair |
| `ipsae_binder_top2` | mean of the two highest `ipsae_max_<x>` | Equals `_max` with one target chain; an unreached chain contributes 0.0, hence `_top2 ≈ _max / 2` for a one-protomer binder |
| `ipsae_binder_mean` | mean of `ipsae_max_<x>` over all target chains | **More diluted than `iptm_binder_mean`**: it averages a hard 0.0 (a single-protomer binder at 0.60 reports 0.60/0/0 → 0.20). **Never gate on it** |
| `ipsae_binder_to_<x>` | ipSAE, frame on the binder, chain *x* scored | NA when no binder→*x* pair is under the cutoff, never 0.0 |
| `ipsae_<x>_to_binder` | ipSAE, frame on chain *x*, binder scored | Other direction; identical values on every row would mean the raw matrix was not used |
| `ipsae_max_<x>` | larger of the two directions | The reference's "max" row; NA when both are NA |
| `ipsae_n0res_binder_to_<x>`, `ipsae_n0res_<x>_to_binder` | SUPPORT: `n0res` of that direction | How many chain-*x* residues the single best binder residue is confident about; read beside the score; whole number stored as float; NA with the score |
| `ipsae_d0res_max_<x>` | SUPPORT: d0 of the winning direction, Å | Floored at 1.0 (hit for every `n0res <= 26`); the other direction's d0 is `max(1.0, 1.24*(max(26, n0res) - 15)^(1/3) - 1.8)` |
| `ipsae_nres_binder_<x>`, `ipsae_nres_<x>` | SUPPORT: binder residues with ≥ 1 sub-cutoff PAE to *x* (reference `nres1`) / chain-*x* residues reached (`nres2`) | Real counts; `0` is exactly when `ipsae_binder_to_<x>` is NA |

**NA is never 0.** Every `error:` and `missing` row is NA in every numeric column; `ipae_*` are NA when there is nothing to average; every `ipsae_*` score, `n0res` and `d0res` is NA when that pair has nothing under the cutoff (a 0.0 would sort to the top of every "interface PAE below X" filter as if measured). Real zeros: `n_if_pairs`, `n_if_res_binder`, `ipsae_nres_*` and the binder-level `ipsae_binder_*` aggregates.

**Reading NA as 0 for `ipsae_*` is legitimate only on `OK` rows** (NA there means no pair was confident enough, which is what a non-interaction looks like); on an `error:` or `missing` row the cell is empty because the measurement failed, and reading it as 0 turns a failure into "does not interact". Mask on status first (storage stays NA, and `ipsae_nres_binder_<x>` separates a real zero from a low score):

```python
scores = df.loc[df["boltzconf_status"] == "OK", "boltzconf_ipsae_binder_max"].fillna(0.0)
```

## Filters

Cuts are placeholders in the right shape: take them from your own batch's distribution (`table-analyst` gives quantiles).

```python
# filters/ipsae_gate.py  -- rank within the batch, support beside the score
def apply_filter(df):
    ok = (
        (df["boltzconf_status"] == "OK")
        & (df["boltzconf_stem_matches_name"] == 1)          # trust columns first
        & (df["boltzconf_input_source"] == "own table")
        # binder_n_res catches a WRONG --binder-chains (the status never will: scoring the
        # target chain succeeds). Use the DESIGNED length, not a floor: >= 40 passes a target copy.
        & df["boltzconf_binder_n_res"].between(40, 60)
        & (df["boltzconf_n_if_pairs"] > 0)                  # it touched something
        & (df["boltzconf_ipsae_pae_cutoff"] == 15.0)        # one cutoff per distribution
    )
    good = df[ok]
    return good[
        (good["boltzconf_ipsae_binder_max"] >= good["boltzconf_ipsae_binder_max"].quantile(0.75))
        & (good["boltzconf_ipsae_n0res_binder_to_a"] >= 8)  # NA rows drop out automatically
        & (good["boltzconf_ipsae_nres_binder_a"] >= 5)
        & (good["boltzconf_binder_plddt"] >= 70.0)
    ]
```

Spanning two protomers (uses `_top2`, never `_mean`): `ok & ((df[["boltzconf_ipae_a_mean", "boltzconf_ipae_b_mean", "boltzconf_ipae_c_mean"]] <= 12.0).sum(axis=1) >= 2) & (df["boltzconf_iptm_binder_top2"] >= 0.4)`; the same shape works on ipSAE by counting `boltzconf_ipsae_max_<x> >= <bar>` (NA counts as False) and gating on `ipsae_binder_top2`. For an ipAE-led ranking, take `ipae_mean <= quantile(0.25)` together with `iptm_binder_max >= quantile(0.75)`.

## What it is blind to

- **Whether the binder binds.** Confidence is the predictor's belief; an orthogonal predictor or the wet lab settles it.
- **Where it bound**: `epitope` (hotspot recall, footprint, steric overlap). **Whether the predicted target landed where the template put it**: `usalign` against the template, excluding rows that did not land. **Neighbouring-protomer occupancy**: `ringfit`. **Interface area, shape, energy**: `cms` (sc, sc_area), `pyrosetta` (`if_dG`, `if_dSASA`).
- **The other Boltz models** (only `model_0` is read) and **the reference's other scores** (not collected; adding one is an edit to this tool).
- **How an ipSAE here compares with a published one.** The reference's thresholds come from AF2/AF3 models; this is the same quantity on Boltz's PAE, comparable within this workbench at one `ipsae_pae_cutoff`, not against the paper's numbers.
- **Whether the sequence is designable or expressible.**
