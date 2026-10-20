---
name: propka
description: >-
  Does a designed histidine resist protonation IN THE COMPLEX more than in the free binder (the pKa shift a pH switch needs), and how much binding does that cost at a chosen pH? Runs PROPKA3 on the complex and on the binder chain(s) alone; writes `pka_complex`, `pka_free`, `dpka`, the SWITCH RATIO (`switch_ph6p5_all` = K_D(6.5)/K_D(7.4)) and `stab_*` fold-stability loss. Load before a propka run or reading any `dpka`, F or `stab_*`: a negative shift means nothing unless `pka_free` is near the release pH, and the free state is a RIGID unbinding.
---

# propka

**Custom tool** (`tools/propka/`, not bundled with prosapia). **`action: update`**: annotates the table it ran against and mints no child. **`-t`, `-i/--input-column`, `--binder-chains` and exactly one centre source (`--resnums` or `--auto-centres`) are required** (no honest default for any of them). CPU only: the builder forces `gpus_per_task = 0`, so `-g 0` is never needed. Runs on **both executors**: `modal_image.py` ships (`debian_slim` + `propka==3.5.1`, nothing else), and on vib `SAPIA_ACTIVATE_PROPKA` in `.env` points at the lab's shared `propka.sh`. Resources: 1 CPU, 4 GB, 2 h per task. `propka` is not installed in this local workspace.

It answers two questions. **(1) For each designed histidine, what pKa does PROPKA3 predict in the complex, what for the same residue in the free binder, and how big is the shift, and what weakening of binding does that buy at a pH you care about? (2) For a buried core histidine with `dpka` ~ 0, does protonation destabilise the binder's own fold as the pH drops (`stab_*`)?** Campaign sequencing and shortlist rules live in `ph-switch`, the run loop in `running-a-step`, commissioning in `verifying-a-tool`.

## Premise, and when its numbers are meaningless

One structure file per design holding the **binder chain(s)** (`--binder-chains`, one or more) and, normally, a partner (the target), plus the designed histidines: named by a spec (`--resnums`, whose `{expr}` islands can read table columns) or taken as every His of the binder chains (`--auto-centres`). Two PROPKA3 runs per design, on temporary copies of that one file:

| run | what PROPKA sees |
| --- | --- |
| complex | every chain of the file, hydrogens stripped |
| free | **only the binder chains' records of the same file**, same coordinates, hydrogens stripped |

`dpka = pka_complex - pka_free`.

**The free state is a RIGID UNBINDING**: the binder chain cut out of the complex, not a relaxed or re-predicted apo binder (`free_state_model == chain_only_same_coordinates` on every row). If the real free binder rearranges (a side chain that repacks, a helix that frays once the target leaves), `pka_free` is the pKa of a conformation that does not exist and `dpka` inherits the error. The tool cannot see it.

**Several binder chains are one free state, and a file with no partner is refused.** `--binder-chains A,B` cuts both chains out together, so a homo-oligomeric binder is compared with the whole oligomer alone, not one protomer. Centres are `(chain, residue number)`: with one binder chain a bare number is enough, with several write `A:45,B:45` (`--resnums`). **A file with no ATOM chain outside the binder chains is an `error` row** (`chain_ids_found` shows it): the complex run would be the free run and every `dpka` 0. **`--free-only` is the explicit single-run mode**: one PROPKA run on the binder chains alone, a binder-only file allowed, `free_only` = 1 on every row. Everything that needs the complex is then **NA by construction, not missing data**: `centre<k>_pka_complex` / `_dpka` / `_buried_complex` / `_n_coupled_complex`, `min_dpka`, `n_centres_coupled`, `per_centre_pka_detail`, `term_groups_match`, `binder_term_groups_complex`, `n_propka_warnings_complex`, the F, switch, release, dips and `ph_release_saturated` columns and `ph_half_complex_*`. What stays: `centre<k>_pka_free` / `_buried_free` / `_n_coupled_free`, `binder_his_pka` (complex slot `NA`), `f_ceiling_*`, `ph_half_free_*` and the whole `stab_*` block (free state only anyway).

**PROPKA is an empirical predictor, run here on predicted and relaxed coordinates.** Its usual accuracy is about **one pK unit** (the tool measures none of its own), on top of a model input. **A `dpka` of -0.3 is noise; a ranking over many designs is worth more than any single row.** **It does not assign the tautomer**: PROPKA protonates the ring itself and reports one pKa per histidine; which nitrogen carries the proton comes from the pose and is `hisbond`'s `centre<k>_tautomer`. The two tools answer different halves and must both be read. **It is not a titration and not a binding measurement**: the real test is binding at pH 7.4 against 5.5-6.0, and everything here is a prior over which designs deserve it.

## The reading rule for `dpka`

For one titratable residue, `dG_bind(pH) - dG_bind(high pH) = RT ln[ (1 + 10^(pKa_free - pH)) / (1 + 10^(pKa_complex - pH)) ]`, so **in the direction bind at neutral pH, release in acid**:

- **`dpka` must be NEGATIVE**: the complex disfavours the protonated His more than the free binder does, so protonating it costs binding energy. A positive `dpka` binds *better* in acid. The sign convention reverses for a binder meant to bind in acid; nothing in the columns knows which direction you want.
- **`pka_free` must be near or above the release pH (about 5.5-6.0).** A shift from 4.0 to 3.0 does nothing at pH 5.5: the His is protonated in neither state, the log term is ~0 and the switch never fires. **Read `centre<k>_dpka` and `centre<k>_pka_free` together, always.**

Rough size at 298 K: `RT ln 10 = 1.36 kcal/mol` per pK unit, so a `dpka` of -1.0 with `pka_free` well above the release pH is about 1.4 kcal/mol, roughly a 10-fold weakening, and one pK unit is also PROPKA's own error bar.

## How much binding is lost at a pH: the F columns and the switch ratio

For one histidine the two-state independent-sites weakening factor relative to the high-pH limit is `F_i(pH) = (1 + 10^(pKa_free_i - pH)) / (1 + 10^(pKa_complex_i - pH))`, and independent sites multiply (`F = prod_i F_i`). It is **arithmetic on the two pKa values already collected, with no extra PROPKA run** (`f_model == two_state_independent_sites`). Three properties decide how to read it:

- **F -> 1 at high pH and -> `10^(-dpka)` at low pH.** `10^(-dpka)` is the most the switch can ever deliver (`dpka` -1.0 saturates at 10-fold, -0.5 at ~3.2-fold), approached asymptotically.
- **`F_i(pH) <= 1 + 10^(pKa_free_i - pH)` for any complex pKa.** That bound, the `f_ceiling_*` columns, is set by the free binder's own protonation. If the ceiling at your release pH is 2, no `dpka` can give a 3-fold release.
- **For a complex His that stays neutral, `F(pH) = X` at `pH = pKa_free - log10(X - 1)`** (3-fold at about `pKa_free - 0.3`, 10-fold at about `pKa_free - 0.95`). The release pH is set mainly by `pka_free`; `dpka` only decides how much of the ceiling is reached. For a release at pH 6.5 with X = 3 you need `pka_free >= 6.8` and a complex pKa roughly a unit below that pH.

**The switch ratio is what a pH switch is judged on**: `K_D(pH) / K_D(pH_ref) = F(pH) / F(pH_ref)`. A binder is released from the pH it was holding at, not from the high-pH limit, so F at the release pH alone overstates the switch whenever the free binder is already part-protonated at the reference (a His with `pKa_free` 7.2 and `pKa_complex` 5.0 has `F(7.4) = 1.62` already spent, so `F(6.5) = 5.83` is really a **3.59-fold** switch). `switch_ph<tag>_designed` and `switch_ph<tag>_all` are that ratio with `--reference-ph` (default 7.4) as denominator. **Gate on the ratio at the release pH, not on F.** A ratio below 1 means binding is stronger at that pH than at the reference. The reference pH has no ratio column (it would be 1).

**Two site sets, always side by side.** `*_designed` uses only the centres named by the centre flags; `*_all` uses **every** histidine of the binder chains (the residues of `binder_his_pka`, at full precision). An undesigned His still titrates, so `_all` is the honest one for "will this binder release" and `_designed` isolates the designed switch. **If any binder His has no pKa in either state, every `_all` column is NA as a block** and `n_binder_his_no_pka` says how many (a refusal, never a product over a silently smaller set). A binder carrying more histidines than designed centres is common, and an extra His titrating the other way can cancel the switch, invisible in any per-centre column (hence `binder_his_pka`, `n_binder_his`).

**Assumptions that can all be wrong at once**: a rigid unbinding; independent sites (doubtful for two close histidines, which `n_centres_coupled` flags); only these histidines titrate (no acid, terminus or target residue is in the model); PROPKA's ~1 pK error, which is a 10-fold error in the low-pH limit; and F is a predicted factor, not a measured K_D ratio.

## Trust columns: read these before any pKa

| column | read it as |
| --- | --- |
| `status` | `OK`, `error: ...` or `missing`; the only proof of success. Every science column is NA on an error row; only the echo and trust facts established before the failure survive. **NA is never 0.** |
| `propka_version` | The release behind every number. Parameters ship with the package: two leaves on different versions must not be pooled. |
| `binder_chains` | The chain(s) that *were* the free state (echo, comma-joined). A wrong letter makes every `dpka` another molecule's, and nothing else in the row says so. |
| `binder_chains_n_res` / `expect_binder_len` | ATOM residues in those chains (total), against what `--expect-binder-len-col` supplied. `expect_binder_len` is NA when the flag was not used: **pass it**, so a disagreement becomes an error row. |
| `binder_chains_n_het_res` | HETATM residues on the binder chains. `0` for a clean binder. They are copied into both PROPKA inputs but never looked up as a centre. |
| `binder_chains_sequence` | The one-letter sequence of the binder chains, `/` between chains (non-standard residues `X`). **The only post-hoc evidence of WHICH molecule the letters named**, whether or not a length was expected: compare it with the designed sequence (`atomiumph_sequence`). It is what catches a target chain mistaken for the binder under `--auto-centres`, where no residue number anchors the choice. |
| `chain_ids_found` | Every chain ID in file order. Catches an absent or mis-cased `--binder-chains` and a binder-only input. |
| `centres_source` | Where the centres came from (echo): `resnums:<spec>` or `auto`. |
| `free_only` | `1` when only the free state was run (`--free-only`): every complex-dependent column is NA by construction; `0` otherwise. |
| `centres_n` / `centres_reported` | Centres resolved for this design / how many got per-centre columns (`min(centres_n, 4)`). **Below `centres_n` means the row has more centres than column sets**: every other column, `per_centre_pka_detail` and `binder_his_pka` included, still covers all of them. |
| `centres_his_ok` | Centres that were a **histidine**, uniquely numbered, in the binder chains. **Always equals `centres_n` on an OK row** (a non-His or repeated number fails the design); on an error row, how many were confirmed before the failing one. |
| `centres_pka_ok` | Centres that got **both** pKa values. Equals `centres_n` on an OK row; on an error row, how many centres' pKa stage completed before the failure (it localises it). |
| `binder_term_groups_complex` / `_free` | The terminal groups PROPKA built on the binder chains in each run, sorted and comma-joined (`C-78D,N+1D`); empty when none. Kept on error rows (complete before the centre loop that produces the likeliest error). |
| `term_groups_match` | `1` when those strings are identical, `0` when not, NA only if the row failed before both runs returned. **Gate on it.** `C-` comes from an OXT atom and is symmetric; a `0` means the complex and free state got different charged termini on the binder, a file artefact that puts a spurious shift on any N-terminal-proximal centre. |
| `n_h_stripped` / `n_h_stripped_guessed` | H/D atoms removed from the copies (non-zero for a Rosetta pose, 0 for AF3/Boltz; the input is never modified) / how many were identified by the **name** fallback rather than the element columns 77-78. Stripping is belt and braces on a standard file. `guessed` is `0` on any file with element columns; non-zero means the tool is guessing from names (ATOM records only, so a HETATM metal named `HG`/`HO`/`HF`/`HE`/`HS` with a blank element is never deleted). |
| `n_his_renamed` | HID/HIE/HIP/HSD/HSE/HSP residues renamed to HIS (PROPKA only knows `HIS`; with hydrogens gone the variants are the same heavy atoms). |
| `n_propka_warnings_free` | WARNING-and-above records in the **free** run, a property of *this* binder. **Gate on it** (expected 0; usual messages `Unexpected number of atoms in residue`, `His group does not seem to contain a ring`). Warning text goes to the task `.err`, tagged with design and run. |
| `n_propka_warnings_complex` | The same for the complex run, which prechecks the whole assembly. **Do not gate on it**: one truncated residue in a shared target makes it non-zero in every design and a `== 0` clause returns an empty frame. A constant across the table is the target's. |
| `ph_values` / `release_factors` / `reference_ph` | The pH list, factor list and reference pH this leaf used, in the canonical token form the column names are built from. Leaves that differ are not comparable column for column. |
| `f_model` / `free_state_model` | Constants (`two_state_independent_sites`, `chain_only_same_coordinates`), written so nobody has to infer them. |
| `n_binder_his_no_pka` | Binder histidines with no pKa in one run. `0` normally; non-zero makes every `*_all` column NA for the row. |

## Invocation and flags

```bash
sapia run propka <run_dir> -t <table> -i pyrosetta_path \
  --binder-chains D \
  --resnums '{hisbond_relax_centre1_resnum},{hisbond_relax_centre2_resnum}' \
  --expect-binder-len-col atomiumph_seq_len -C 40

# no upstream columns needed: every His of the binder chain, or a literal list
sapia run propka <run_dir> -t <table> -i <structure column> --binder-chains D --auto-centres
sapia run propka <run_dir> -t <table> -i <structure column> --binder-chains D --resnums 45,52

# the free state alone (a binder-only file is allowed)
sapia run propka <run_dir> -t <table> -i <structure column> --binder-chains A --auto-centres --free-only

sapia collect propka <run_dir> -t <table>
```

| flag | default | meaning |
| --- | --- | --- |
| `-i/--input-column` | **none** (`"not applicable"` sentinel) | The structure column of the **complex**, resolved up the lineage (the builder raises and lists the table's `*_path` columns if omitted). **CIF is accepted**: converted to PDB at submit time with prosapia's `ensure_pdb` (gemmi, auth chain letters and numbers kept), cached under `run_dir/.cif_to_pdb`, because **PROPKA 3.5.1 reads `.pdb` only** (`Unknown input file type`) and the tool image holds propka alone. `<leaf>_path` records the file PROPKA read. |
| `--binder-chains` | **required** | One or more chains (`D`, `A,B`, `A:C`). Together, alone, they are the free state and the centres are looked up in them. |
| `--resnums` | one of two | The centres as a spec: `45,52` (one binder chain) or `A:45,B:45`, in `auth_seq_id` numbering. `{expr}` islands are resolved per design up the lineage (integers, column names, `+ - * //`; the nearest non-null value wins), so table columns are `'{hisbond_relax_centre1_resnum},{hisbond_relax_centre2_resnum}'` and a multi-chain binder `'A:{colA},B:{colB}'`. **A design whose island resolves nothing (an NA, a missing column) is an error row**: a varying number of centres per design cannot be expressed, use `--auto-centres` for that. |
| `--auto-centres` | one of two | Centres = **every histidine** of the binder chains, in file order; no upstream tool needed. More than 4 is fine: the first 4 get per-centre columns (`centres_reported`), every other column covers all of them. **It has no residue number to anchor the chain choice**: a target chain named as the binder still gives a clean `OK` row of the target's histidines, so pass `--expect-binder-len-col` and check `binder_chains_sequence`. |
| `--free-only` | off | Run the free state alone; see the premise above. |
| `--expect-binder-len-col` | none | Column holding the expected total residue count of the binder chains: an integer (`atomiumph_seq_len`) or a one-letter sequence (`/` breaks allowed). A mismatch is an error row, and **a check asked for and unresolvable is also an error row**. |
| `--designs-per-task` | `10` | Designs per task. Each is **two** PROPKA runs and **the complex run is the slow one** (its non-covalent coupling scan is O(n^2) in the titratable groups of the whole assembly, so cost follows the target; a ~500-residue complex is expected to take seconds to under a minute per run). The default is not tuned: time the complex and free runs separately on a verification leaf before a big fan-out. |
| `--ph-values` | `7.4,6.5,6.0,5.5` | pH values at which F is reported; each from 0 to 14 inclusive, **at most 6**, no duplicates (refused at submit and in the worker). Each adds **ten columns, five at the reference pH**: `f_`, `f_ceiling_` and `switch_` each in `designed` / `all`, plus `stab_dg_fold_`, `stab_loss_` and its `_designed` / `_other` shares, named with the pH (`6.5` -> `f_ph6p5_designed`; the decimal point becomes `p`, one decimal is always kept: `6.0` -> `ph6p0`); at the cap of 6 that is up to 58 pH columns. No extra PROPKA run. |
| `--reference-ph` | `7.4` | The denominator of the switch ratio: the pH the binder is meant to hold at. **Must be one of `--ph-values`** (refused otherwise). Every other pH gains `switch_ph<tag>_designed` and `switch_ph<tag>_all`. **The column name does not carry the reference**, so `switch_ph6p0_all` is `F(6.0)/F(7.4)` from one run and `F(6.0)/F(6.5)` from another: changing it on a collected leaf means a new `-l` unless every row is re-run. |
| `--release-factors` | `3,10` | Weakening factors X for the release-pH columns; each **> 1**, at most 4, no duplicates. Each adds four columns, `ph_release_x<X>_designed` / `_all` and `f_dips_below_x<X>_designed` / `_all` (`2.5` -> `x2p5`). |
| `--keep-propka-input` | off | Keep the two PDBs PROPKA read as `<out_dir>/<name>_propka/complex.pdb` and `free.pdb`. Under `--free-only` only `free.pdb` is written. **Turn it on for a verification leaf**: they are what `propka3 <file>` must reproduce by hand. |
| `-l/--dir-label` | empty | Score two structure columns onto one table (`-l relax`, `-l af3`); every column becomes `propka_<label>_*`. **Changing `--ph-values`, `--release-factors` or `--reference-ph` on a collected leaf needs a new `-l` unless every row is re-run**: columns are a union, old-pH columns go stale (orphans, below) and nothing but each row's own echo distinguishes "never asked" from "asked and failed"; for `--reference-ph` it is worse because the `switch_*` and `stab_loss_*` columns keep their names while changing meaning. |

**The centre flags take RESIDUE NUMBERS, not 1..L indices.** Exactly one of `--resnums`, `--auto-centres` is required (zero or two is refused at submit). They coincide only when the chain is numbered from 1, which `hisbond` records per design as `centres_idx1_eq_resnum`: check it is 1 before taking `hisbond_<leaf>_centre<k>_resnum` as the source on a new table.

**What the tool does to the file.** The input is never modified. PROPKA sees a temporary copy under the task's `TMPDIR` with exactly four edits, all counted or constant: (1) every H/D atom dropped (element columns 77-78; when blank, and on ATOM records only, the first alphabetic character of the atom name); (2) HID/HIE/HIP/HSD/HSE/HSP renamed HIS; (3) a bare `TER` at every chain change (PROPKA re-arms its N-terminal search at a `MODEL`, a `TER` and an OXT atom, so chains ending in OXT already get an `N+` each and this edit is defensive: it covers chains ending without OXT, and a chain written just before the binder that ends on the *same residue number* the binder starts on, which suppresses the binder's `N+` in the complex run while the free run always has one; `term_groups_match` reports what PROPKA actually did); (4) first model only. HETATM records are copied as they are: PROPKA ignores `HOH H2O SO4 PO4 PEG EPE TRS`, anything else (a ligand, a metal) becomes a group and **can shift a pKa**, and `binder_chains_n_het_res` is the only hint (binder chains only).

## What it collects

Leaf-prefixed `propka_*` (`propka_<label>_*` with `-l`). **Spelling is British throughout**: `propka_center1_pka_complex` is a `KeyError`. **`hisbond` spells the two counts the American way** (`hisbond_centers_n` / `hisbond_centers_his_ok`) while `propka` has `propka_centres_n` / `propka_centres_his_ok`: a filter reading both leaves needs both spellings.

**Per centre k = 1..4, in the order the centres were named (file order under `--auto-centres`):**

| column | meaning |
| --- | --- |
| `centre<k>_chain` / `centre<k>_resnum` | The chain and residue number that centre's pKa belongs to. Under `--auto-centres` they are read off the residues found; under `--resnums` they are **the spec echoed back**, so `propka_centre1_resnum == hisbond_..._centre1_resnum` proves nothing about the mapping (the real identity check is that the residue is a histidine). NA past `centres_reported` and on error rows. |
| `centre<k>_pka_complex` / `_pka_free` | PROPKA's pKa for that histidine in the complex / in the free binder. |
| `centre<k>_dpka` | `pka_complex - pka_free`. **The filter column**; negative is the wanted sign. |
| `centre<k>_buried_complex` / `_buried_free` | PROPKA's buried fraction for that group, **percentage 0-100**. A centre that buries on binding is the usual cause of a negative shift. `0` is a measurement (fully exposed); NA is an absence. |
| `centre<k>_n_coupled_complex` / `_n_coupled_free` | **Read before `dpka`.** Groups PROPKA found **non-covalently coupled** to that His in each run (the ones it marks `*`). Non-zero means its pairwise-independent model is breaking down for that group, and `dpka` is a difference of two such values so the error does not cancel; two close designed histidines are the canonical trigger. `0` is wanted and a measurement. Distinct from *covalent* coupling, which withholds the pKa and errors the row. The list is copied from the first conformation only (pKa terms are averaged over all), so the count is exact for single-conformation files (every Rosetta and AF3 output) and approximate with altlocs. |

**Per row:**

| column | meaning |
| --- | --- |
| `n_centres_coupled` | Designed centres with a non-zero coupling count in either run: the one-clause form of the two columns above. `0` wanted. Counted inside the centre loop, so NA on an error row (unlike the terminal-group columns). |
| `min_dpka` | The most negative `dpka` over the designed centres, the one-number **ranking** key. It says nothing about *which* His moved or whether it is protonated at the release pH: sort with it, never gate on it alone. |
| `n_binder_his` | Histidines on the binder chains in total, designed centres included. |
| `binder_his_pka` | **Every** His of the binder chains, file order: `;`-joined `<chain>:<resnum>:<pka_complex>:<pka_free>`, two decimals. `NA` in a slot means no pKa for that residue in that run (incomplete side chain, or a group discarded as coupled); this does not fail the row, whose guarantee is about the designed centres. |
| `per_centre_pka_detail` | One `;`-joined entry per centre, **complex run only**: `<chain>:<resnum>:model=<v>:buried=<pct>:desolv_reg=<v>:desolv_re=<v>:sc_hb=<v>:bb_hb=<v>:coulomb=<v>`. The pKa is `model` plus the five terms, each the sum over PROPKA's determinants of that kind: it says *which physics* moved the pKa, not which partner. |

**The pH dependence** (names carry the pH, `ph6p5` is 6.5, and the site set `designed` / `all` always comes last). All floats, NA on an error row, never 0:

| column | meaning |
| --- | --- |
| `f_ph<tag>_designed` / `f_ph<tag>_all` | The weakening factor at that pH over the **designed** centres / over **every** binder His: `F = prod_i (1 + 10^(pKa_free_i - pH)) / (1 + 10^(pKa_complex_i - pH))`. `F = 1` is no effect; `F = 5` is 5-fold weaker than at high pH. `_all` is NA as a block when `n_binder_his_no_pka > 0` and is usually the one to gate on. |
| `switch_ph<tag>_designed` / `switch_ph<tag>_all` | **The headline for a switch**: `F(pH)/F(reference_ph)`, the predicted `K_D(pH)/K_D(reference)`, for every pH except the reference. `switch_ph6p5_all = 4` means 4x weaker at 6.5 than at 7.4. **Below 1 is the wrong direction.** `switch_*_all` is NA as a block like `f_*_all`. |
| `f_ceiling_ph<tag>_designed` / `f_ceiling_ph<tag>_all` | `prod_i (1 + 10^(pKa_free_i - pH))`: the largest F any complex pKa could produce, set by the free binder alone. A ceiling below your target X says the free binder is not protonated enough at that pH. Always `>= f_*`. |
| `ph_half_free_<set>` / `ph_half_complex_<set>` (`<set>` = `designed`, `all`) | The **titration midpoint** of the site set: the pH at which half the histidines are protonated (`sum_i 1/(1 + 10^(pH - pKa_i)) = n/2`, bisection on [0, 14] to 1e-6), "the pKa of the whole binder". **For one site it equals that site's pKa** (with one designed centre `ph_half_free_designed` is `centre1_pka_free`: a free cross-check of the site set); **for several it is a summary that hides them** (5.0 and 8.0 give the same midpoint, 6.5, as 6.5 and 6.5). NA when the half-point is outside 0-14. Read F, not the midpoint, to decide a release. |
| `ph_release_x<X>_designed` / `_all` | The **highest pH in [3, 9] at which F reaches X**, on a 0.01 grid scanned downwards (never optimistic by more than the grid). It does **not** say F >= X everywhere below, since F is not monotonic in pH. **A ranking key, not a gate.** |
| `f_dips_below_x<X>_designed` / `_all` | `1` when F falls back below X somewhere between that release pH and pH 3, `0` when it stays at or above X all the way down, NA exactly when the release pH is NA. A site with positive `dpka` contributes a factor below 1 that keeps falling with pH, so a mixed binder can cross X on the way down and cross back (sites free 6.5 / complex 5.0 and free 4.0 / complex 6.0 peak near pH 5.5 and have collapsed by pH 4). **`1` means the natural reading of the release pH is false for that row.** |
| `ph_release_saturated` | `;`-joined `<factor>_<set>` tags whose F is already >= X at pH 9.0 (e.g. `x3_all`). **Empty is normal and is a measurement, but arrives as NA**: test `df[f"{L}_ph_release_saturated"].fillna("") == ""`, never `== ""` alone. It distinguishes the two NAs of `ph_release_x*`: NA **with** the tag means the answer lies above the window (weakened X-fold even at pH 9, usually a sign something is wrong); NA **without** it means F never reached X in [3, 9]. If `x10_<set>` is listed, `x3_<set>` must be too. |

**Fold stability vs pH: the `stab_*` block**, a different question: not "does protonation weaken binding" but "does protonating the designed histidines destabilise the BINDER'S OWN FOLD as the pH drops", the only mechanism left for a buried core His whose `dpka` is ~0 and switch ratio ~1. It is PROPKA's own folding free energy (reference `neutral`) on the **FREE run only**, evaluated directly at each pH of `--ph-values` (no third run). **Sign: PROPKA's `dG_fold` is LOWER = MORE STABLE** (its optimum is the minimum); the loss columns are oriented so that **positive = worse**.

| column | meaning |
| --- | --- |
| `stab_dg_fold_ph<tag>` | `dG_fold` of the whole free binder at that pH, kcal/mol, written at every pH including the reference. **Lower is more stable.** |
| `stab_loss_ph<tag>` | `dG_fold(pH) - dG_fold(reference pH)`, the reference being `--reference-ph` (7.4 by default; **the name does not carry it**). **POSITIVE = LESS STABLE at the lower pH.** No column at the reference. **The column to rank on.** |
| `stab_loss_ph<tag>_designed` / `stab_loss_ph<tag>_other` | The same difference summed over the **designed centres only** (each group's own folding term) / over **everything else titratable** in the binder (Asp, Glu, Lys, Arg, Tyr, Cys, termini **and any non-designed histidine**). The split is designed-vs-rest, **not His vs non-His**. |
| `stab_share_check` | max abs(his + other - total) over the reported pHs; ~0 (1e-15) means the shares reconstruct PROPKA's total. It cannot see a centre that matched no group (the next column does). |
| `stab_n_his_groups` | Designed centres that matched a histidine group. **MUST EQUAL `centres_n`**: the worker raises into `stab_error` when it does not, so a designed share of 0 is never published for a centre that matched no group. |
| `stab_n_groups_titrated` | Titratable groups PROPKA built for the free binder. |
| `stab_profile_check_maxdiff` | max abs(direct evaluation - PROPKA's `get_folding_profile`) over pH 5, 6, 7, 8, 9. **Plumbing only**: both sides call the same function. |
| `stab_opt_ph` / `stab_opt_dg` | The pH of the profile's optimum (its **minimum**) and the `dG_fold` there. NA when PROPKA does not determine them. |
| `stab_error` | Empty (-> NA) when the block succeeded, otherwise why not. **The block is NON-FATAL**: a failure leaves status `OK` with every `stab_*` NA, so a `stab_loss_*` filter drops the row silently. Collect prints `propka: stab_error: N of M OK rows carry a stab_error` (always, including 0): **read it before concluding "no stability loss"**. |

**What `dG_fold` is blind to**: one fixed conformation with no response to protonation; an unfolded reference of isolated groups at their model pKa (His 6.50); PROPKA's empirical electrostatics only; a small difference of larger terms each carrying a ~1 kcal/mol class error. **Rank by `stab_loss_*`; never read an absolute value as a folding free energy, and absolute `dG_fold` is not comparable across designs.** `stab_*` from relaxed coordinates (`-i pyrosetta_path`) and from AF3 coordinates are not comparable in one ranking: split the table on the structure column first. Typical: on core-His binders of ~80 residues, `stab_loss` at pH 6.5 / 6.0 / 5.5 is about +1.1...+1.2 / +2.0...+2.6 / +3.2...+4.0 kcal/mol, 84-97 % of it from the designed histidines (a His whose folded-state pKa is depressed costs stability against the unfolded reference, model pKa 6.50, as soon as the pH approaches that reference, so a loss at pH 6.5 appears even for a centre whose own pKa is ~5).

## Filters a caller can write

**The leaf prefix is `propka_` only when the run used no `-l`**; with `-l relax` every column is `propka_relax_*` and a pasted filter raises `KeyError`. Set `L` once.

```python
L, PH = "propka", "ph6p5"          # PH is the tag of the RELEASE pH


def apply_filter(df):
    # THE GATE: the outcome asked for (release at that pH), with the trust clauses.
    # NA fails every clause, so error rows drop silently: COUNT THE STATUSES.
    assert df[f"{L}_reference_ph"].dropna().nunique() == 1, "mixed reference pH"
    ok = df[f"{L}_status"] == "OK"
    mapped = (df[f"{L}_centres_his_ok"] == df[f"{L}_centres_n"]) & (
        df[f"{L}_centres_pka_ok"] == df[f"{L}_centres_n"])
    # Tool's own echo of the expected length: NA without --expect-binder-len-col and
    # NaN != NaN, so this admits NOTHING then.
    right_chain = df[f"{L}_binder_chains_n_res"] == df[f"{L}_expect_binder_len"]
    # FREE run only: the complex run prechecks the whole target and would empty the frame.
    quiet = df[f"{L}_n_propka_warnings_free"].fillna(1) == 0
    uncoupled = df[f"{L}_n_centres_coupled"].fillna(1) == 0
    same_termini = df[f"{L}_term_groups_match"] == 1
    # THE SWITCH RATIO (K_D(6.5)/K_D(7.4)), not F; 'all', not 'designed'.
    switches = df[f"{L}_switch_{PH}_all"] >= 5.0
    return df[ok & mapped & right_chain & quiet & uncoupled & same_termini & switches]
```

A two-proxy variant gates one switch centre K on `centre<K>_dpka <= -0.5` **and** `centre<K>_pka_free >= 6.0` with the same trust clauses (the thresholds are decisions, not table properties). **Gate on `f_ph<tag>_<set>` or the switch ratio; rank with `ph_release_*` and `min_dpka`.** `ph_release_x3_all >= 6.5` looks like the same question and is not (the release pH is only the highest pH at which F reaches X; `f_dips_below_x3_all == 1` is the row that passes it and is below X at 6.0). **Before believing either, look at `f_ceiling_<tag>_all`**: below your X, the design cannot release and `dpka` is beside the point.

**Columns for a pH no longer in `--ph-values` are ORPHANS.** A re-collect writes only the columns the new run produced, so a stale `propka_f_*_ph6p4` keeps its old values, valid for that pH as computed then but not this leaf's output and flagged only by `ph_values`. Drop them or keep them knowingly; never mix them with current-pH columns in one comparison. **An empty frame is not a result**: no design passed, every row errored, or a clause compared against NA (the `right_chain` line without `--expect-binder-len-col` is the easy one). Always report `df[f"{L}_status"].value_counts()` with the survivor count.

## Traps, each with its signature

- **PROPKA's own proton filter is better than it looks, but incomplete.** It guesses the element from columns 13-14 and takes only the first character when the name fills all four (`propka/atom.py`), so ` HA `, `HD11`, `1HB`, `2HD1` are all read as H and dropped. What survives is a left-justified two- or three-character hydrogen name such as `HZ1 ` (read as `Hz`, kept as a pseudo-atom) and **deuterium** (element `D`). This tool strips from the element columns, so stripping is harmless and idempotent on a well-formed file. *Action:* run `propka3` on the raw file and on `--keep-propka-input`'s `complex.pdb` and report whether they differ.
- **A residue number is not a 1..L index.** `--resnums` reads `auth_seq_id`; handing it `hisbond`'s `centre<k>_idx1` on a chain not numbered from 1 looks up the wrong residue. *Signature:* `error: centre 1: residue 12 of chain D ... is ALA, not a histidine` if lucky, and a plausible pKa for the wrong His if the file happens to hold one there. Check `hisbond_centres_idx1_eq_resnum == 1` before equating them.
- **The binder's terminal groups can differ between the two runs: a file artefact, not chemistry.** Do not predict it from the file; read `term_groups_match`. *Signature:* `term_groups_match == 0` with `binder_term_groups_complex` short of an `N+` entry that `_free` has.
- **A centre PROPKA cannot place fails the WHOLE row.** An incomplete side chain, or a group PROPKA discards as **covalently** coupled (a path an ordinary His essentially never takes), gives `error: centre 2 (residue 57): propka discarded ...` with every science column NA. The coupling that does fire for histidines is **non-covalent**, suppresses nothing, and is reported as `centre<k>_n_coupled_*` / `n_centres_coupled` to be filtered on. `OK` therefore means every designed centre has both pKa values; the cost is that one bad centre hides a good one (the status text names the residue).
- **A row whose own structure cell is blank is never collected.** `filter_ready` drops it, so it keeps an EMPTY `propka_status` (or whatever an earlier run left, a stale `OK` included) instead of `error:`; the builder prints `could not resolve ...` at submit. Count non-empty statuses against `len(df)`, and gate on `propka_status == "OK"`.
- **A converted CIF is cached by path, not by content.** `ensure_pdb` keys `run_dir/.cif_to_pdb/<stem>.<digest>.pdb` on the source path alone, so a predictor re-run `--force` into the same CIF path is read from the SUPERSEDED cached PDB with status `OK`. *Signature:* `propka_path` under `.cif_to_pdb/` older than the input CIF. Delete the cached PDB (or use a new run label) after regenerating structures; the cache is shared with `pyrosetta`.
- **A multi-chain binder needs chain-qualified centres and a total length.** With several binder chains every `--resnums` token must read `<chain>:<resnum>` (hisbond writes bare numbers, so use `'A:{hisbond_relax_centre1_resnum},B:{hisbond_relax_centre1_resnum}'`), and `--expect-binder-len-col` is compared with the TOTAL over the chains. Both fail loud.
- **Two histidines at one residue number (insertion codes) are refused.** PROPKA drops the insertion code, so two groups would share a label: a designed centre there errors, a non-designed His there gets `NA:NA` in `binder_his_pka`.
- **Waters and ligands are not stripped.** A ligand in the complex is absent from the free run unless it sits on a binder chain: a real difference between the states, not an error. *Signature:* `binder_chains_n_het_res > 0`, or a `dpka` no burial change explains.
- **A two-character chain ID is refused**: PROPKA reads the chain from column 22 only, and a file with something in column 21 errors with `carries a chain ID in column 21`.
- **A leaf written by an older version of the tool is not readable.** Columns were renamed (the site set now comes last in the pH columns, `binder_chain*` became `binder_chains*`, `stab_his_*`/`stab_other_*` became `stab_loss_*_designed`/`_other`) and `centre<k>_chain`, `centres_reported`, `centres_source` and `free_only` were added; such a row collects as `error: OK row missing ...`. Re-run it.
- **Tool discovery blast radius**: `spec.py` imports `run_propka.py`, which imports `propka_worker.py` (and the collector imports it for the shared parsers and `MAX_CENTRE_COLUMNS`), so a syntax error or module-level exception in the worker breaks discovery for every `sapia` command. The worker must import only the standard library at module level (`import propka` stays inside functions) because the same file runs in the tool image, which has no prosapia.
- **Submit proof per executor** is in `running-a-step`; an exit code proves nothing.

## Invariants a worker can check

Run on the first real leaf, under a `verification` label, with `--keep-propka-input`.

**Real checks:**

1. **The free run is reproducible by hand.** `propka3 <out_dir>/<name>_propka/free.pdb` gives the same pKa as `centre<k>_pka_free` (PROPKA prints 2 dp, the column carries more). Independently, write a chain-only, hydrogen-stripped PDB of the same input with any tool and run `propka3`: same number. **The complex run likewise**, against `complex.pdb` and `centre<k>_pka_complex`.
2. **Hydrogens change nothing.** Run on the Rosetta PDB and on a hydrogen-stripped copy: every pKa identical and `n_h_stripped` 0 on the second. Separately run `propka3` on the raw file and on the kept `complex.pdb` and report whether they differ.
3. **The decomposition adds up.** Parse `per_centre_pka_detail`: `model + desolv_reg + desolv_re + sc_hb + bb_hb + coulomb` equals `centre<k>_pka_complex` to within **0.01** (PROPKA's own definition of a pKa), the cheapest check that the per-centre mapping reads the right group.
4. **The free file is what it claims**: the distinct `(resnum, icode)` among ATOM records of `free.pdb` equals `binder_chains_n_res`, and no record has a chain other than the binder.
5. **Terminal groups agree**: `term_groups_match == 1` on every OK row; report how many are 0 and which centres sit near the binder's N-terminus.
6. **Physics sanity**: a surface His with no target contact has `dpka` ~ 0 and `buried_complex` ~ `buried_free`; a centre that buries at the interface shows both a burial rise and a shift. Large `dpka` with equal burials deserves suspicion.
7. **The F columns recompute by hand** from `centre<k>_pka_free` / `_pka_complex` (`f_ph<tag>_designed`) and from `binder_his_pka` (`f_ph<tag>_all`; agrees only to rounding, because `binder_his_pka` is 2 dp and the column is full precision, so a few borderline designs cross a threshold). **The dips flag agrees with F**: a row with `f_dips_below_x3_all == 0` and a release pH at or above 6.5 must have `f_ph6p5_all >= 3` (the grid scan and the direct evaluation are separate paths); a row with the flag at 1 need not, which is the point of the flag. **One-site sanity**: with one designed His, `ph_half_free_designed` equals that site's `pka_free` to 1e-6, and `ph_release_x<X>_designed` equals `pka_free - log10(X - 1)` to within the 0.01 grid when the complex pKa is far below (grid answer at or just below the exact crossing, never above).
8. **The stability block against PROPKA's own print**: on every OK row `stab_n_his_groups == centres_n` and `stab_error` is empty; on two designs run `propka3` on the kept `free.pdb` and compare its "Free energy of folding" rows at pH 7.00 and 6.00 with `stab_dg_fold_*` from a run with `--ph-values 7.0,6.0`: they agree to the 2-dp print (<= 0.005).
9. **Statuses**: report `value_counts(dropna=False)` and the row count; rows whose own structure cell is blank never reach collect (below) and only the NaN bucket counts them.
10. **The non-designed share of F, recomputed from the text column.** On every `OK`, non-`--free-only` row with `centres_n <= 4` and `n_binder_his_no_pka == 0`: parse `binder_his_pka` into `{(chain, resnum): (pka_complex, pka_free)}`; every `(centre<k>_chain, centre<k>_resnum)` must appear in it with pKa values matching `centre<k>_pka_complex` / `_pka_free` to 0.005; `len(binder_his_pka) == n_binder_his`; and for each pH tag `f_ph<tag>_all / f_ph<tag>_designed` must equal the product over the REMAINING entries of `(1 + 10^(pka_free - pH)) / (1 + 10^(pka_complex - pH))` to about 1e-2 relative (the text is 2 dp). The three products come from separate loops, so it fails if the centre loop and the all-His loop ever key on different residues or the `_all` set is smaller than `n_binder_his`.

**Plumbing and tautologies** (pass by construction; they test the implementation, not the physics): `f_ceiling_* >= f_*` (the denominator is >= 1); `F` is 1 at pH 14 and approaches `prod_i 10^(-dpka_i)` at pH 0 (asymptotes, agreeing to ~1e-5 relative); `ph_release_x3_* >= ph_release_x10_*` where both are non-NA; `switch_ph<tag>_<set> == f_ph<tag>_<set> / f_ph<ref_tag>_<set>` (same products) and ratios multiplying like F; no ratio column for the reference and `reference_ph` in `ph_values` on every row; `stab_share_check <= 1e-6` and `stab_profile_check_maxdiff <= 1e-6`.

**Break-it controls:** a wrong `--binder-chains` (`A` on a `D` binder) must give an error row for every design with `chain_ids_found` filled; a non-His residue number must give `error: centre 1: residue N of chain D ... is <RESNAME>, not a histidine`; a wrong `--expect-binder-len-col` must error; `--ph-values 15` and `--ph-values 7.4,7.4` must be refused at submit and in the worker; `--release-factors 1` must be refused; `--reference-ph 6.0` with a list lacking 6.0 must be refused at submit naming the list; `--ph-values 6.5 --reference-ph 6.5` must collect cleanly with only `ph6p5` columns and no `switch_*`. A binder-only file without `--free-only` must give the `holds no chain outside the binder chains` error row, and with it an `OK` row with `free_only` 1 and the complex columns NA. `--binder-chains A,B --resnums 102` must error asking for `<chain>:<resnum>`. Two centre flags, or none, must be refused at submit.

## Standalone worker

The worker runs on one structure with no table, manifest or prosapia, the fastest way to check a single design or a deliberately bad input. It needs `propka` importable (run it in a container of the tool's image, or any env with `pip install propka==3.5.1`):

```bash
python tools/propka/propka_worker.py --structure /path/to/complex.pdb --binder-chains D \
  --resnums 45,52 --name mydesign --keep-propka-input --out-dir /tmp/propka_check
```

It writes `/tmp/propka_check/mydesign.tsv` with exactly the columns the collector reads plus `mydesign_propka/{complex,free}.pdb`; `--expect-binder-len <int-or-sequence>` adds the length check, `--resnums auto` takes every His, and `--free-only` runs the free state alone.

## What it is blind to

- **Conformational change on unbinding**: the free state is rigid by construction.
- **The tautomer**, and so whether the contact a protonated His would break is a donation or an acceptance: `hisbond`'s `centre<k>_tautomer`, `per_centre_role_target` and `any_switch_candidate`.
- **Which partner moved the pKa**: `per_centre_pka_detail` gives sums per interaction type, never determinant labels.
- **Cooperativity between titratable groups.** PROPKA discards a member only for **covalent** coupling; for the non-covalent coupling two nearby histidines produce it still reports a pKa and merely flags the group (`n_coupled_*`), which this tool cannot correct. A two-His switch acting together is outside the single-residue equation and outside this tool.
- **Everything about binding that is not protonation** (affinity, kinetics, expression, aggregation, the target's own titration) and **the target's histidines and acids** (only the binder chains' histidines are reported; `hisbond`'s `_n_acidic_env` / `_n_basic_env` is the target-side context).
- **Everything the F model leaves out** (rigid unbinding, independent sites, only these histidines titrate) and **whether PROPKA is right**: its ~1 pK error becomes a 10-fold error in the low-pH limit of F.
