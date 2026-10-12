---
name: binder-campaign
description: >-
  How to sequence and judge a de-novo binder campaign: defining a bindable target, the order to apply the gates, which tool answers which question, and what each metric is blind to. Load at the start of any binder design against a specific target and before interpreting any interface number. Tool flags live in the per-tool skills.
---

# Running a binder campaign

The per-tool skills tell you how to invoke things. This tells you **in what order, and what the numbers mean**.

## The one-paragraph version

Define a target a binder could physically reach, and **delete** everything else rather than avoiding it if necessary (membrane domains for example). Design, then judge in a fixed order: **did the target land → did the binder fold → did it stay in its designed pose (and at the epitope) → is the interface any good**. "Folded" and "stayed in its designed pose" are different gates — a binder can fold perfectly and dock in the wrong place, so never let a fold score or a hotspot-contact count stand in for the pose gate (`usalign --mm 1` vs the rfd3 backbone; see Gate 3). Never read a later number before the earlier gate passed, because broken predictions produce the most attractive interfaces in the batch. Keep spend late: everything before the co-fold should be cheap. 

---

## 1. Define the target

1. **Establish the biological assembly before anything else.** How many chains, identical or distinct, what ligands, what is membrane-embedded. A brief naming "subunits A and B" may mean two identical protomers of a homo-oligomer, which turns the job into a composite epitope across a seam.
2. **Trim to what a binder can reach**, and  *delete* the rest — do not merely keep hotspots away from it. Identify a membrane belt from direct lipid contacts or other ligands like glycans and require the methods to agree.
3. **Cutting manufactures new decoy surface. Cap it.** Slicing protomers out of an oligomer exposes faces that do not exist in reality, and a predictor docks to them as readily as to real ones. Take this into account when choosing epitopes.
4. **Prove the remainder is still a domain.** Count heavy-atom contacts and backbone H-bonds between the segments you kept vs. what you removed. If the retained parts only pack through the deleted piece, the target is a fiction. (A good sign: the trimmed unit is *more* compact than the original.)
5. **A gap you create is a real gap.** If trimming leaves two segments per chain separated in space, they must be given to any predictor as **two chains**. Fusing them into one sequence makes the predictor close the gap by distorting the domain — silently, with confident output.

## 2. Choose the epitope

- Work out the **approach vector** first. A "side chain points radially outward" test is wrong for a flat-bottomed particle whose accessible face points along the symmetry axis; residues there can read as "inward" while having rel_sasa 0.98.
- Score candidate patches on **balance** — per-chain exposed area, `min/max` — not total area, or you will pick a one-protomer patch.
- Single-linkage clustering **percolates** on a continuous exposed surface (one cluster, 73  residues, 49 Å across). Use seed-centred footprint patches with non-max suppression.
- **Reject on physics even when the score is good.** A patch whose centroid plus the binder's  footprint reaches the membrane plane, or that abuts a surface your own trimming created, is  disqualified regardless of balance.
- Check the epitope **fits**: patch longest dimension vs. binder size, and the arc of target   available before the next symmetry-related site.
- Present the pool to the user with reasoning and let them choose.

## 3. The pipeline for plain binder design

Plain binder design is a simple process: design a backbone, design a sequence, and judge the result.

```
rfd3/bc2_traj          backbones, hotspots on BOTH target chains          → table0
epitope                hotspot recall/avoidance on the diffusion pose     → table0
-> Gate 1: did the diffusion land? (hotspot recall)
chainsel               extract designed binder chain                      → table0 (update)
proteinmpnn/atomium    binder chain only, target chains fixed             → table1 (child)
mkcomplex              rebuild the complex: target chains + binder        → table1
boltz                  co-fold, forced template                           → table1
epitope                hotspot recall/avoidance on the PREDICTED pose     → table1
usalign                --mm 1 co-fold vs rfd3 backbone: POSE gate         → table1
-> Gate 2: did it stay in its designed pose? (TM-score, hotspot recall)
[boltz]                binder alone — cheap fold pre-screen               → table1 (Optional)
[usalign]              --mm 0 vs designed binder                          → table1 (Optional)
-> Gate 2.1: did it fold? (TM-score)
-> Further analysis: pyrosetta, csm...
```

The `usalign --mm 1` pose gate and the `epitope`-on-prediction gate are **both** Gate 3 and
**both required** — the first asks whether the binder stayed in its designed 3D pose, the second
which target residues it actually ends up touching. They answer different questions (see Gate 3).

**`mkcomplex` is not optional.** `proteinmpnn_sequence` is the binder monomer alone; feeding it to a predictor folds the binder in isolation and answers a different question with no error.

**The binder-alone screen is optional**. It is *standard practice* as a self-consistency test. Its unique value is being the **undiluted** fold measurement. Its an pseudo-orthogonal check on the co-fold, it can be used to further rank the best binders after the co-fold.

## 4. The gate order — do not reorder this

### Gate 1 — did the diffusion land?

Rfdiffusion is a stochastic generator. Even when setting hotspots, it can produce a diffusion pose that is far from the target. The first gate is to check whether the diffusion pose landed on the target at all, before doing any expensive co-folding or scoring.

### Gate 2 — did the binder stay in the designed pose after prediction?

"Did it stay" (is the binder in the 3D pose the generator designed, in the target frame) is NOT the same as "where" (which target residues does it contact). A binder can clip two chosen hotspots while sitting well off its designed pose, and epitope-contact recall will pass it. Run the pose gate; do not let hotspot recall stand in for it.**

**The `usalign --mm 1` route**: predicted complex vs the rfd3 backbone complex, `--col-a boltz_<label>_path --col-b rfdiffusion3_path` (lineage-resolved). Both hold the *same* target, so the superposition locks onto it and the binder's displacement shows up in that frame. **Do NOT threshold on TM or RMSD here — both are target-dominated and the RMSD actively lies.** `--mm 1` keeps only the residues that structurally align: a binder that drifted drops out, leaving a *target-only* alignment that reports a **beautiful low RMSD at ~0.4 Å**. **Read `Lali − L_target`** — binder residues that actually aligned in the target frame (`_Lali` over all chains, minus the fixed target length): near the binder length = it stayed, near 0 = the binder is gone. *Measured on `binder_nohis` (818 rows):* `Lali − L_target` had **median 1, max 67 against 62–88-residue binders** — the binder was off its designed pose for essentially the whole batch. There, `Spearman(RMSD, Lali−L_target) = +0.945` and `Spearman(Lali−L_target, hotspot_recall) = +0.671`, so **low RMSD coincided with the drifted (target-only) rows** — i.e. ranking by low RMSD selects the *worst* designs. `Lali − L_target` is the honest column; cross-check with `epitope`-on-the- prediction hotspot recall, which it tracked at +0.671.

### Gate 2.1 — did it fold?

Binder chain only, `usalign --mm 0` against the designed backbone. `--mm 1` on a complex is
dominated by the large fixed target and its TM/RMSD read ~1.0 regardless; `--ter 2` isolates the
*first* chain while the binder is usually last. Use `chainsel` to extract. **This answers only
"did it fold" — not "did it fold *where it was designed to*".

## 5. What each metric cannot see

- **Self-consistency cannot see solubility.** A 4.5%-charged, 70%-TSVIG β-sandwich folds at
  pLDDT 0.96. Check composition directly: % charged, net charge, % aromatic, % β-branched.
- **CMS and SC cannot see chemistry.** Geometry only.
- **`buried_unsat` is mostly not the interface.** *Measured:* totals 25–44 while
  `if_delta_unsat` (the interface share) was only 4–16. A prediction that composition problems
  would show up as interface buried-unsat was **falsified**; they showed up as packing quality
  (`sc` 0.36–0.53, `packstat` < 0.6).
- **Whole-complex confidence cannot see the binder.** See the `boltz` skill: `iptm` averages
  N(N−1)/2 pairs, nearly all template-forced target–target.
- **The one cheap honest signal:** binder mean pLDDT **alone vs. in complex**. A real interface
  raises it. *Measured:* all 27 designs *lost* 7.5–23.3 points — a confidently-folded domain the
  predictor is confidently unsure where to put.

## 6. Judgment

- **Keep the spend late.** Recon, trimming, epitope choice, backbones, sequences and a fold
  screen should all be cheap. Pilot the co-fold on 2 designs before scaling.
- **Pilot the target definition itself**, not just the parameters — run two target definitions
  head to head on 2 designs and let a measurement pick, rather than arguing.
- **Verify chain mappings empirically.** Superpose over candidate permutations; the right answer
  separates by orders of magnitude (0.050 Å vs ≥16 Å). Generators renumber from 1.
- **Demand the check that would have failed** for any number gating a spend.
- **State expectations before results arrive.** In this campaign one advance prediction was right
  (the decoy) and two were wrong (β designs would fail — they were the best; per-pair ipTM would
  re-rank — it did not). The wrong ones were the more useful.
- **Write the campaign up.** See `campaigns/`.
