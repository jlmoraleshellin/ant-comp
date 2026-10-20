---
name: binder-campaign
description: >-
  How to sequence and judge a de-novo binder campaign: defining a bindable target, the order to apply the gates, which tool answers which question, and what each metric is blind to. Load at the start of any binder design against a specific target and before interpreting any interface number. Tool flags live in the per-tool skills.
---

# Running a binder campaign

The per-tool skills tell you how to invoke things. This tells you in what order, and what the numbers mean. For binders that must release in mild acid, add `ph-switch` after the gates below.

## The one-paragraph version

Define a target a binder could physically reach, and delete everything else rather than avoiding it (membrane domains, for example). Design, then judge in a fixed order: did the target land, did the binder fold, did it stay in its designed pose (and at the epitope), is the interface any good. "Folded" and "stayed in its designed pose" are different gates: a binder can fold perfectly and dock in the wrong place, so never let a fold score or a hotspot-contact count stand in for the pose gate (`usalign --mm 1` against the rfd3 backbone, Gate 2). Never read a later number before the earlier gate passed, because broken predictions produce the most attractive interfaces in the batch. Keep spend late: everything before the co-fold is cheap.

## 1. Define the target

1. **Establish the biological assembly first**: how many chains, identical or distinct, what ligands, what is membrane-embedded. "Subunits A and B" may be two identical protomers of a homo-oligomer, which makes the job a composite epitope across a seam.
2. **Trim to what a binder can reach and delete the rest**, do not merely keep hotspots away from it. Identify a membrane belt from direct lipid contacts or other ligands (glycans) and require the methods to agree.
3. **Cutting manufactures decoy surface.** Slicing protomers out of an oligomer exposes faces that do not exist in reality and a predictor docks to them as readily as to real ones. Take that into account when choosing epitopes.
4. **Prove the remainder is still a domain**: count heavy-atom contacts and backbone H-bonds between the kept segments against what you removed. If the kept parts only pack through the deleted piece the target is a fiction.
5. **A gap you create is a real gap.** If trimming leaves two segments per chain separated in space, give a predictor two chains; fusing them makes it close the gap by distorting the domain, silently and with confident output.

## 2. Choose the epitope

- Work out the **approach vector** first: a "side chain points radially outward" test fails on a flat-bottomed particle whose accessible face points along the symmetry axis.
- Score patches on **balance** (per-chain exposed area, `min/max`), not total area, or you pick a one-protomer patch; single-linkage clustering percolates on a continuous surface, so use seed-centred footprint patches with non-max suppression.
- **Reject on physics even when the score is good**: a patch whose centroid plus the binder footprint reaches the membrane plane, or that abuts a surface your own trimming created. Check the epitope fits (patch size against binder size, arc before the next symmetry-related site), present the pool to the user with reasoning, and let them choose.

## 3. The pipeline for plain binder design

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
-> Gate 2: did it stay in its designed pose? (Lali − L_target, hotspot recall)
[boltz]                binder alone, cheap fold pre-screen                → table1 (optional)
[usalign]              --mm 0 vs designed binder                          → table1 (optional)
-> Gate 2.1: did it fold? (TM-score)
-> Further analysis: pyrosetta, cms, ...
```

The `usalign --mm 1` pose gate and the `epitope`-on-prediction gate are both Gate 2 and both required: the first asks whether the binder stayed in its designed 3D pose, the second which target residues it actually touches. **`mkcomplex` is not optional**: `proteinmpnn_sequence` is the binder monomer alone, and feeding it to a predictor folds the binder in isolation, answering a different question with no error. The binder-alone screen is optional but standard: the undiluted fold measurement and a pseudo-orthogonal check on the co-fold.

## 4. The gate order, do not reorder

**Gate 1, did the diffusion land?** RFdiffusion is stochastic: even with hotspots it can produce a pose far from the target. Check that the diffusion pose landed on the target at all before any expensive co-folding or scoring.

**Gate 2, did the binder stay in the designed pose after prediction?** "Did it stay" (the binder in the pose the generator designed, in the target frame) is not "where" (which target residues it contacts): a binder can clip two chosen hotspots while sitting well off its designed pose and epitope recall will pass it. Run the pose gate; do not let hotspot recall stand in for it. The route is `usalign --mm 1`, predicted complex against the rfd3 backbone complex, `--col-a boltz_<label>_path --col-b rfdiffusion3_path` (lineage-resolved). Both hold the same target, so the superposition locks onto it and the binder's displacement shows in that frame. **Do not threshold on TM or RMSD: both are target-dominated and the RMSD lies.** `--mm 1` keeps only residues that structurally align, so a binder that drifted drops out and leaves a target-only alignment with a beautiful low RMSD (about 0.4 A). **Read `Lali − L_target`** (binder residues that aligned in the target frame: `_Lali` over all chains minus the fixed target length): near the binder length means it stayed, near 0 means the binder is gone. In a batch where most binders had drifted it sat near 0 to 1 against 60-90 residue binders and correlated strongly and positively with RMSD (Spearman about +0.9), so ranking by low RMSD selected the worst designs; cross-check against `epitope` recall on the prediction.

**Gate 2.1, did it fold?** Binder chain only, `usalign --mm 0` against the designed backbone (`--mm 1` on a complex is dominated by the fixed target and reads about 1.0 regardless; `--ter 2` isolates the first chain while the binder is usually last; use `chainsel` to extract). This answers only "did it fold", not "where".

**Interface numbers come last, and only on a target that is intact.** A large interface area or a high confidence on a prediction whose target has collapsed or deformed (target RMSD to its reference high) is not a result; check the target against its reference before reading any interface number.

## 5. What each metric cannot see

- **Self-consistency cannot see solubility**: a binder with a few percent charged residues and a majority of TSVIG residues can fold at pLDDT 0.96. Check composition directly (% charged, net charge, % aromatic, % beta-branched).
- **CMS and SC cannot see chemistry**; they are geometry only.
- **`buried_unsat` is mostly not the interface** (typical totals 25-44 against `if_delta_unsat` 4-16); composition problems show up as packing quality (`sc` 0.36-0.53, `packstat` < 0.6).
- **Whole-complex confidence cannot see the binder** (see `boltz`): `iptm` averages N(N-1)/2 pairs, nearly all template-forced target-target.
- **The one cheap honest signal** is binder mean pLDDT alone against in complex: a real interface raises it, while a confidently folded domain the predictor is unsure where to put loses several points up to 20 or more.

## 6. Comparisons across rounds

- **Use one reference model for every pose comparison across rounds, and say which** (the backbone it was designed on, or one named prediction). Comparing round n to round n-1 lets drift accumulate unseen.
- **Relax before reading side-chain-dependent metrics** (pKa, contacts, tautomers, hydrogen-bond roles), and never mix raw and relaxed values in one ranking: split the table by provenance first.
- **A rank across runs needs a key unique across runs**: names repeat between run_dirs, so key on (run, name).

## 7. Judgment

- **Keep the spend late.** Recon, trimming, epitope choice, backbones, sequences and a fold screen are cheap. Pilot the co-fold on 2 designs before scaling, and pilot the target definition itself (two definitions head to head on 2 designs, and let a measurement pick).
- **Verify chain mappings empirically**: superpose over candidate permutations; the right one separates by orders of magnitude (0.05 A against 16 A or more). Generators renumber from 1.
- **Demand the check that would have failed** for any number gating a spend, and **state expectations before results arrive**: the wrong ones are the most useful.
- **Write the campaign up** in `campaigns/`.
