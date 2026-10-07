# Recon: human vs mouse TNF-alpha, cross-reactive epitope selection

**Date:** 2026-10-07
**Goal:** binders that bind both human and mouse TNF-alpha.
**Inputs:** `projects/TNFalpha_binders/1TNF.cif` (human), `2TNF.cif` (mouse).
**Method:** gemmi; Shrake-Rupley SASA (probe 1.4 A, Tien 2013 reference maxima);
CA superposition; all distances are heavy-atom minima.

## 1. The two structures

| | 1TNF | 2TNF |
| --- | --- | --- |
| organism (from file metadata) | *Homo sapiens* | *Mus musculus* |
| method / resolution | X-ray 2.6 A | X-ray 1.40 A |
| chains | A/B/C homotrimer | A/B/C homotrimer |
| modelled residues (chain A) | 152 (6-157) | 148 (9-157) |
| SEQRES length | **157** | **156** |

Both verified from `_entity_src_gen` / `_struct.title`, not from the filenames.

**Backbone superposition (chain A, 148 CA): RMSD 1.11 A.** The two proteins are
structurally near-identical. Cross-reactivity is a *sequence* problem here, not a
fold problem.

## 2. Sequence relationship

**80.4% identity** over the 148 shared modelled positions (29 substitutions).

The PDB numbering is already aligned between species: human 9 = mouse 9, etc.
Identity computed by residue number is 80.4%, which would be ~5% if the numbering
were offset, so the alignment holds throughout.

### Residue 73 is deleted in mouse

Mouse SEQRES is 156, human is 157. Mouse chain A jumps 72 -> 74. This is a genuine
one-residue deletion, **not** disorder. Human 73 is His.

> **Numbering trap.** Because of this deletion, position *n* in the real mouse
> sequence does not equal position *n* in the mouse PDB for n > 73. The 2TNF file
> preserves human-aligned numbering by skipping 73. Any hotspot list handed to a
> tool must state which numbering it is in.

### Substitutions (human -> mouse, by residue number)

```
 20 P>H    22 A>V    24 G>E    27 Q>E    30 N>S    31 R>Q
 41 V>M    42 E>D    44 R>K    52 S>A    53 E>D    58 I>V
 71 S>D    72 T>Y    80 I>V    83 I>F    85 V>I    89 T>E
 97 I>V   102 Q>P   103 R>K   104 E>D   111 A>L   131 R>Q
136 I>V   138 R>L   140 D>K   143 L>D   154 I>V
```

## 3. The colleagues' candidate patches, scored

Relative SASA computed in the **trimer** context (chain A surrounded by B and C).
"Exposed" = rel-SASA >= 0.25. A buried difference does not affect binding, so the
column that matters is the last one.

| patch | size | diffs | exposed | **exposed + different** | verdict |
| --- | --- | --- | --- | --- | --- |
| 65-67 | 3 | 0 | 2 | **0** | clean |
| 143-149 | 7 | 1 | 3 | **0** | clean |
| 31-33 | 3 | 1 | 3 | **1** (R31Q) | usable |
| 86-92 | 7 | 1 | 6 | **1** (T89E) | usable |
| 72-79 | 8 | 2 | 1 | **1** (T72Y) | mostly buried, poor epitope |
| 20-25 | 6 | 3 | 5 | **2** (P20H, G24E) | risky |
| 135-140 | 6 | 3 | 4 | **2** (R138L, D140K) | risky |

Notes:

- **72-79 is almost entirely buried** (only residue 72 is exposed). It cannot serve
  as an epitope regardless of its conservation.
- **135-140's two differences are charge reversals** (R138L loses a positive charge,  D140K reverses it). Charge reversals at an interface are the substitutions most
  likely to break cross-reactivity.
- **20-25** carries P20H and G24E, both exposed. G24E also adds a charge.

## 4. Spatial layout: which patches form one surface

A binder lands on a contiguous surface, not on a list of ranges. Minimum
heavy-atom distance between the *exposed* residues of each patch, chain A:

```
             20-25    31-33    65-67    72-79    86-92  135-140  143-149
    20-25        -      3.6      9.0     20.7     12.5      5.0      3.6
    31-33      3.6        -     13.8     29.7     21.5     15.3      2.8
    65-67      9.0     13.8        -     11.4     21.5      4.1      5.3
    72-79     20.7     29.7     11.4        -     24.7      8.8     22.1
    86-92     12.5     21.5     21.5     24.7        -      3.1     18.3
  135-140      5.0     15.3      4.1      8.8      3.1        -      7.5
  143-149      3.6      2.8      5.3     22.1     18.3      7.5        -
```

Two clusters fall out:

- **Cluster 1 — 31-33 / 143-149 / 65-67 / 20-25.** Mutually 2.8-9.0 A apart: one
  continuous surface. Contains both "clean" patches.
- **Cluster 2 — 86-92 / 135-140 / 72-79.** Also contiguous (3.1-8.8 A), but built
  from the two risky patches plus the buried one.

### Which protomer interface each patch sits on

| patch | -> chain B | -> chain C | |
| --- | --- | --- | --- |
| 20-25 | 20.1 | 13.7 | single-protomer face |
| 31-33 | 21.1 | **5.8** | A/C interface |
| 65-67 | 14.1 | **6.4** | A/C interface |
| 72-79 | **3.0** | 10.6 | A/B interface |
| 86-92 | **3.2** | 20.7 | A/B interface |
| 135-140 | 9.6 | 13.0 | single-protomer face |
| 143-149 | 18.0 | **5.1** | A/C interface |

Cluster 1 is a **composite epitope straddling the A/C seam** — the same geometry as
the receptor-binding groove, and the same situation as the earlier SlyB campaign.

## 5. Recommended epitope

**Cluster 1, A/C seam, trimmed to its conserved core.**

Exposed *and* conserved residues available there:

- chain A: **32, 33, 65, 67, 113, 144, 145, 147**
- avoid: A31 (R>Q), A20 (P>H), A24 (G>E)

Chain C contributes the other half of the seam. Within 10 A of the chain-A core,
chain C is largely conserved — C82, C87, C91-C96, C123-C126 all identical — with
four exceptions to route around:

| chain C residue | change | distance to core |
| --- | --- | --- |
| C97 | I>V | 3.7 A (closest) |
| C102 | Q>P | 6.4 A |
| C73 | **deleted in mouse** | 7.9 A |
| C80 | I>V | 9.5 A |

I>V is conservative (both small hydrophobics) and both sit at the rim. **Q102P and
the C73 deletion are the two to keep outside the footprint.**

### Open issue: the footprint is currently too wide

The chain-A core as listed spans **32.6 A** at its widest. A 50-120 residue binder
covers roughly 20-30 A. The epitope needs trimming to a ~25 A core before it is
handed to a generator — otherwise the binder will cover only part of it, and which
part is left to chance.

## 6. LPC's built-in TNF target cannot be used as shipped

`configs/targets/targets_dict.yaml`, target `38_TNFalpha`:

```yaml
    target_input: "A12-157,B12-157,C12-157"
    hotspot_residues: ["A113", "C73"]
    binder_length: [50, 120]
```

- **A113 is fine** — conserved (Pro in both), and 2.9 A from chain C, so genuinely
  at the seam.
- **C73 is deleted in mouse.** A binder steered onto C73 is being steered onto a
  residue that does not exist in the mouse protein.

So the target entry is reusable but **the hotspot list must be replaced**. That is
a one-line config change, and it is the single most important decision in the
campaign.

Separately: do **not** use the `38_TNFalpha_FIX` variant. Its PDB
(`1tnf_cropped_fixed.pdb`) is committed as a 0-byte file upstream, with the real
content in an accidental editor backup (`1tnf_cropped_fixed.pdb~`) beside it. Use
`38_TNFalpha` or `38_TNFalpha_REPACK`.

## 7. What is not yet done

- Trim the epitope to a ~25 A footprint and pick the final hotspot set.
- Build the mouse target entry. 2TNF is a clean 1.40 A trimer and needs its own
  `target_input` (mouse chain A is modelled 9-157, so `A9-157,B9-157,C9-157`) plus
  hotspots **in mouse numbering**.
- Decide the binder length range.

## 8. Reproducing this

Scripts are in the session scratchpad and should be copied here if this is to be
rerun:

```
campaigns/tnf/scripts/
  inspect_cif.py  file metadata, organism, chains
  seqs.py      extract chain A sequences
  align.py     identity + substitution list
  recon.py     SASA + per-patch conservation table
  geom.py      patch adjacency + interface proximity
```

Run with the repo venv, e.g.:

```bash
cd ~/ant-comp
.venv/bin/python -I campaigns/tnf/scripts/recon.py \
    ~/projects/TNFalpha_binders/1TNF.cif ~/projects/TNFalpha_binders/2TNF.cif
```
