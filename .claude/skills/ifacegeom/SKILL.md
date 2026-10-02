---
name: ifacegeom
description: How to run the custom ifacegeom tool — recording a binder's interface footprint on its target (the epitope residue list, as a column that survives lineage), the target-side histidine count, and where the binder's N- and C-termini sit relative to that interface. Covers the two selection methods and why the vector one is the default, the chain-selection trap, the residue-label format every downstream tool parses, and every column it collects. Load before composing an ifacegeom run or interpreting its columns.
---

# ifacegeom

**Custom tool** (lives in `tools/ifacegeom/`, not bundled with prosapia).
**`action: update`** — annotates the table it reads, in place. `-t` is required.
**CPU-only**; the manifest builder forces `gpus_per_task = 0`, so you do not pass `-g 0`.
**No default input column** — `-i` is effectively required (see Traps).

Measures, for each binder/target **complex**:

1. **Which residues form the interface, on both sides.** The binder-side list is the
   *epitope footprint*; the target-side list carries the histidine count.
2. **Where the termini sit relative to that interface** — a signed projection onto the
   binder-COM → interface-COM axis, plus two distances per terminus.

## Why it exists

When the binders are inherited rather than designed here, the epitope is not known —
it has to be discovered before anything downstream can respect it. And every later
step needs it:

- the rotamer graft needs to know *which* binder residues to transplant,
- the H-bond network design needs to know which positions it must not overlap,
- the C-terminal tag requirement needs to know where the C-terminus is.

The reason this is a **tool** and not a side script is the word *column*. An epitope
list in a sidecar file does not survive the generation boundary into the dock and
network tables, which is exactly where it is needed. As columns, the lists ride
lineage for free and can be re-filtered at any cutoff later without recomputation.

`cms` computes per-residue interface contributions but writes them to a **file**
(`cms_path`), so they cannot be filtered on or carried through lineage.
`pyrosetta --interface` gives scalar `if_*` only — no residue list, no terminus
geometry. Neither replaces this.

## Invocation

```bash
sapia run ifacegeom <run_dir> \
    -t table0 \
    -i input_path \
    --binder-chains A --target-chains auto

sapia collect ifacegeom <run_dir> -t table0
```

Every input file must hold **both** binder and target — there is no interface in a
monomer. Run it on whichever column holds the complex: the ingested BindCraft2
output, or a `boltz_path` to measure a prediction's interface instead of the designed
pose (label the second with `-l` so the columns do not collide).

### Flags

| Flag | Default | What it does |
| --- | --- | --- |
| `--binder-chains` | `A` | Binder chain IDs, comma-joined. |
| `--target-chains` | `auto` | Target chain IDs. `auto` = every protein chain that is not a binder chain — right for a file holding one binder and one target of however many chains. |
| `--method` | `vector` | `vector` or `heavy` — see below. |
| `--contact-cutoff` | `5.5` | Heavy-atom distance (Å) for the contact test. |
| `--cb-dist-cut` | `11.0` | CB–CB distance beyond which a pair is not considered at all. |
| `--vector-dist-cut` | `9.0` | CB–CB distance within which the pointing test applies. |
| `--vector-angle-cut` | `75.0` | Max angle (°) between CA→CB and the direction to the partner's CB. |
| `--designs-per-task` | `200` | Designs per task. The work is milliseconds per design; this exists to amortise the container cold start. |

### The two methods

`vector` (default) reproduces the **idea** of Rosetta's
`InterGroupInterfaceByVector`. A residue is at the interface either because it
touches the other side, or because its side chain *points at* it:

```
nearby   any heavy atom of i within --contact-cutoff of any heavy atom of j
         -> both i and j selected
vector   CB(i)-CB(j) within --vector-dist-cut AND the angle between i's CA->CB
         direction and the CB(i)->CB(j) direction is below --vector-angle-cut
         -> i selected (and symmetrically for j)
```

The pointing test is applied **per residue, not per pair**, so each side's list
stands on its own: binder residue *i* is listed because *i* points at the target,
whether or not the target residue points back.

`heavy` keeps the contact test alone — the plain heavy-atom contact definition.

The vector test is what excludes residues that are merely nearby while facing away.
Those are exactly the rotamers that would be wasted work to graft back later, so
`vector` is the right default for this campaign. `vector` is always a **superset** of
`heavy` (verified: 24 vs 21 binder residues on barnase/barstar).

This is a reimplementation of the concept, **not a port of Rosetta's code** — there
is no PyRosetta in this image and the numbers are not bit-identical. The defaults are
Rosetta's documented ones.

## Traps

- **`-i` is effectively required.** `default_input_column` is the `"not applicable"`
  sentinel (as in `cms`, `usalign`, `chainsel`): there is no honest default complex
  column, and a wrong one fails silently, so the builder refuses the run and lists
  the `_path` columns the table actually has.
- **A chain named but absent is an error for that design, never a smaller
  selection.** `--binder-chains A,B` against a file holding only `A` fails that row
  rather than quietly measuring a half binder. Check `ifacegeom_binder_chains` and
  `ifacegeom_target_chains` after collect — with `--target-chains auto` they tell you
  what was actually used.
- **`warn: no interface residues` almost always means the wrong chains**, not a real
  non-contact. The row keeps `binder_len`/`target_len`/the chain names (which is what
  makes the mistake obvious) and every geometry column is NA, because without an
  interface there is no axis. The `warn:` prefix keeps it out of any `== "OK"`
  selection.
- **`--contact-cutoff` is an ATOM-to-atom distance.** The ~8–11 Å numbers usually
  quoted alongside `InterGroupInterfaceByVector` are its **CB–CB** cutoffs
  (`--vector-dist-cut`, `--cb-dist-cut`), not this one. Setting `--contact-cutoff 8`
  does not "match Rosetta" — it selects most of the binder.
- **The COMs are in each file's own frame.** `ifacegeom_binder_com` /
  `_iface_com` are recorded for auditability, not for comparison across designs.
  Anything downstream that works in a different frame (a dock, a superposition) must
  recompute them there.
- **Residue labels are `chain:resnum[icode]`** — `A:12`, not `A12`. This is the
  format every downstream consumer of `binder_res` parses. Note `ringfit --hotspots`
  uses the *other* convention (`A47`); do not paste one into the other.
- **The numbering is this file's numbering.** The epitope list is only meaningful
  against a structure with the same residue numbers. It has to survive three later
  renumberings (chainsel's extraction, rpxdock's forced 1..N-per-chain dump,
  hbdesigner's offset), and this tool does not do that remapping — the consumer does.

## Columns collected (`ifacegeom_` prefix)

`<leaf>_status` is `OK`, `warn: no interface residues`, `error: …` or `missing`.
`<leaf>_path` is the per-residue TSV (side, label, chain, resnum, resname, selected,
selected_by, min_heavy_dist, min_cb_dist) — every residue near the other side,
selected or not: the evidence behind the lists.

| Column | Meaning |
| --- | --- |
| `method` | `vector` or `heavy` — which criterion produced the lists |
| `binder_chains`, `target_chains` | the two sides actually used (`auto` is resolved here) |
| `binder_len`, `target_len` | amino acids per side |
| `binder_res` | **the epitope**: binder interface residues as `A:12,A:15,…` |
| `n_binder_res` | |
| `binder_res_seq` | their one-letter codes, in the same order — lets a later graft be checked by residue identity |
| `binder_his`, `n_binder_his` | histidines already in the binder's own epitope |
| `target_res`, `n_target_res` | target interface residues, same format |
| `target_his`, `n_target_his` | **histidines in the target epitope** — the residues whose contacts change when the pH drops |
| `binder_com` | mass-weighted centre of the whole binder, `x,y,z` |
| `binder_iface_com` | centre of the binder-side interface residues = the **epitope COM** |
| `iface_com` | centre of both sides' interface residues = the axis endpoint |
| `axis_len` | `|iface_com − binder_com|` |
| `nterm_res`, `cterm_res` | which residues the termini are (`A:1` / `A:115`), so every projection is checkable by hand |
| `nterm_proj`, `cterm_proj` | **signed** projection (Å) of the terminus CA on the `binder_com → iface_com` axis, from `binder_com`: **> 0 = interface side** of the plane, **< 0 = far side** |
| `nterm_iface_dist`, `cterm_iface_dist` | terminus CA to `iface_com` (Å) |
| `nterm_iface_min_dist`, `cterm_iface_min_dist` | terminus CA to the **nearest epitope CA** (Å) — the number that actually answers "would the tag sit in the epitope". `0.0` means the terminus *is* an interface residue. |
| `seconds` | wall time of the design's measurement |

NA (empty) rather than 0 marks a field that did not apply: the design failed, or the
geometry was undefined.

## Reading the numbers

**The C-terminus requirement is two-sided, and neither side is filtered on here.**
The C-terminus must stay clear of the target epitope (the GFP/split-strep tag would
otherwise sit in the binding site) *while still being reachable* from the partner
monomer's N-terminus across the dimerization interface. So `cterm_proj` wants to be
negative, but not so negative that the linker cannot reach. Both numbers are
**recorded**; the thresholds are the user's call, from the distribution.

`cterm_iface_min_dist` is the sharper of the two C-terminus numbers —
`cterm_iface_dist` is to the interface *centroid* and so grows with binder size,
while the min-distance is to the nearest epitope residue and means the same thing on
every design.

`n_target_his` is a **record, not a gate**. A target epitope with histidines is one
whose contacts change with pH, which the strategy would rather avoid — but the pool
is inherited and finite, so exclude nothing before seeing the distribution.

## The inherited hEGFR pool — facts you need before running it there

Measured on `sapia-runs-toon:inputs/bc2_output_for_anthony/` (run
`outputs/20261002_143419_dimer_phase2`):

- **The binder is chain B; the target is chain A.** The opposite of this tool's
  default. Pass `--binder-chains B --target-chains A` or every number is wrong.
  The filename's `l<N>` is the binder length, and it matches `ifacegeom_binder_len`
  on all 37 designs.
- **Chain A is human EGFR 311–503, renumbered 1–193.** The mapping is
  `local resnum = human EGFR resnum − 310`, verified by the residue identity of all
  five design-document hotspots (L325→A:15 L, P349→A:39 P, F412→A:102 F,
  V417→A:107 V, I467→A:157 I). `ifacegeom_target_res` is in **local** numbering;
  add 310 before comparing to anything written in EGFR numbering.
- **Use the top-level `*_hEGFR.cif`, not `renum/`.** `renum/` differs only in that
  the binder chain continues the target's numbering (B:194… rather than B:1…); the
  target numbering is identical, and per-chain numbering is what `chainsel`/`rpxdock`
  expect downstream.
- Each design has a paired `*_mEGFR.cif`. Measuring the mouse complexes is a second
  `ifacegeom` run — seed them as their own rows, or use `-l` on a second table.

## Verification done

On barnase/barstar (1BRS, chains D=binder / A=target), checked against the
literature and against brute-force computation written independently of the tool:

- `target_his` = `A:102` — barnase **His102**, the known interface histidine, and the
  only one of barnase's two at the interface.
- `binder_res` contains barstar Tyr29, Asp35, Asp39, Trp44, Glu76 — the classic hot
  spots; `binder_res_seq` is in register at all of them.
- `--method heavy` reproduces a brute-force 5.5 Å heavy-atom contact set **exactly**
  (21/21 binder, 25/25 target residues).
- `binder_len` / `target_len` match an independent residue count.
- On a synthetic two-chain case with a known answer, the axis is exactly `+x`, the
  C-terminus projects positive and the N-terminus negative.
- Absent chain, same chain on both sides, and a missing file each give an `error:`
  status with every metric NA, and the task still exits 0.

## Verified on real data (Modal)

`outputs/20261002_143419_dimer_phase2` on `sapia-runs-toon`, 37 inherited
BindCraft2 hEGFR complexes, one batched task, **37/37 `OK`**, empty stderr:

> **Cite a run_dir that still exists.** This section previously cited
> `outputs/20261002_140341_ifacegeom_hegfr_test`, which was deleted before the
> campaign proper began. The numbers were genuine and reproduced exactly in the
> run above — but for a while this page asserted results nobody could re-read,
> and a campaign decision was nearly taken on them. **A commit message is not a
> measurement.** If the run backing this section is ever cleaned up, re-point
> this section or delete the claim.

- `ifacegeom_binder_len` matched the `l<N>` parsed from every filename (37/37) — an
  invariant the tool did not compute.
- **All five design-document hotspots are contacted by 37/37 binders**
  (L325, P349, F412, V417, I467). The epitope the tool finds is exactly the one
  BindCraft2 was aimed at — independent evidence that the selection is right.
- The whole pool contacts **His409** (37/37) and most of it **His346** (29/37), so
  **0/37 have a histidine-free target epitope**. His409 sits two residues from the
  F412 hotspot, i.e. the intended hydrophobic patch has a histidine built into it.
  This is the §1.1 "record, don't filter" principle earning its keep: filtering on
  histidine-free up front would have emptied the pool at step 1.
