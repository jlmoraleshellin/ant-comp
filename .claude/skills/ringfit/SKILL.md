---
name: ringfit
description: >-
  Scores whether a binder straddles two protomers of an oligomeric target and clashes with the rest of the assembly or its lipid belt. Load before composing a ringfit run or interpreting its columns. Traps: the target/binder chain mapping, hotspot numbering and resnum matching.
---

# ringfit

**Custom tool** (lives in `tools/ringfit/`, not bundled with prosapia).
**`action: update`** — annotates the table it reads, in place. `-t` is required.
**CPU-only**; the manifest builder forces `gpus_per_task = 0`, so you do not pass `-g 0`.

Superimposes each design's *target* chains onto the matching chains of a reference
full assembly, then asks two separate questions:

1. **Does the binder straddle the seam?** buried surface area per target chain,
   `bridge_ratio`, contacted residues per protomer, hotspot recall.
2. **Would it survive in the real assembly?** clashes against the reference chains
   that were *not* part of the design, and against the lipid belt — surface that is
   artificially exposed whenever you design against an isolated protomer pair.

## Why it exists

Designing against two protomers cut out of a ring leaves two lies in the model: the
faces where the neighbouring protomers used to be, and the lipid-facing belt. A binder
is free to occupy either. ringfit puts the design back into the full assembly and
measures how much of that freedom it used.

## It needs exactly two target chains

Its premise is **two adjacent protomers cut from a larger assembly**, plus a reference
structure to put them back into. **`--target-chains` must name exactly two chains**: one
is a hard error (`ringfit_worker.py:494-498`, `--target-chains needs exactly 2 chains`).

**Never work around that with `--target-chains A,A`.** It passes every guard and scores
chain A against itself, so `bsa_t1 == bsa_t2` by construction and
**`bridge_ratio` comes back `1.0` for every design** — a perfect straddle of a seam that
does not exist. The signature is a whole table of `ringfit_bridge_ratio == 1.0`.

**A single protomer cut out of an oligomer is a real question, but ringfit is not the
tool for it.** When the design has one target chain but the reference has neighbours
(chain A of a homotrimer, say), the failure mode `n_clash` detects **can** occur: the
binder can sit on surface the other protomers bury. That is a different measurement
from ringfit's two-protomer one, and bending ringfit into it only manufactures a
number.

**"Did my binder hit the epitope I chose, on any target" is the `epitope` tool** — same
hotspot idea, no assembly premise, and a hotspot the target does not have is a hard
error. It is **not** a substitute for the neighbour-space question: it knows nothing
about chains that are absent from the design file. A matching `default_input_column`
(`rfdiffusion3_path`) is not a reason to reuse this one.

## Invocation

```bash
sapia run ringfit <run_dir> \
    --table table1 \
    --input-column boltz_path \
    --ref-structure inputs/7ojg_assembly.cif \
    --target-chains A,B --ref-target-chains A,B \
    --hotspots A35,B30 --hotspots-in-ref-numbering
sapia collect ringfit <run_dir> --table table1
```

| Flag | Default | Notes |
| --- | --- | --- |
| `--ref-structure` | **required** | Full assembly, PDB or mmCIF. **Must be a path the workstation can see under `/runs`** — a file under `<run_dir>/inputs/` is on the runs volume and is fine. The builder `raise`s if it cannot see it, and `volume_path()` rewrites it for the task container. |
| `--input-column` | `rfdiffusion3_path` | The structure to score. `boltz_path` to score predictions, `rfdiffusion3_path` to score the designed backbone. Nothing assumes a column. |
| `--target-chains` | `A,B` | The **design's** target chains, ordered — first is `t1`, second is `t2` in `bsa_t*` / `n_contact_res_t*`. |
| `--ref-target-chains` | `A,B` | The **reference** chains they superimpose onto, same order. Every other reference protein chain becomes the clash check. |
| `--binder-chains` | `auto` | `auto` = every protein chain that is not a target chain. |
| `--hotspots` | `""` | e.g. `A47,B42`. Empty → `hotspot_recall`/`hotspot_hits` collected as NA. |
| `--hotspots-in-ref-numbering` | off | **See the trap below.** |
| `--resnum-match` | `auto` | `ordinal` \| `resnum` \| `auto`. See below. |
| `--clash-cutoff` | `2.5` | Å, heavy-atom, binder vs. rest of assembly. |
| `--contact-cutoff` | `5.0` | Å, heavy-atom, defines a contacted target residue. |
| `--lipid-resnames` | `PLM,LPP,L8Z` | Reference HETATM comp-ids marking the belt. Empty disables the check (NA). |

## Traps

**Hotspot numbering.** Generators such as rfdiffusion3 renumber every chain from 1, so
reference residue `A35` may be design residue `A18`. A hotspot list copied off the
reference and passed *without* `--hotspots-in-ref-numbering` silently scores the wrong
residues — you get a meaningless `hotspot_recall`, not an error. Rule: **if you read the
hotspots off the reference structure, set the flag.** Then check
`ringfit_resnum_offset` — `0` means the numberings agreed and the flag changed nothing.

**Residue matching drives every number.** `ordinal` pairs the i-th design residue with
the i-th reference residue (right after renumbering); `resnum` pairs equal residue
numbers (right when numbering was preserved, and gap-tolerant); `auto` picks ordinal
when the CA counts match, resnum otherwise. Getting this wrong does not error — it
produces a plausible superposition of the wrong residues. **Always read
`ringfit_seq_match_frac` before trusting anything else:** `1.0` = correct pairing,
low = register shift or wrong chains. Below `0.9` the status becomes `warn: ...`, which
keeps the row out of a `== "OK"` selection while still reporting the numbers so you can
judge the suspicion yourself.

**`align_rmsd` is the trust metric, not a result.** It is the CA RMSD of the
target-chain superposition. If it is large, the design's targets did not match the
reference and every ring/lipid/BSA number downstream is meaningless.

**Lipid comp-ids are structure-specific.** The default `PLM,LPP,L8Z` is a guess. Check
the actual HETATM comp-ids in your reference and pass them, or the lipid check silently
measures nothing and collects NA.

## Columns collected (`ringfit_` prefix)

**Trust first:** `align_rmsd`, `seq_match_frac`, `resnum_offset`.

**Straddle:** `bsa_t1`, `bsa_t2`, `bsa_total` (Å², per target chain and combined,
t1/t2 in `--target-chains` order); `bridge_ratio` = `min(bsa_t1,bsa_t2)/max(...)`,
so `1.0` is an even straddle of the seam and `~0` is a one-protomer binder;
`n_contact_res_t1`, `n_contact_res_t2`; `hotspot_recall`, `hotspot_hits`.

**Assembly viability:** `n_clash`, `n_clash_res`, `min_dist_ring` (binder vs. the
reference's other protein chains); `lipid_clash`, `min_dist_lipid` (vs. the belt).

**Bookkeeping:** `binder_len`, `binder_chains`. Plus `ringfit_status` and
`ringfit_path` (the design superposed into the reference frame — load this in a viewer
to see what actually happened).

**NA (empty) means not applicable, not zero** — no hotspots given, no lipid resnames,
or no reference chains/lipids to measure against. Do not read a blank as a pass.

## Forking

It is an `update` tool, so run it twice with different `--dir-label` / `--table-label`
to score two different structure columns (e.g. the designed backbone and the Boltz
prediction) onto the same table and compare leaf-keyed columns side by side. That
comparison — designed BSA/hotspot recall vs. predicted — is the "did the pose hold"
readout.
