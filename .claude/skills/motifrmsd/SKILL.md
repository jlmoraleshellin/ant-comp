---
name: motifrmsd
description: >-
  Did rfdiffusion3 rebuild a grafted motif in its reference 3D geometry? Takes an rfd3
  motif-scaffolding backbone column and writes motif CA RMSD after a Kabsch fit ON THE MOTIF,
  plus a per-segment RMSD from that one shared fit. Load before composing a motifrmsd run, and
  before gating a motif-graft campaign on any number: it reads rfd3's OWN per-design
  diffused_index_map sidecar, so it works on rfd3 motif output and errors on everything else —
  it is not a general two-structure comparator (that is usalign). Also writes
  design_fixed_positions — the SAME motif segments in the DESIGN's own output numbering,
  ready to paste into a downstream ProteinMPNN --fixed-positions so the graft is held fixed
  while the rest of the binder is redesigned. In a BINDER-SCAFFOLDING contig (a target held
  fixed plus a graft on the binder chain) rfd3 maps the whole target too, so that column is
  always NA unless you pass --motif-chain <reference chain letter> — read that section before
  any such run. Covers the --ref-structure vs specification.input policy, why a count mismatch
  is an error rather than a smaller RMSD, and why seq_match_frac must be read before the RMSD.
---

# motifrmsd

**Custom tool** (lives in `tools/motifrmsd/`, not bundled with prosapia).
**`action: update`** — annotates the table it reads, in place. `-t` is required.
**CPU-only**; the manifest builder forces `gpus_per_task = 0`, so you do not pass `-g 0`.

For each rfdiffusion3 motif-scaffolding design: read rfd3's own per-design sidecar JSON,
take the residue correspondence it recorded (`diffused_index_map`), lift the CA atom of
every mapped pair out of the reference and the design, **superpose the design's motif CA
onto the reference's motif CA (SVD/Kabsch)**, and report the RMSD of those same atoms
after that fit. That is the literal "low motif RMSD after alignment on the motif" gate.

## Premise, and when its numbers are meaningless

**The premise is: the motif correspondence is supplied externally, by rfd3, per design.**
The tool never derives it. It does not parse the contig, and it does not read
`metrics.max_ca_deviation` (a maximum, not an RMSD, in a reference frame nothing in this
repo can confirm).

That makes the scope narrow and the refusals clean:

| Input | What happens |
| --- | --- |
| rfd3 **motif-scaffolding** output (a sidecar with `diffused_index_map`) | the measurement works |
| rfd3 **de-novo** output, Boltz/AF3 predictions, ProteinMPNN children, any non-rfd3 structure | **`error: rfd3 sidecar JSON could not be located…`** or `error: … has no 'diffused_index_map'`. There is **no** fallback to a whole-structure alignment |
| a design whose mapped residues are not all present in both files | **`error: matched N of M motif residues: …`**, every numeric column NA |
| a `diffused_index_map` that is not one-to-one (two reference residues sent to one design residue, or two key spellings of one reference residue) | **`error: diffused_index_map in … is not one-to-one: …`**, every numeric column NA, **even with `--allow-partial`** |
| the rfd3 `.cif_to_pdb` cache pruned from the run_dir | **`error:`** — the `.src` record is how the path column is traced back to the sidecar |

**Unlike most tools here, it does not return a number anyway.** Every path that cannot
resolve the correspondence writes an `error:` status and NA, never a plausible RMSD over a
smaller or different atom set. That is deliberate: this tool exists because the column it
replaces (`rfdiffusion3_rfd3_ca_rmsd_to_input`) was a passthrough of a key that does not
exist, collected as `OK`.

**What it is NOT scoped to:** comparing two structures that have no externally supplied
correspondence. That is `usalign`, which finds its own alignment — a different premise and a
different question.

**"The motif" is every residue in the design's `diffused_index_map` — full stop.** The tool
does not know which of those residues are the loops you care about. On this campaign's
current contigs (motif segments only, no target) the map is exactly the grafted CDRs. On a
**target-bearing or partial-diffusion rfd3 run** the map can also list fixed context or
target residues; those sit at ~0 deviation by construction, so they **dilute
`motifrmsd_rmsd` toward zero** and add their own entries to `n_segments`, `seg_ranges` and
`seg<i>_rmsd`. Before using this tool on such a run, check `seg_ranges` and `n_expected`
against the residues you meant: if they include anything that is not the motif, the RMSD
columns answer a different question. (Scope caveat for future runs, not a current bug.)

## Trust columns — read these before any RMSD

| Column | Expected | How to read it |
| --- | --- | --- |
| `motifrmsd_status` | `OK` | **The only proof of success.** `warn:` = numbers computed but the correspondence or the reference is suspect. `error:` = nothing computed, every number NA. `missing` = no output file |
| `motifrmsd_n_expected` | the motif size, e.g. 34 | residues listed in `diffused_index_map` |
| `motifrmsd_n_matched` | **equal to `n_expected` on every `OK` row** | pairs whose CA was found in both files. A shortfall is an `error` by default; with `--allow-partial` it is a `warn:` row with the counts in both columns |
| `motifrmsd_seq_match_frac` | **unknown until real data — see below** | fraction of matched pairs whose residue identity agrees between reference and design |
| `motifrmsd_ref_source` | `sidecar` or `flag` | which reference was measured against |
| `motifrmsd_ref_matches_spec` | `1` | `1` = byte-identical to what rfd3 diffused. `0` = a different file (status is `warn:`). **`NA` = the question is open, never read it as 1** — and **`NA` does not downgrade `status`**: with `--ref-structure` given and the sidecar's `specification.input` unreadable from the task, the row is `OK` with `ref_matches_spec` NA. That is deliberate (the flag's file is what was measured; the failed check is a secondary cross-check against rfd3's provenance record), but it means **`status == "OK"` does not imply agreement with what rfd3 diffused**. If you need that, filter on `ref_matches_spec == 1` explicitly |
| `motifrmsd_ref_used` / `motifrmsd_spec_input` | — | the file actually read, and the one rfd3 recorded. The pair is what makes a `ref_matches_spec` of 0 diagnosable |
| `motifrmsd_ref_chains` / `motifrmsd_design_chains` | e.g. `D` / `A` | chains the map addressed. **rfd3 renumbers every output chain from 1**, so these two are routinely different letters over different number ranges — that is normal here, and it is exactly why the map has to be read rather than guessed. **Multi-letter (`A,B,C,D`) is normal and correct on a binder-scaffolding run**, where rfd3 maps the held-fixed target as well as the graft — it is not a bug, but it does mean `design_fixed_positions` needs `--motif-chain` |
| `motifrmsd_fixed_positions_chain` | e.g. `D`, or NA | **read before pasting `design_fixed_positions` anywhere.** The design-side chain letter that spec addresses — the `--chains-to-design` value. NA exactly when the spec is NA. Never substitute `design_chains` for it |
| `motifrmsd_ref_motif_seq` / `motifrmsd_design_motif_seq` | — | the one-letter motif sequences, in segment order. Read both whenever `seq_match_frac` is not 1 |
| `motifrmsd_sampled_contig` | — | the resolved contig rfd3 built for this design (lengths are drawn per batch, so this is genuinely per design) |

### `seq_match_frac` is reported, never gated — and its expected value is unknown

Whether rfd3 writes the motif's real residue identities into a scaffolded backbone, or a
poly-glycine, **is not established in this repo**. So the status is deliberately not gated on
this column. Read it as a pattern across the table, not per row:

- **uniform and low across every design** → it reflects rfd3's output convention (check
  `design_motif_seq`: all `G` means poly-glycine). Not a mapping failure.
- **1.0 on most rows and low on a few** → those rows' mapping is wrong. Treat them as failures
  even though the status says `OK`.

Establishing which regime holds is part of the real-data verification, and the number should
be written into the campaign log once it is known.

## The question it answers, and the threshold

**"Did rfd3 rebuild the grafted motif in the reference's 3D geometry?"** — and nothing else.

There is **no house threshold here**. The MotifBench-style convention for a solved motif is
sub-angstrom CA RMSD, but that is a benchmark convention on their motif definitions, not a
precedent measured on this repo's data. Pick a bar from the distribution of the first real
run, and record it in the campaign log with the run it came from.

Use `worst_seg_rmsd`, not `rmsd`, when the motif has more than one segment: a global RMSD of
0.9 Å can hide one loop at 2 Å averaged against three at 0.5 Å, and for a paratope graft the
bad loop is the whole point.

## Invocation and flags

```bash
sapia run motifrmsd outputs/20261008_085450_5M2J_vhh2_motif_graft \
    -t table0_seg2 \
    -i rfdiffusion3_seg2_path \
    --ref-structure /runs/inputs/5M2J_D_vhh.pdb

sapia collect motifrmsd outputs/20261008_085450_5M2J_vhh2_motif_graft -t table0_seg2
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `-i/--input-column` | **`rfdiffusion3_path`** | the backbone column. **A labelled rfd3 run writes `rfdiffusion3_<label>_path`, so this default is wrong for any run that used `--dir-label`** — pass `-i` explicitly or the run submits nothing |
| `--ref-structure` | `""` | the reference the motif was lifted from, in **its own** chain ids and numbering (the left side of `diffused_index_map`, e.g. chain D of 5M2J). Empty falls back per design to the sidecar's `specification.input` — what rfd3 actually diffused. **When given, the flag wins**, and its *content* (sha256) is compared against the sidecar's input → `ref_matches_spec`. (Modal: must be a path on the runs volume) |
| `--allow-partial` | off | measure the matched subset when some mapped residues are absent, with a `warn:` status. **Off on purpose**: a count mismatch means the correspondence did not resolve, and an RMSD over a quietly smaller atom set is the failure this tool exists to prevent. `design_fixed_positions` is **NA** on such a row — it needs a complete match |
| `--motif-chain` | `""` (unrestricted) | restrict **`design_fixed_positions` alone** to the motif residues lifted from ONE **reference** chain, e.g. `D`. **Required for any binder-scaffolding contig** — see the section below. The letter is the **reference/key** side of `diffused_index_map` (the same frame as `--ref-structure` and the contig), **not** the design's output chain. Every other column keeps measuring the whole mapping, flag or no flag. A malformed value (`D1`, `A,B`, `D24-31`) raises at submit; a well-formed letter absent from a design's map NAs that one design's column with a reason |
| `-l/--dir-label` | `""` | suffix for the output dir and every column, for same-tool variants on one table |

Resources: 2 CPU, 4 GB, 20 min per task; the fit is milliseconds, so the cold start dominates.
**(Modal only)** the image is `debian_slim + gemmi + numpy`, builds in seconds.
**(vib)** `SAPIA_ACTIVATE_MOTIFRMSD` is wired to `activation/vib/prosapia_python.sh`, the same
script `ringfit`/`ifacegeom`/`dimerfit` use — gemmi and numpy come from prosapia's own deps.

## Gotchas, each with its signature

- **`Submitting 0 designs` / `No designs to submit.`** — the input column is wrong for a
  labelled rfd3 run. `rfdiffusion3_path` does not exist in `table0_seg2`;
  `rfdiffusion3_seg2_path` does.
- **Every row comes back `error: rfd3 sidecar JSON could not be located…`, and the submit log
  printed `WARNING: N of N design(s) have no resolvable rfd3 sidecar JSON`** — the input
  column is not rfd3 motif-scaffolding output (a prediction, a de-novo backbone, a child
  table's structure). The tool is refusing, correctly; it is not broken.
- **`error: … has no 'diffused_index_map' (top-level keys: [...])`** — rfd3 output, but from a
  *de novo* run with no motif. The key list in the message is verbatim from the file.
- **`error: matched 33 of 34 motif residues: 1 absent from reference …(D37)`** — the reference
  file you passed is not the one the motif was taken from, or it is a trimmed copy. **The
  counts live in the status string, not in columns**, because the numeric columns on an error
  row are all NA by design. Fix the reference rather than reaching for `--allow-partial`.
- **`warn: --ref-structure X differs in content from the sidecar's specification.input Y`** —
  you are measuring against a different file than rfd3 diffused. The numbers are kept so you
  can judge it, but the row will not pass `== "OK"`.
- **`error: … has duplicate residue keys (insertion codes or a repeated chain)`** — the
  reference is a crystal structure with insertion codes. The map's `<chain><resnum>` keys
  cannot address those unambiguously; renumber the reference first.
- **`error: diffused_index_map in … is not one-to-one: design residue(s) mapped more than once: A23 <- D24, D25`**
  — two reference residues point at the same design residue (or, for `reference residue(s)
  keyed more than once: D24 <- D24, D024`, two key spellings parse to one reference residue).
  Every count check would still have passed, and one CA would have entered the fit twice.
  This is a broken correspondence, not a missing residue, so `--allow-partial` does not
  override it. Not observed on real rfd3 output yet; if it appears, look at the sidecar.
- **`n_segments` is one higher than the number of motif ranges you wrote in the contig** —
  two ranges that are consecutive in reference numbering (`D26-32` then `D33-40`) were
  placed apart in the design, so they are scored as two loops. Correct behaviour since
  2026-10-08; before that they were merged and `worst_seg_rmsd` blended them.
- **`motifrmsd_seg2_rmsd` raises `KeyError`** — no design in this table has three segments.
  The per-segment columns exist only for indices some row actually has; `worst_seg_rmsd` is
  always there.
- **`ValueError: could not convert string to float: '5:12,30:37'` at submit time** — you
  passed `--fixed-positions "{motifrmsd_design_fixed_positions}"` to proteinmpnn. The `{…}`
  island resolves integers only; this column is a whole position spec. Substitute it
  caller-side, per design — see the hand-off section.
- **`motifrmsd_design_fixed_positions` is NA on EVERY row and `motifrmsd_design_chains` is
  `A,B,C,D`** — the signature of a binder-scaffolding contig, not of a failure. rfd3 maps the
  whole held-fixed target as well as the graft, so no single chain group can describe the
  motif. The `.err` says `the motif maps into 4 design chains (A, B, C, D); a single
  --fixed-positions chain group cannot describe that. Pass --motif-chain …`. Re-run with
  `--motif-chain <reference chain letter>` (`--force` to redo the already-`OK` rows).
  Measured on `20261009_130629_5M2J_vhh_redmotif_lowtemp_scaleup`, all 5 designs.
- **`--motif-chain A` and the column is still NA on every row**, with
  `no reference residue in chain A is part of this design's motif (diffused_index_map names
  reference chain(s) …)` — you passed a **design-frame** letter. The flag takes the
  **reference/key** side of the map; the `.err` message lists the letters that would work.
- **`ValueError: --motif-chain 'D1' is not a chain id`** at submit, nothing queued — the flag
  takes one chain id (letters only), not a residue, a range or a list. Caught before
  submission on purpose; note it still leaves a phantom registry row (see CLAUDE.md).
- **`motifrmsd_design_fixed_positions` is NA on a row whose `rmsd` is fine** — not a graft
  failure. Either the design chain is not numbered `1..N` in file order (so its residue
  numbers are not ProteinMPNN's 1-indexed chain positions), the motif spans several design
  chains (see above), or the row is a `warn:` `--allow-partial` row with
  `n_matched < n_expected`. The task's `.err` carries
  `design_fixed_positions not emitted (NA): …` with which. NA across the *whole* table with a
  single-chain `design_chains` means rfd3's output numbering is not what this column assumes.
- **ProteinMPNN redesigned the motif and the chain letter was `D`** — `D` is the *reference*
  chain from `seg_ranges`/`ref_chains`, not the design's. Take `--chains-to-design` from
  `motifrmsd_fixed_positions_chain` on that row. (On a binder-scaffolding run the design
  letter happens to be `D` too — check `fixed_positions_chain`, do not reason from the
  coincidence.)
- **The number of fragments in `design_fixed_positions` is smaller than
  `motifrmsd_n_segments`** — expected whenever `--motif-chain` is used. `n_segments` counts
  the whole map's segments; the spec counts the restricted chain's. They are different
  measurements and neither is wrong.
- **ProteinMPNN redesigned the motif anyway, at exit 0** — check for a stray trailing `/` in
  the value you passed with a single `--chains-to-design` chain. One group too many is not
  always loud; a group of the wrong *length* silently frees the residues it misses.
- **(Modal only) a half-pruned `.cif_to_pdb` cache** makes previously-fine designs error. The
  cached `.src` record is the only unambiguous route from the `_path` column to the sidecar —
  two rfd3 runs in one run_dir produce **identical file stems** in different output dirs
  (measured: `shard_0_..._model_0.cif.8ecf9c29.src` → `table0_seg2`,
  `...ac344f9c.src` → `table0_seg3`), so a stem glob would silently score against the wrong
  run's motif map. Do not "fix" a resolution error by adding a glob fallback.
- **An empty, undecodable or path-less `.src` record** (e.g. a truncated write) costs **that
  one design only**: the submit log prints `<name>: no rfd3 sidecar JSON -- cache record … is
  empty` (or `… unreadable: …`, `… holds no file path ('/')…`), the design is still submitted,
  and it collects as `error: rfd3 sidecar JSON could not be located…`. Before 2026-10-08 one
  such file raised `ValueError: PosixPath('.') has an empty name` and aborted the whole submit
  — if you see that traceback, the tool on the executor is stale.

## What it collects

Columns are leaf-prefixed (`motifrmsd_…`, or `motifrmsd_<label>_…` with `-l`).

| Column | Meaning | How to read it |
| --- | --- | --- |
| `motifrmsd_status` | **trust** | `OK` is the only proof of success. `warn:`/`error:`/`missing` as above |
| `motifrmsd_path` | the design moved into the reference frame, `<name>_onmotif.pdb` | load it next to the reference in PyMOL to see the graft. Empty on an error row |
| `motifrmsd_n_expected` | **trust** — residues in `diffused_index_map` | the motif size rfd3 recorded |
| `motifrmsd_n_matched` | **trust** — pairs actually used | equals `n_expected` on every `OK` row |
| `motifrmsd_seq_match_frac` | **trust** — identity agreement | read the pattern across the table, not the row (see above) |
| `motifrmsd_ref_source` | **trust** — `flag` \| `sidecar` | which reference won |
| `motifrmsd_ref_matches_spec` | **trust** — `1` \| `0` \| NA | NA means open, not agreement |
| `motifrmsd_ref_used` | **trust** — file read | — |
| `motifrmsd_spec_input` | **trust** — file rfd3 recorded | — |
| `motifrmsd_ref_chains` | **trust** — e.g. `D` | reference chains the map addressed |
| `motifrmsd_design_chains` | **trust** — e.g. `A` | design chains the map addressed |
| `motifrmsd_ref_motif_seq` | **trust** — one-letter, segment order | — |
| `motifrmsd_design_motif_seq` | **trust** — same positions in the design | all `G` ⇒ poly-glycine backbone |
| `motifrmsd_sampled_contig` | **trust/provenance** — resolved contig | per design; lengths are drawn per batch |
| `motifrmsd_rmsd` | CA RMSD (Å) of the motif after the fit **on the motif** | the headline gate. Lower is better. `0` is a real measurement; NA is the absence of one |
| `motifrmsd_max_dev` | largest single-CA deviation (Å) after that same fit | an RMSD hides one badly-placed residue; this does not |
| `motifrmsd_n_segments` | runs of the motif contiguous on **both** sides of the map — consecutive in the reference **and** consecutive in the design | 2 for a CDR1+CDR2 graft, 3 for CDR1+2+3. Two pieces that abut in reference numbering but were placed apart by the contig count as **two** |
| `motifrmsd_worst_seg_rmsd` | max over the per-segment RMSDs | **the filter-safe per-segment scalar** — present on every `OK` row whatever the segment count |
| `motifrmsd_seg_ranges` | e.g. `D24-37/D48-67` | segments in reference order (sorted by chain then resnum, **not** contig order). Labels are **reference** ranges, so a design-side split shows as touching ranges such as `D26-32/D33-40` — that is two segments, not a typo |
| `motifrmsd_seg_rmsds` | e.g. `0.412/0.883` | same order as `seg_ranges`; complete regardless of the column count |
| `motifrmsd_seg0_rmsd`, `_seg1_rmsd`, … | one filterable column per segment | **from the SAME single global superposition, never a per-segment refit** — a loop that is locally perfect but globally dislocated must show up, and a refit would zero it out. Exists only for indices some row has |
| `motifrmsd_design_fixed_positions` | e.g. `28:35,50:59,98:103` — **the same segments, printed from the design side**, as a ProteinMPNN `--fixed-positions` value for a **single** chain group | inclusive `start:end`, comma-joined, **no trailing slash**. Fragment *i* is segment *i*, same order as `seg_ranges`/`seg_rmsds`, so the three read side by side. See the hand-off below. **NA — never a partial string — when the motif (after any `--motif-chain` restriction) spans more than one design chain, when that chain is not numbered `1..N` in file order, or when fewer residues matched than the map names (an `--allow-partial` row).** **In a binder-scaffolding contig it is NA on every row unless `--motif-chain` is given.** Read `motifrmsd_fixed_positions_chain` with it — **that**, not `design_chains`, supplies the `--chains-to-design` letter |
| `motifrmsd_fixed_positions_chain` | **trust** — e.g. `D` | the **design-side** chain letter `design_fixed_positions` addresses, i.e. the value to pass as `--chains-to-design`. NA exactly when `design_fixed_positions` is NA. **Not `design_chains`**, which keeps listing every chain the RMSD covers (`A,B,C,D` in a binder-scaffolding run) and is therefore the wrong column to hand to ProteinMPNN |

On disk, per design: `<run_dir>/<table>/motifrmsd/<name>.tsv` (one row) and
`<name>_onmotif.pdb` (the whole design in the reference frame).

### Binder scaffolding against a held-fixed target: you need `--motif-chain`

**This is not an edge case — it is the normal shape of a binder campaign**, and the flag is
mandatory for it.

In a contig that holds a target fixed *and* grafts a motif onto the binder chain, rfd3's
`diffused_index_map` records correspondence for **every held-fixed residue** — the entire
target as well as the graft. The motif therefore lands on several design chains, the
single-chain refusal fires (correctly, per the column's own contract), and
`design_fixed_positions` is **NA on every row of the table**.

**Measured, 2026-10-09** on `outputs/20261009_130629_5M2J_vhh_redmotif_lowtemp_scaleup`:
`motifrmsd_design_chains` came back `A,B,C,D` on all 5 designs — chains A, B, C are the
three held-fixed target protomers and D is the binder carrying the CDR grafts — and
`design_fixed_positions` was NA on all 5. Nothing was wrong with the designs; the column
simply had no single chain to describe.

`--motif-chain <reference chain letter>` is the fix. It filters the pairs feeding
`design_fixed_positions` to the residues lifted from **one reference chain** before the
segment grouping runs, then applies the same checks to the restricted result:

```bash
sapia run motifrmsd <run_dir> -t table0 -i rfdiffusion3_path \
    --ref-structure /runs/inputs/5M2J.pdb \
    --motif-chain D          # D = the native VHH chain in the REFERENCE
```

Three things to hold on to:

- **The letter is in the REFERENCE frame.** It is the chain id on the **left/key** side of
  `diffused_index_map` — the same frame `--ref-structure`, `seg_ranges` and the contig use.
  It is **not** the design's output chain, which rfd3 relabels and renumbers. Getting this
  backwards gives NA with `no reference residue in chain <X> is part of this design's motif`,
  never a wrong number.
- **It changes exactly one column** (plus the letter reported in `fixed_positions_chain`).
  `rmsd`, `max_dev`, `seg_ranges`, `seg_rmsds`, `n_segments`, `seq_match_frac`,
  `design_chains` and the rest keep measuring the **whole** mapping — so the RMSD still
  covers the full held-fixed set, and `design_chains` still honestly says `A,B,C,D`.
- **The restriction is not a free pass.** The `1..N`-in-file-order check still runs, on
  whichever single design chain the restricted residues landed on, and the per-segment span
  check still runs. A restricted motif can still be NA, for the same reasons as before.

Verified by the suite: the real campaign's shape in miniature (reference chains A and B in
one map, the graft on B) gives NA with no flag, `11:14,45:47,87:91` + chain `B` with
`--motif-chain B`, and NA with a reason for a letter that is in no key.

### Handing the graft to ProteinMPNN: `design_fixed_positions`

`seg_ranges` is in **reference** numbering (`D24-37/D48-67`) and is useless as a flag.
`design_fixed_positions` is the *same* segments in the **design's own output numbering** —
the numbering `rfdiffusion3`'s backbone file itself uses. rfd3 renumbers every output chain
from 1 and draws the free-span lengths per design, so this correspondence genuinely differs
design to design; that is why it has to be a column and cannot be a constant on the command
line.

The intent is to redesign the whole binder chain **except** the graft. **Both** values
come from that design's own row — the positions from `motifrmsd_design_fixed_positions`,
and the chain letter from `motifrmsd_fixed_positions_chain`:

```bash
# the motifrmsd run, on a binder-scaffolding table:
sapia run motifrmsd <run_dir> -t table0 --ref-structure … --motif-chain D

# then, per design row:
--chains-to-design <row's motifrmsd_fixed_positions_chain> \
  --fixed-positions "<row's motifrmsd_design_fixed_positions>"
# e.g. for a row with fixed_positions_chain=D, design_fixed_positions=11:18,45:54,87:92:
--chains-to-design D --fixed-positions "11:18,45:54,87:92"
```

**The chain letter must come from `motifrmsd_fixed_positions_chain` — never assumed, never
copied from `seg_ranges`, `ref_chains` or the contig, and (since `--motif-chain` exists)
never from `motifrmsd_design_chains` either.**

- `seg_ranges` / `ref_chains` are the **reference** frame (typically `D` for a VHH graft);
  the design's own chain is relabelled and renumbered by rfd3.
- `motifrmsd_design_chains` lists **every** design chain the RMSD covers. On a
  binder-scaffolding run that is legitimately four letters (`A,B,C,D`) and there is nothing
  in it to tell you which one the positions belong to. It was the right source before
  `--motif-chain` existed, only because the column refused to exist whenever it was
  ambiguous.
- `motifrmsd_fixed_positions_chain` is the design-side letter the spec actually addresses,
  and is NA exactly when the spec is.

Pair a reference letter, or the wrong one of four design letters, with design-side positions
and ProteinMPNN either redesigns the wrong chain in full (if that letter is not in the
design) or applies the motif positions to the wrong chain (if it is) — at exit 0, with no
status to catch it.

**The positions are 1-indexed FILE-ORDER positions in the exact structure file motifrmsd
read** — the `_path` cell of the input column at the time it ran — not residue numbers in
any other frame, and they are valid only if the downstream ProteinMPNN run's `-i` is **that
same column**. Point `-i` at a different structure (a merged complex, a relaxed or
renumbered file, a chainsel output) and every position silently shifts by however many
residues were inserted or removed before the motif, freeing part of the graft with no status
flag anywhere to catch it.

**No trailing slash, and that is load-bearing.** In the position mini-language
(`prosapia/utils/positions.py`) `/` breaks chains and the groups map **one-to-one onto
`--chains-to-design`, in order**. For a single designed chain a bare group is the whole
value; a trailing `/` declares a *second, empty* chain group, which for a one-chain
`--chains-to-design` is a group too many. Measured against the parser, not read from prose:

```
parse_positions("5:12,30:37")   -> [[5..12, 30..37]]        # 1 group  — correct
parse_positions("5:12,30:37/")  -> [[5..12, 30..37], []]    # 2 groups — wrong here
```

(The trailing slash is right only in the *two*-chain idiom
`--chains-to-design A,B --fixed-positions 1:609/`, where the empty second group leaves
chain B free.)

#### `{motifrmsd_design_fixed_positions}` does NOT work — measured

The obvious one-liner is a trap:

```bash
--fixed-positions "{motifrmsd_design_fixed_positions}"   # raises at submit time
```

`{…}` is prosapia's **expression** island and it resolves **integers only** — `utils/expr.py`
accepts an int literal, a bare column name, and `+ - * //`, then coerces. A column holding
`5:12,30:37` dies with:

```
ValueError: could not convert string to float: '5:12,30:37'
```

The value is a whole position *spec*, not a number *inside* one, so the substitution has to
happen **caller-side**, before `sapia run` sees it. The designs sharing one spec share one
run:

```python
ok = df[(df["motifrmsd_status"] == "OK")
        & df["motifrmsd_design_fixed_positions"].notna()]
# group on the chain letter TOO: the flag pair is (chain, spec), and the letter
# comes from the row, never from the reference frame.
for (chain, spec), group in ok.groupby(
    ["motifrmsd_fixed_positions_chain", "motifrmsd_design_fixed_positions"]
):
    # one proteinmpnn run per distinct pair: a -f filter selecting group.index,
    # plus --chains-to-design "<chain>" --fixed-positions "<spec>"
    ...
```

Expect **many** distinct specs — the span lengths are drawn per design, so designs from
different rfd3 batches will rarely share one. Group first and count before assuming this is
one extra run; if it is one run per design, that is the honest cost of per-design numbering,
not a mistake.

#### Why it can be NA on an otherwise-`OK` row

ProteinMPNN's `--position_list` is **positional**: each number is a 1-based index into that
chain's parsed sequence, *not* a PDB residue number (`make_fixed_positions_dict.py`; prosapia
says the same — "positions are 1-indexed within each chain"). Residue number and position
coincide **only** when the design chain is numbered `1..N` in file order with no gaps, which
is what rfd3 is understood to write.

The tool **checks** that rather than assuming it. The column is NA when:

- the design chain is not numbered `1..N` in file order (gaps, or not starting at 1) — the
  check is on the **whole chain as it stands in the file**, an absolute property, so it is
  unchanged by `--motif-chain` and is applied to whichever chain the restricted residues
  landed on;
- the motif — **after** any `--motif-chain` restriction — maps into **more than one** design
  chain, so no single chain group can say it. Without the flag this fires on every row of a
  binder-scaffolding run; see that section; or
- `--motif-chain` names a chain no reference key in that design's map mentions. That is a
  caller mistake in the flag, not a measurement failure: the `.err` names the chains that
  *are* in the map, and nothing else about the row changes.

In both cases **every other column, including `motifrmsd_status`, is untouched** — the row
can still be `OK` with a good `rmsd`. The reason is printed to that task's `.err`
(`design_fixed_positions not emitted (NA): …`). So an NA here next to a fine `rmsd` is a
statement about the *numbering*, not about the graft. Read `motifrmsd_design_chains`
alongside it. An NA across the *whole* table means rfd3's real output numbering is not
`1..N`, and the premise of this column needs revisiting before it is used.

A third NA case: **an incomplete match**. Under `--allow-partial` the other columns
(`rmsd`, `seg_ranges`, `seg_rmsds`, …) are measured on the matched subset as before, but
this column requires `n_matched == n_expected` and is NA whenever it is short — a spec over
the matched subset would leave the unmatched motif residues *free to be redesigned*, the
exact opposite of what the column is for. The `.err` carries
`design_fixed_positions not emitted (NA): only 15 of 16 motif residues matched …`. So a
`warn:` row never hands you a partial fixed-positions string to paste.

## A filter you can write

```python
# tools/filters/motif_faithful.py  — pass with -f
def apply_filter(df):
    return df[
        (df["motifrmsd_status"] == "OK")
        & (df["motifrmsd_worst_seg_rmsd"] < 1.0)
    ]
```

To also require that the reference measured against is byte-identical to what rfd3 itself
diffused — which `status == "OK"` does **not** imply when `ref_matches_spec` is NA — add
`& (df["motifrmsd_ref_matches_spec"] == 1)`. NA compares False, so this drops the rows where
the question could not be settled; count them before deciding that is what you want.

Two things to do alongside it, every time:

```python
# NaN compares False, so an error row is DROPPED without saying so.
# Report the three counts together or "sidecar unreadable" reads as "motif not kept".
n_ok   = (df["motifrmsd_status"] == "OK").sum()
n_err  = df["motifrmsd_status"].astype(str).str.startswith("error:").sum()
n_warn = df["motifrmsd_status"].astype(str).str.startswith("warn:").sum()
```

The `1.0` is a placeholder until a real distribution exists — see the threshold section.

To select the rows whose graft can actually be handed to ProteinMPNN — a faithful motif
**and** a usable fixed-positions string:

```python
def apply_filter(df):
    return df[
        (df["motifrmsd_status"] == "OK")
        & df["motifrmsd_design_fixed_positions"].notna()
        & (df["motifrmsd_fixed_positions_chain"] == "D")   # the binder chain
    ]
```

`notna()` rather than a comparison: the column is a string, and NA means the numbering check
did not pass (or the motif spans several chains), not that the graft is bad. The
`fixed_positions_chain` clause makes the hand-off explicit — every surviving row's spec
addresses the same design chain, so one `--chains-to-design` serves the whole group. Count
the dropped rows and read their `motifrmsd_design_chains` before accepting the loss; if the
count is the whole table and `design_chains` is multi-letter, you forgot `--motif-chain`.

## What it is blind to

- **Whether the motif is *bindable*.** It measures geometry against a reference, with no
  target present. Whether a scaffolded paratope still binds needs the target back —
  `epitope` for the footprint, `cms`/`pyrosetta` for the interface, and a prediction first.
- **Whether the scaffold clashes with the target** once the design is placed back in the
  target frame. No existing tool does that placement for a motif graft; `ringfit` does the
  structurally similar thing for an oligomeric assembly but with a two-protomer premise.
- **Whether the sequence folds to this backbone.** This is a backbone-vs-backbone measurement.
  Self-consistency is ProteinMPNN → Boltz/AF3 → `usalign`.
- **Non-CA atoms.** It is a CA RMSD. Side-chain rotamers and backbone O/N placement are
  invisible; a motif whose CA trace is perfect can still present the wrong side chains.
- **Non-motif residues.** The scaffold is never scored. A design with a perfect motif and a
  nonsense scaffold scores perfectly here — `pyrosetta` and the predictors judge the rest.
- **Whether ProteinMPNN will parse the chain the same way the tool did.**
  `design_fixed_positions` assumes ProteinMPNN's parsed chain is the structure's CA residues
  in file order. A residue ProteinMPNN drops for an incomplete backbone would shift every
  position after it. This tool indexes on CA only and cannot see that.
- **Whether `--motif-chain` names the right chain.** It is taken on trust: the tool checks
  only that the letter appears among that design's reference keys, never that it is the
  *grafted* chain rather than a held-fixed target chain. Pass `A` on a binder-scaffolding run
  and you get a perfectly valid spec for the target copy, which ProteinMPNN will happily use
  to redesign the binder's free positions while fixing a chain nobody meant to touch. Check
  `fixed_positions_chain` and the fragment ranges against the contig before the hand-off.
- **Whether `rfd3`'s own `diffused_index_map` is correct.** It is taken as ground truth. If
  rfd3 recorded the wrong correspondence, this tool faithfully measures the wrong thing; the
  only check available is `seq_match_frac`.

## Status

**Built 2026-10-08; `design_fixed_positions` added 2026-10-09; `--motif-chain` and
`fixed_positions_chain` added 2026-10-09 after real-data verification. Tests pass (255
checks, up from 177), which is not verification.** A `tool-reviewer`
code review (2026-10-08) found two gaps, both fixed the same day: the map is now refused
unless one-to-one (a duplicated design residue, or two key spellings of one reference
residue), and segments now break on design-side as well as reference-side discontiguity. A
bad `.src` cache record now costs one row instead of the whole submit. Still outstanding: a
`tool-reviewer` confirmation pass on those fixes, and a real-data run under a `verification`
label on `table0_seg2`/`table0_seg3`. **Verified on real data 2026-10-09** (5 designs of
`20261009_130629_5M2J_vhh_redmotif_lowtemp_scaleup`): rfd3 output chains ARE numbered `1..N`
in file order, with zero exceptions — the one unknown `design_fixed_positions` rested on. The
same run found the multi-chain problem this flag fixes, and the expected per-design D-chain
segments were derived independently (e.g. `scaleup0_diff_0` → `11:18,45:54,87:92`).

**Still outstanding for `--motif-chain` specifically:** a `tool-reviewer` pass on the new
code, and a real-data re-run against those same 5 designs checking
`motifrmsd_design_fixed_positions` and `motifrmsd_fixed_positions_chain` against the known
expected values above. Until that is recorded, treat the column as built-and-tested, not
verified. Do not gate a campaign on these columns until the campaign log records what
`seq_match_frac` and `ref_matches_spec` actually came back as.
