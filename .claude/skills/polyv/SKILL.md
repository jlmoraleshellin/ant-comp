---
name: polyv
description: >-
  Sequence-blind backbone triage with PyRosetta — erase the binder to poly-valine, repack only
  its side chains against a frozen target, and score the interface geometry (contact molecular
  surface, dSASA, shape complementarity) of a backbone that has no sequence yet. It RECORDS
  scores for you to read — a floor for discarding backbones that barely touch the target, never
  a ranker for picking the good ones, and this skill suggests no threshold. Load before scoring
  rfd3/RFdiffusion backbones ahead of sequence design, and before comparing any
  polyv_ number with anything else: its CMS is NOT comparable with the cms tool's and its dsasa
  is NOT comparable with pyrosetta's if_dSASA. Covers the chain-validation contract, why
  polyv_chains_found must be read before every science column, why dsasa is summed over BOTH
  partners (about 2x a one-sided BSA), and why a backbone that does not touch the target gives
  dsasa 0 with cms NA rather than 0.
---

# polyv

**Custom tool** (lives in `tools/polyv/`, not bundled with prosapia). Written entirely
with PyRosetta. **Modal-only and CPU-only**: the manifest builder forces
`gpus_per_task = 0`, so `-g 0` is never needed.

**`action: update`**: a geometric score is a property of a backbone that already exists,
so it annotates the table it reads, in place. **`-t` is required.**

## Premise, and when its numbers are meaningless

polyv answers exactly one question:

> **Does this backbone present a large, well-fitted surface to a fixed target, before any
> sequence exists?**

Per design it loads the structure, erases the binder chain(s) to **poly-valine**, repacks
**only the binder's** valine side chains against a **completely frozen** target, and measures
interface geometry. **Nothing is minimised or relaxed, anywhere** — speed is the point, and a
moved backbone is no longer the backbone that was triaged.

It is **not** an interface-quality score for a designed sequence, and it does **not** replace
`cms`.

| Do not compare `polyv_*` with | Why |
| --- | --- |
| `cms_target` / `cms_binder` (the `cms` tool) | Different implementation (Rosetta's `ContactMolecularSurface` filter vs cms-cuda) **and** a different pose (poly-valine vs the as-built sequence). |
| `pyrosetta`'s `if_dSASA` | Different pose, no FastRelax, and `pack_separated` is **off** here (rigid-body separation, not a repacked unbound state). |
| another polyv run with a different `--distance-weight` | The weight is inside the CMS sum. |
| the same design in `.pdb` and in `.cif` | PDB stores 3 decimals; measured ~0.1 % drift (see Gotchas). |

**Compare a polyv column only with the same polyv column, from the same run, on the same
target.** And see "What the columns are for" before comparing at all: polyv is a floor for
discarding backbones that barely touch, not a score for picking the good ones.

Its numbers are meaningless when:

- the input is not a binder/target complex **in one frame** — a binder docked somewhere else,
  or a structure where the "binder" chain is actually part of the target;
- the chain letters you passed do not mean what you think. polyv validates that they are
  *present*, never that they are *right*. **Read `polyv_chains_found` first.**
- you wanted chemistry. There is no energy, no dG, no hydrogen bond, no electrostatics and no
  buried-unsat count here, by design — and poly-valine has no polar side chain at all;
- **there are heteroatoms on a named chain that you did not mean to score.** Selection is purely
  by chain letter, so **everything Rosetta accepts on those chains is scored** — a ligand, a
  cofactor, a lipid or a modified residue sitting on the target chain contributes to `_dsasa`,
  `_cms_target` and `_n_res_target`. Name them in `--exclude-resnames` (and check
  `_n_res_dropped_resname`), or accept that they are part of the measurement. The flag
  takes **non-polymer** residues only — see the Gotchas for why "score it ignoring the
  histidines" is not something polyv will do.

## Determinism: one input file, one answer

Rosetta's packer is **stochastic**, and its RNG state carries across the designs inside one
task. Measured before the fix: the same byte-identical input file scored `cms_target` 402.4235
at one position of a task and 402.1243 at another, with `dsasa` 1597.13 vs 1599.69 and `sc`
0.4017 vs 0.3972 — i.e. `--designs-per-task` silently changed the science.

polyv now **reseeds Rosetta's RNG per design from a SHA-256 hash of the design name**, before
anything stochastic runs. Consequences you can rely on, each covered by a test:

- a design's numbers are **independent of its position in the batch, of `--designs-per-task`
  and of row order**;
- a rerun with `--force` reproduces the previous numbers exactly;
- the seed is collected as **`polyv_seed`**, so any number can be traced to the state that
  produced it.

**Does the seed *value* change the answer? Measured: no, on every input tested.** Three real
BindCraft2 backbones (binder lengths 86, 115, 122) scored under three different names each — nine
runs, pairwise-distinct seeds — came back **identical to every emitted decimal** on `cms_target`,
`cms_binder`, `dsasa` and `sc`; two synthetic fixtures (binder lengths 20 and 40) × three seeds
did the same. Valine has a single chi and the target is frozen, so the annealer has very little
to get wrong. **`polyv_seed` is therefore kept for traceability, not because it is expected to
matter** — and two rows with different names on one file have, so far, always agreed.

Scope that claim honestly: *seed-independent on every input tested* (5 inputs, one chain
configuration, binder lengths 20–122), **not** "the packer is deterministic in general". The
reseeding stays regardless, because it is what makes a rerun reproducible and removes
`--designs-per-task` from the set of things that can change a number — **confirmed on real data
on Modal**: the same 10-row table scored as one task of 100 and as ten tasks of 1 gave 10/10
names identical on all four float columns, seeds matching per name (see Verified).

### Dated note — 2026-10-07: every earlier number is non-comparable

The per-design reseeding **changes the float columns of every row collected before it**. Rows
from the pre-fix verification run cannot be compared with rows collected afterwards, and they
cannot be reproduced. Re-run those designs with `--force` before using their numbers.

Column renames landed at the same time. **`polyv_n_input_nongly` is gone**, replaced by
`polyv_n_input_restypes` (a count of distinct residue names, which also catches a poly-ALA or
poly-VAL backbone); and **`polyv_n_res_dropped` is gone**, split into
`polyv_n_res_dropped_chain` and `polyv_n_res_dropped_resname`. A table collected earlier keeps
the old columns and will never have the new ones filled; re-collect with `--force` if you need
them.

## Trust columns — read these before the science

| Column | Expected | What a surprise means |
| --- | --- | --- |
| `polyv_chains_found` | exactly the chains you expect, e.g. `A,B` | The chain IDs **as loaded**. If a `.cif` was remapped on the way in, or a predictor split a protomer across two letters, you see it here. A chain that appears here but is in neither `--binder-chains` nor `--target-chains` was **silently dropped** from the calculation. A blank chain ID is shown as `_` and can never be selected. |
| `polyv_n_res_binder` | the binder length you designed | Residues **scored on** the binder side, i.e. **after** `--exclude-resnames`. A wrong-but-present chain scores without complaint; this is the check. |
| `polyv_n_res_target` | the full target, **all protomers**, plus any heteroatom riding on those chains | Residues **scored on** the target side, i.e. **after** `--exclude-resnames`. Half of this means you forgot a target chain; one or two more than expected usually means a ligand or ion came in on the chain. |
| `polyv_n_mutated` + `polyv_n_kept_gly_pro` | `== polyv_n_res_binder` | Enforced on **both** paths — with `--keep-gly-pro` off, `n_kept_gly_pro` is 0 by construction, so a residue that failed to mutate cannot hide in the gap. An `OK` row that violates it is impossible; the error names the residue types still present. |
| `polyv_n_input_restypes` | `1` for a uniform backbone; `>5` for an already-sequenced pose | Distinct residue **names** on the binder in the input. `1` means poly-GLY *or* poly-ALA *or* poly-VAL — any uniform backbone, whatever letter the generator used. A real sequence reads 8–20 and means polyv **erased a sequence you may have thought you were measuring**. (Measured on real sequenced poses with the older non-GLY count: 82/86, 109/115, 120/122 — i.e. essentially fully sequenced.) |
| `polyv_n_res_dropped_chain` | `0` unless you expected a chain to go | Residues deleted because their chain was in neither selection. **The complementary half of the `n_res_*` check, not a duplicate of it:** `n_res_binder`/`n_res_target` say what was *scored*, this says what was *removed*, and only together do they account for everything in `chains_found`. Non-zero means a whole chain went. |
| `polyv_n_res_dropped_resname` | `0` with no flag; exactly what you meant to remove with one | Residues deleted because `--exclude-resnames` named them. A `0` here against a non-empty flag is impossible — a resname matching nothing is an error row, not a quiet no-op. |
| `polyv_seed` | any 31-bit int; `sha256(design name)` and nothing else | **Derived, not observed.** It is a function of the row's name alone, so it can never surprise you and can never flag a problem on its own — **do not write a filter against it.** Its value is traceability: it lets you say which RNG state produced a number, and lets you check that two rows were scored from the same state. |
| `polyv_seconds` | seconds | Wall time for this design. |

**`polyv_status` is the only proof of success.** Failures exit 0 (errors are recorded as data).

## What the columns are for: a floor, not a ranking

**polyv records scores. It does not choose designs.** Read this before you reach for a cutoff.

**The informative end of every polyv column is the bottom.** A near-zero interface is strong
evidence that a backbone is not worth sequencing: it barely touches the target, and no sequence
design can fix a backbone that is not there. A *large* polyV interface is weak evidence of
anything — it says the geometry could in principle present a surface, not that a designed
sequence will be able to exploit it. Valine is not a design; it is a probe.

So the instrument is **a floor for discarding the backbones that barely touch**, not a ranker
for picking the good ones. Use it to remove the obvious failures and carry everything else
forward; the decision about what is actually good belongs downstream, to sequence design and
prediction, and to the person running the campaign.

Two consequences:

- **This skill does not suggest a threshold, and will not.** There is no published bar for
  poly-valine CMS, polyv's values are not the `cms` tool's, and the right floor depends on the
  target, the epitope and how much sequence design you can afford. Any number you use is yours;
  record it in the campaign log so later generations are treated the same way.
- **The one unambiguous discard is `polyv_nres_int == 0`** (equivalently `polyv_dsasa == 0`,
  with `polyv_cms_target` NA): the binder is not touching the target at all. That is a property
  of the row, not a cutoff anyone chose.

### The polyV scale, measured on 8 real verification designs

Context for reading a number, not a bar to apply:

| Column | Observed range | How to read it |
| --- | --- | --- |
| `polyv_sc` | **0.1811 – 0.4017** | **Literature SC for natural and designed interfaces is 0.5–0.75. polyv never reaches it and is not supposed to.** Valine cannot fill a pocket the way a designed sequence does, so `polyv_sc` lives on its own scale — **never compare it with a published SC number, or with `cms_sc` / `pyrosetta`'s `if_sc`.** |
| `polyv_cms_target` | **75.6 – 402.4 Å²** | A >5× spread across 8 designs, so the column does separate backbones rather than flattening them — which is what makes the low end readable at all. |

Those are the ranges from one target and 8 designs: enough to recognise what "near the bottom"
looks like on this scale, not enough to be a threshold anywhere else.

### Two observations about what the columns measure

Neither is advice about how to combine them.

**`sc` ran roughly opposite to interface size.** On those 8 rows the **second-largest**
`cms_target` (343.7) had the **lowest** `sc` in the set (0.1811), while the two smallest
`cms_target` (129.2 and 75.6) had among the highest `sc` (0.3653 and 0.3229). Large-and-sloppy
versus small-and-tight. They are measuring different things and they do not agree about which
backbone is better — which is one more reason not to treat either as a selector.

**`dsasa` tracked `cms_target` closely** — one inversion in 8 rows. It comes out of different
machinery (InterfaceAnalyzerMover's SASA, not the CMS filter's surface dots), so a row where the
two disagree is a signal that one of them is wrong. It is also the only one of the two that is
an honest `0` rather than NA when nothing touches, which is exactly the case the floor is for.

## Invocation

```bash
sapia run polyv <run_dir> -t table0 -i rfdiffusion3_path \
    --binder-chains B --target-chains A,C -l bb
# poll .exit as usual (Modal), then:
sapia collect polyv <run_dir> -t table0 -l bb
```

| Flag | Default | Notes |
| --- | --- | --- |
| `-i/--input-column` | **effectively required** | Sentinel default `"not applicable"` (the cms/usalign/chainsel convention). The builder **raises** unless the table has the column. Reads `.pdb` or `.cif`. |
| `--binder-chains` | **required** | Comma-joined. These chains are mutated to VAL and repacked. |
| `--target-chains` | **required** | Comma-joined; the target may be multimeric (`A,C`). Frozen entirely. Must be disjoint from the binder — checked at submit time. |
| `--keep-gly-pro` | **off** | Off means GLY and PRO on the binder are mutated to VAL too, which is what makes the probe uniform. It does **not** weaken the trust contract: the check becomes `n_mutated + n_kept_gly_pro == n_res_binder`, recounted from the pose and exactly as strong, so a residue that failed to mutate still cannot hide in the gap. |
| `--extra-rotamers` | `ex1` | `none` \| `ex1` \| `ex1ex2`. Valine has one chi, so `ex1ex2` buys little for a lot of time. |
| `--distance-weight` | `0.5` | `distance_weight` of the CMS filter (the Baker-lab convention; **Rosetta's own XML default is 1.0**). Changing it makes `cms_*` non-comparable with every other run. |
| `--exclude-resnames` | `""` (**nothing excluded**) | **Non-polymer** residue names dropped before scoring, comma-joined (`HEM,ZN`) — a ligand, ion, lipid or cofactor. Without it polyv scores every residue Rosetta accepted on the named chains, heteroatoms included. **A polymer residue is refused** and **a name matching nothing is an error** (see Gotchas). Removals show up in `_n_res_dropped_resname` and in `_n_res_binder` / `_n_res_target`. |
| `--designs-per-task` | `100` | See batching below. |
| `-l/--dir-label` | `""` | **In practice required.** See the traps. |
| `-g/--gpus-per-task` | forced to `0` | CPU-only; do not pass it. |

Resources (`RESOURCES` in `modal_image.py`, and the `#SBATCH` lines): 2 CPU, 8 G, 4 h per task.

### Batching: one task holds many designs

`pyrosetta.init()` costs ~3–4 s per container, so the builder packs `--designs-per-task`
designs into a sub-manifest (`<out_dir>/polyv_tasks/task_<t>.tsv`) and submits **one task per
chunk**. **`n_tasks` in `polyv_modal.json` is the number of chunks, not the number of
designs** (Modal only). With the default of 100, 250 designs make 3 tasks. Since the per-design
reseeding it is a **scheduling knob only**: changing it cannot change a single number.

Per-design cost measured only on 46–66-residue toy helices: **0.5–2.9 s**. A real
100-residue binder on a few-hundred-residue target will be slower and is **unmeasured** —
treat the first real run as a timing shakedown before setting `--designs-per-task` for a big
table.

### Image (Modal only)

`tools/polyv/modal_image.py` is **byte-identical, layer for layer**, to the bundled
`pyrosetta` tool's image, so Modal reuses the ~1.5 GB PyRosetta wheel this workspace already
built instead of downloading it again. **Keep the two in step** — changing a character forks
the build cache and costs another multi-minute build.

**Not runnable on vib.** There is no `SAPIA_ACTIVATE_POLYV` entry; the task script calls
`sapia_activate SAPIA_ACTIVATE_POLYV`, so on the cluster every task dies at activation. The
signature is **`sapia: set SAPIA_ACTIVATE_POLYV in your .env to a tool activation script` on
stderr** (the prelude writes it before exiting) — *not* the empty-`.out`-and-`.err` signature,
which CLAUDE.md reserves for a death before the first statement. Same state as the bundled
`pyrosetta` tool.

## Gotchas, each with its signature

**An absent chain is an error; a wrong chain is not.** A chain named but missing gives
`polyv_status = "error: chain X not in <file> (have: A,B,C)"`. A chain that is present but
*wrong* — the target named as the binder, one protomer forgotten — scores without complaint.
Signature: `polyv_n_res_binder` is the target's length, or `polyv_n_res_target` is half what
you expect.

**A chain in neither list is silently dropped.** InterfaceAnalyzerMover splits the pose into
exactly two partners, so anything else is deleted before scoring. Signature:
`polyv_chains_found` lists a letter you did not pass. This is correct behaviour for a ligand
or a spectator chain and a **bug in your command** for a forgotten target protomer.

**`cms_target` / `cms_binder` NA is not 0.** Rosetta's `ContactMolecularSurface` filter
**raises a message-less `RuntimeError`** when the two selections do not touch, instead of
returning 0. polyv records NA and keeps `polyv_status = "OK"`. Signature: `polyv_cms_target`
empty while `polyv_dsasa == 0` and `polyv_nres_int == 0` — the binder is not touching the
target. That is an honest, very bad result, not a failure.

**`sc` NA is not 0.** Rosetta reports `sc_value = 0.0` when the SC calculation fails (it can,
on a tiny interface). polyv converts any non-finite or ≤ 0 value to NA. Signature: `polyv_sc`
empty on a row whose `polyv_dsasa` is positive — the interface was too small or degenerate for
Lawrence–Colman SC.

**`dsasa_per_res` is NA when `nres_int` is 0**, because it is undefined, not because it is
zero.

**`.pdb` and `.cif` of the same pose differ by ~0.1 %, and the cause is coordinate precision.**
Measured with the seed pinned: `cms_target` 214.4614 vs 214.733 (0.13 %), `dsasa` 862.3364 vs
861.3117, `sc` identical at 0.4246, every integer column identical. The difference also survives
when the seed is *varied*, so it is not an RNG artefact — PDB stores coordinates to 3 decimals
and mmCIF keeps more, and that is the whole of it. Do not mix formats inside one comparison.

**Several models in one file are refused.** A chain that reappears after another chain gives
`error: chain X is not contiguous (several models, or a repeated chain ID, in one file)`.
Without this every count would silently double.

**`dsasa` is summed over BOTH partners.** InterfaceAnalyzerMover's delta SASA counts the
surface buried on the binder *and* on the target, so it is roughly **2×** the single-side buried
surface area (BSA) quoted in the binder literature. A 1600 Å² `polyv_dsasa` is an ~800 Å²
one-sided interface. Nothing derived from it is wrong — `dsasa_per_res` and every comparison use
the same convention — but a threshold carried over from a paper will be out by a factor of two.

**`--exclude-resnames` is empty by default**, so a ligand, ion or cofactor on a named chain is
part of the measurement until you say otherwise. Signature: `polyv_n_res_target` one or two
above the target length you expect, with `polyv_n_res_dropped_resname` at 0.

**`--exclude-resnames` will not exclude an amino acid, by design.** "What does this interface
look like ignoring the histidines" is refused, at submit time and again in the worker:
`error: --exclude-resnames names polymer residue HIS; the flag is for non-polymer heteroatoms`.
Deleting an interior polymer residue leaves a hole that exposes backbone and neighbour atoms it
previously occluded, so the SASA, CMS and SC of every **remaining** residue change for reasons
that have nothing to do with the question — believable numbers on a mutilated structure. The
submit-time list covers the canonical amino acids and nucleotides; the worker's authoritative
test is Rosetta's own `is_polymer()`, so an exotic polymer residue is still caught, just later.

**A resname that matches nothing is an error, not a no-op.** Mirrors the absent-chain contract:
`error: --exclude-resnames UNL matched no residue on chains A,B; a requested resname that is
absent is an error, not a smaller selection`. Without it, typing `LIG` for `UNL` or `HEME` for
`HEM` — or naming something Rosetta already dropped under `-ignore_unrecognized_res` — would
leave the thing you meant to exclude in the measurement with `status OK`. Matches are counted
only on the **scored** chains: a hit on a chain that was dropped anyway excluded nothing.

**`_` is refused as a chain ID**, at submit time and in the worker. It is how polyv *renders* a
blank chain ID in `chains_found`, and it is also the separator in InterfaceAnalyzerMover's
partner string. Signature: `'_' is polyv's rendering of a BLANK chain ID …`. Give the structure
real chain IDs with `chainsel` first.

**Disulfides are not detected** (`-detect_disulf false`). A binder about to become poly-valine
cannot keep them, and nothing polyv reports is energetic, so this changes no number — but it
means the packed pose is not a pose you should score for energy.

**`-l/--dir-label` is effectively required.** The leaf `polyv` names both the output dir and
every column, so triaging `rfdiffusion3_path` and then a predicted complex on the same table
with no label overwrites the first run's columns. Use `-l bb` and `-l pred` to get
`polyv_bb_*` and `polyv_pred_*`. **Pass the same `-l` to `collect`.**

**Failures exit 0.** Errors are recorded as data, so `.exit` files that are all `0` (Modal) do
not mean every design succeeded. Only `polyv_status` is reliable. A worker crash writes
`error: worker crashed` for every design in that chunk.

**Rows can drop out of the manifest silently.** A design whose input file does not exist is
printed as `MISSING … (skipping)` and left out; its row collects as `missing`. Signature:
`Submitting N designs` is below the table's row count.

**Chain letters, not residue numbers.** rfd3 renumbers every output chain from 1, which does
not affect polyv (it never names a residue number) — but generators do change chain *letters*
between steps. `polyv_chains_found` is the only thing that catches that.

**On tiny poses `nres_int` counts everything.** On two bare 20–26-residue helices
InterfaceAnalyzerMover called all 46 residues interface residues, because a bare helix changes
SASA everywhere. On a real target it is a subset. Do not calibrate `dsasa_per_res` on toys.

## What it collects (`polyv_` prefix, or `polyv_<label>_` with `-l`)

Trust columns first; read them before any measurement.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `_status` | **trust** — `OK`, `error: <reason>`, `missing` | The only proof of success. Filter on it first, always. |
| `_chains_found` | **trust** — every chain ID in the pose as loaded, comma-joined, in pose order | Must contain exactly the letters you passed. Extras were dropped from the calculation. Blank IDs show as `_`. |
| `_n_res_binder` | **trust** — binder residues **scored** (after `--exclude-resnames`) | Compare with the binder length you designed. |
| `_n_res_target` | **trust** — target residues **scored** (after `--exclude-resnames`) | Compare with your target, **all protomers** — plus any heteroatom riding on those chains. |
| `_n_mutated` | **trust** — binder residues that are VAL in the packed pose | With `_n_kept_gly_pro` it must sum to `_n_res_binder`; enforced on both paths, so an `OK` row that violates it is impossible. |
| `_n_kept_gly_pro` | **trust** — binder GLY/PRO deliberately left alone | `0` unless `--keep-gly-pro`. This is what makes "the gap is GLY/PRO" checkable instead of asserted. |
| `_n_res_dropped_chain` | **trust** — residues deleted because their chain was in neither selection | `0` is the normal case; filter on it. With `_n_res_binder`/`_n_res_target` it accounts for everything in `_chains_found`. |
| `_n_res_dropped_resname` | **trust** — residues deleted by `--exclude-resnames` | `0` with no flag. A `0` against a non-empty flag cannot occur: that is an error row. |
| `_n_input_restypes` | **trust** — distinct residue **names** on the binder in the input | `1` = a uniform backbone, whatever letter (poly-GLY / poly-ALA / poly-VAL all read 1). `>5` = a real sequence that polyv erased. |
| `_seed` | **trust** — the RNG seed, `sha256(design name)` | **Derived, not observed — never filter on it.** For traceability only: it says which RNG state produced the row. Measured seed-independent on every input tested, so a differing seed has not changed a result. |
| `_seconds` | **trust** — wall time for this design | Sanity/cost. |
| `_cms_target` | **Contact molecular surface of the TARGET**, weighted toward the polyV binder, Å² | The headline number: large *and* tight. **NA** when the filter raised (no contact). |
| `_cms_binder` | the same with the sides swapped — the binder's surface, Å² | Usually within a few % of `_cms_target`. A large gap means one side is much more enclosed than the other. |
| `_dsasa` | interface dSASA **summed over both partners**, Å² | Interface **size**, unweighted by tightness. **About 2× the one-sided BSA** quoted in the literature — halve it before comparing with a published number. `0` is a real measurement: nothing is touching. |
| `_dsasa_per_res` | `_dsasa / _nres_int` | Interface density. **NA** when `_nres_int` is 0. |
| `_nres_int` | interface residues, both sides | `0` ⇒ the binder is off the target. |
| `_nres_int_binder` | interface residues on the binder | Cross-checked against InterfaceAnalyzer's own side-1 count; a disagreement is an error row. |
| `_nres_int_target` | interface residues on the target | How much of the target is engaged. |
| `_sc` | Lawrence–Colman shape complementarity, 0–1 | Fit quality. **NA** (not 0) when Rosetta's SC calculation failed. |
| `_path` | the **packed poly-valine pose** (PDB), run-relative | Open this to see what was actually scored. It is a polyV binder on an untouched target, with no relax — do not score it for energy. |

**NA (an empty cell) means "could not be computed", never 0.** `0` is a measurement. On disk:
`<run_dir>/<table>/polyv[_<label>]/<name>.tsv` and `<name>.pdb`, plus the sub-manifests under
`polyv_tasks/`.

**No energy column is collected, ever** — no dG, no hbonds, no packstat, no score terms. That
is `pyrosetta`'s job.

## Writing an `-f` filter against these columns

**polyv's job is to record scores, not to choose designs.** The intended use is to read the
columns and decide; this section is here for the mechanics of `-f`, and for the two filters that
are about *correctness* rather than selection.

### 1. The trust pre-check — run this on every new target

This is the filter that is genuinely about the table, not about taste. It selects the rows whose
chain handling went wrong, so **anything it returns is a bug in your command**, not a bad design:

```python
# filters/polyv_trust.py
def apply_filter(df):
    return df[
        (df["polyv_bb_status"] != "OK")
        | (df["polyv_bb_chains_found"] != "A,B")        # the chains you expected
        | (df["polyv_bb_n_res_target"] != 274)          # the target length you expected
        | (df["polyv_bb_n_res_dropped_chain"] != 0)     # no chain should have been deleted
    ]
```

`"A,B"` and `274` are properties of **your** target — record them in the campaign log. If this
selects nothing, the numbers are about the structures you meant.

### 2. The floor — discarding backbones that barely touch

The only discard this skill will state outright is the one that needs no judgement: a backbone
that does not touch the target at all.

```python
# filters/polyv_touching.py
def apply_filter(df):
    return df[(df["polyv_bb_status"] == "OK") & (df["polyv_bb_nres_int"] > 0)]
```

If you want a floor above zero, the shape is this — and **`MIN_CMS` is a placeholder. The number
is the caller's to choose; this skill does not suggest one, and there is no published bar for
poly-valine CMS**:

```python
# filters/polyv_floor.py
MIN_CMS = ...   # <- YOU pick this, from your own batch, and log it. Not a default.

def apply_filter(df):
    ok = df[df["polyv_bb_status"] == "OK"]
    return ok[ok["polyv_bb_cms_target"] >= MIN_CMS]
```

Two things not to do with it:

- **Do not invert it into a selector.** This is a floor for removing the obvious failures, not a
  ranking for picking winners — a high `cms_target` is weak evidence (see "What the columns are
  for"). Everything above the floor should go forward and be judged by sequence design and
  prediction.
- **Do not borrow a bar from the literature.** An earlier version of this skill suggested
  `sc >= 0.55`; on the measured polyV scale `sc` tops out around 0.40, so that filter selects
  **zero rows** on every table polyv will ever write.

## What it is blind to

- **Chemistry.** No hydrogen bonds, no charge complementarity, no buried unsatisfied polars, no
  binding energy. Poly-valine cannot even form a polar contact. → `pyrosetta` (`if_dG`,
  `buried_unsat`), on a real sequence.
- **Whether a sequence exists that realises this surface.** A beautiful polyV interface can be
  undesignable. → `proteinmpnn` / `atomium`, then `boltz`.
- **Whether the backbone folds at all.** polyv scores the pose it is handed. → self-consistency
  with `usalign` against the designed backbone.
- **Where on the target it binds.** polyv reports *how much* contact, never *which* epitope.
  → `epitope` (hotspot recall, contact footprint) or `ifacegeom`.
- **Whether the binder clashes with the rest of an assembly or a membrane.** → `ringfit`.
- **Whether the complex is real.** It is a geometry calculation on whatever pose you gave it;
  prediction confidence is a separate gate entirely. → `boltz`.

## Verified

### On real data (Modal, verification run)

- **10/10 `OK`** on BindCraft2 trajectory **`.cif`** files, target `A,B`, binder `C`.
  ~**2.5–3 s per design**.
- Chain selection correct on real mmCIF: **`.cif` chain IDs came through as auth chains**.
- `n_mutated == n_res_binder` on every row; the packed output was confirmed **all-VAL on the
  binder chain**.
- Interface counts on one design: **16 of 115** binder residues and **24 of 296** target
  residues — so `get_side1_nres()` really does mean *interface* residues, strictly
  `0 < 16 < 115`, not the whole side.
- The deliberately bad-chain run **errored on every row with no numbers**:
  `error: chain Z not in <file> (have: A,B,C)`.
- `n_input_nongly` (the column that `n_input_restypes` replaced) read **82/86, 109/115,
  120/122** on sequenced poses — i.e. those poses arrived essentially fully sequenced, and
  polyv erased that sequence. `n_input_restypes` reports the same situation as a count well
  above 1, and unlike the old column it also catches a poly-ALA or poly-VAL backbone.
- **The polyV scale**, on 8 designs: `sc` 0.1811–0.4017, `cms_target` 75.6–402.4 Å². See
  "Reading the numbers".
- **Seed-independence**, 3 backbones (binder lengths 86, 115, 122) × 3 names each: all nine
  results identical to every emitted decimal on `cms_target`, `cms_binder`, `dsasa` and `sc`,
  with pairwise-distinct seeds confirmed.

### Determinism (the defect found in review, now fixed)

- **Before:** two rows with a byte-identical input file agreed with each other but disagreed
  with a third row holding the same file earlier in the same task — `cms_target` 402.4235 vs
  402.1243, `cms_binder` 377.4512 vs 379.0115, `dsasa` 1597.1302 vs 1599.69, `sc` 0.4017 vs
  0.3972. Integer columns were identical throughout. Rosetta's RNG state carries across designs
  inside a task, so `--designs-per-task` silently changed the science.
- **After:** the RNG is reseeded per design from `sha256(name)`, unconditionally and at the top
  of the per-design loop. Four tests pin it — the same design at position 0 and at position 5 of
  one task, reversed row order, one-per-task vs four-per-task, and a **failing design
  interleaved between two good ones** (a design that errors mid-pack has already consumed RNG;
  its neighbours must score identically to running them alone) — each requiring `cms_target`,
  `cms_binder`, `dsasa` and `sc` identical **to every decimal**. All pass. `polyv_seed` records
  the seed.
- **Confirmed on real data, on Modal:** the same 10-row table scored with
  `--designs-per-task 100` (one task) and with `--designs-per-task 1` (ten tasks) gave
  **10 of 10 names identical on all four float columns to every emitted decimal**, read from the
  raw per-design TSVs, with seeds matching per name across both runs. The local tests show the
  symptom is gone in-process; this shows it is gone across container boundaries and batch
  layouts, which is the claim that actually matters. **`--designs-per-task` cannot change a
  number.**
- **Separately measured, and a different question:** the seed *value* does not appear to change
  the answer at all — see "Determinism" above. The reseeding is kept regardless, because it is
  what makes a rerun reproducible. **No claim is made anywhere about the mechanism** of the
  original variation; the fix is justified by the before/after measurements, not by a theory of
  what carried state across designs.

### Local test suite

**`tests/test_polyv.py`: 23 passed** under the repo venv (no PyRosetta needed).
**`tests/test_polyv_worker.py`: 38 passed under PyRosetta 2026.30**
(`PyRosetta4.Release.python314.ubuntu 2026.30+release.bc091c65`).

⚠️ **Check the count, not the exit code.** `test_polyv_worker.py` opens with
`pytest.importorskip("pyrosetta")`, so on an interpreter without PyRosetta it reports **0 tests
and exits 0**. "All pass" is a real claim only when it is 38; otherwise it is vacuous. Run it as
`<a python with pyrosetta> -m pytest tests/test_polyv_worker.py -q` and read the number. Established on synthetic ideal helices: every binder residue really is VAL
in the packed pose and the target sequence is unchanged; the target's coordinates come back
**bit-identical** (no relax, no repack); `n_mutated + n_kept_gly_pro == n_res_binder` on both
paths; an absent chain, a missing file, a single-chain pose, overlapping chain lists and the `_`
sentinel all error rather than scoring; a non-contacting pose gives `dsasa 0`, `nres_int 0`,
`cms` NA, `sc` NA, status `OK`; a chain in neither list is dropped, counted in `n_res_dropped_chain`,
and the remaining pair's numbers match the same pair scored alone to every decimal;
a real non-polymer ligand (a `HETATM` zinc on the target chain) is scored as part of the target
unless named, and excluding it restores the ligand-free numbers **to every decimal**; excluding
a polymer residue and naming a resname that matches nothing are both error rows;
`n_input_restypes` reads 1 for poly-GLY and poly-ALA binders alike; nothing is written into the
input file's directory.

**The CMS direction is pinned by a test**, not by memory: `cms_target` is reproduced to 1e-2 by
a RosettaScripts `<ContactMolecularSurface target_selector=… binder_selector=…/>`, and the
swapped call is not. So `selector1 == target_selector == the side whose surface is measured`.

### Still unverified

- **Heteroatom handling on a real target.** `--exclude-resnames` is exercised against a real
  non-polymer residue (a `HETATM` ZN on the target chain) only in the local suite; no Modal
  verification run has yet included a ligand, cofactor or lipid on a scored chain.
- **`--keep-gly-pro` on real data** — exercised only on synthetic fixtures.
- **Wall time at scale.** 2.5–3 s per design measured on BindCraft2 trajectories with a
  115-residue binder and a 296-residue target; `--designs-per-task 100` therefore implies ~5 min
  per task, comfortably inside the 4 h timeout, but a much larger target is unmeasured.
- **Modal image layer reuse** with the bundled `pyrosetta` tool is assumed from the identical
  build commands, not observed.
