---
name: linkpath
description: How to run the custom linkpath tool — measuring the shortest route between two chain termini that stays OUT of the protein, i.e. the length a flexible linker would actually have to span. Covers why a straight C-to-N distance is only a lower bound, the obstacle-chain flag that is the whole scientific content of the measurement, the 2.0 Å chain-sized probe vs the prototype's 1.4 Å water probe, the straight_ca_dist bridge to dimerfit's retired link_dist, the "no path found is a RESULT not a failure" contract, the C2 path_asymmetry trust metric, and every column it collects. Load before composing a linkpath run or interpreting its columns.
---

# linkpath

**Custom tool** (lives in `tools/linkpath/`, not bundled with prosapia).
**`action: update`** — annotates the table it reads, in place. `-t` is required.
**CPU-only**; the manifest builder forces `gpus_per_task = 0`, so you do **not** pass `-g 0`.
**No default input column** — `-i` is effectively required (see Traps).
**Fast**: **0.161 s per design measured** (216-residue C2 dimer, both directions,
default grid), so one batched task covers a whole table.

> **Premise.** *Measure the shortest route between two chain termini that stays out
> of the protein — the length a flexible linker would actually have to span.*

A straight-line C-term-to-N-term distance is only a **LOWER BOUND**. If the direct
vector passes through the protein body, the linker has to go around. This tool
computes the shortest path through solvent-accessible space:

1. build a 3D occupancy grid over the chosen **obstacle chains** — a voxel is BLOCKED
   when its centre lies within (vdW + `--probe`) of an atom;
2. pick the two anchor atoms (backbone **C** of the donating terminus, backbone **N**
   of the accepting one, falling back to CA);
3. free a bubble around each anchor so the chain can leave its own terminus, then
   **A\*** over the free voxels (26-connected, true Euclidean edge costs, admissible
   Euclidean heuristic);
4. report the route length, the detour over the straight line, and the residue count
   that implies;
5. optionally route the **reverse** pair too — under C2 symmetry a free trust metric.

## Why it exists

The pH-switch strategy needs the two protomers of a C2 dimer **linked**. Deciding
whether a 20-residue linker reaches means knowing how far a chain actually has to
travel — and between two protomers that chain almost never goes straight: it has to
get around one or both bodies.

Nothing in the workspace could answer that. The collectors were read, not the prose:

- **`dimerfit`** (columns: `align_rmsd`, `n_align_atoms`, `seq_match_frac`,
  `resnum_offset`, `resnum_match`, `dock_chains`, `ref_binder_chain`,
  `ref_target_chains`, `complex_target_chains`, `n_protomers`, `n_res_a/_b`,
  `ref_binder_len`, `n_atoms_b`, `dimer_iface_res`, `dimer_iface_com`, `epitope_com`,
  `epitope_com_dist`, `n_epitope_res`, `overlap_res`, `frac_overlap`, `n_clash`,
  `clash_frac`, `min_dist_b_target`, `occluded_res`, `occluded_frac`) had exactly one
  linker number, the **straight-line CA–CA `link_dist`**, now retired. A straight
  line is the lower bound this tool replaces, not a competitor to it.
- **`ifacegeom`** (`binder_res`, `n_binder_res`, `binder_res_seq`, `target_res`,
  `target_his`, `binder_com`, `binder_iface_com`, `iface_com`, `axis_len`,
  `nterm_res`, `cterm_res`, `cterm_proj`, `nterm_proj`, `cterm_iface_dist`,
  `cterm_iface_min_dist`, `binder_len`) says **where** a terminus is relative to the
  interface — a projection and a distance, both straight-line, both within one body.
- **`ringfit`** (`align_rmsd`, `seq_match_frac`, `resnum_offset`, `bsa_t1/t2/total`,
  `bridge_ratio`, `hotspot_recall`, `n_clash`, `min_dist_ring`, `lipid_clash`) and
  **`cms`** (`target`, `binder`, `sc`, `sc_area`, `sc_median_dist`, `n_atoms_*`) and
  **`pyrosetta`** (`if_dG`, `if_dSASA`, `if_hbonds`, `if_delta_unsat`, `packstat`, …)
  and **`usalign`** (`TM1`, `TM2`, `RMSD`, `ID*`, `L*`) have no path search of any
  kind. None of them has a notion of routing *around* anything.

There is a verified single-file prototype at the repo root, `linker_path.py`. It is
a **script**: it prints to a terminal. This is a **tool** because the answer has to
be a filterable column that rides lineage — "does a 20 aa linker reach?" is a gate,
and a gate has to live in the table.

## Scope — where these numbers are meaningless

- **The obstacle set is whatever chains you point it at.** That choice *is* the
  measurement. Route a dimer through only its own two protomers while a target sits
  in the same file and you get a dishonestly short answer. **Measured**: the same C2
  barnase dimer routed through both protomers gives 52.80 Å and through one protomer
  only 49.60 Å. The resolved list is written back as `obstacle_chains` so a later
  reader can audit it without re-deriving anything.
- **The path is computed on a RIGID structure.** Real linkers are flexible and real
  termini move. This is a **geometric floor**, not a thermodynamic statement: it says
  nothing about whether a linker of that length is stable, soluble, entropically
  affordable, or whether it will fold back onto the interface.
- **The grid path is itself a LOWER BOUND.** A real chain has volume and cannot hug
  the surface exactly; the route is a 26-connected voxel path, not a smooth curve.
  **Add margin.** `n_res_relaxed` already builds some in (see the rises).
- **Existing unresolved or flexible terminal residues already contribute length** a
  designer may not need to add. The anchors are the last/first **resolved** amino
  acid; anything disordered beyond them is free length this tool cannot see. Check
  `from_res`/`to_res` against the construct's real termini before adding residues.
- **No energy, no sequence, no flexibility.** `n_res_min`/`n_res_relaxed` are a
  contour-length division, nothing more. It is not a linker *designer* and not a loop
  closure test.
- **Waters are never obstacles**; hydrogens are not either unless `--include-h`.
- **It does not supersede the straight-line distance, it bounds it from above.**
  `straight_ca_dist` is collected precisely so old and new rows stay comparable.

## Invocation

For this campaign — linking the two protomers of a `dimerfit`-placed C2 dimer:

```bash
sapia run linkpath <run_dir> \
    -t table1 \
    -i dimerfit_path \
    --from-chain auto --to-chain auto \
    --obstacle-chains all \
    --probe 2.0

sapia collect linkpath <run_dir> -t table1
```

`dimerfit_path` is the **transformed dimer without the target** — exactly two
protein chains, so `auto` resolves cleanly and `--obstacle-chains all` means "both
protomers and nothing else". **No `-g 0`**: the builder forces it.

To ask the *other* question — can the linker get around the target too? — run
`dimerfit_complex_path` with a label, so the two answers get their own columns:

```bash
sapia run linkpath <run_dir> -t table1 -i dimerfit_complex_path \
    -l withtarget --from-chain A --to-chain B --obstacle-chains all
```

(`-l withtarget` → `linkpath_withtarget_*`. Name the two protomer chains explicitly
there: with the target present the file has more than two protein chains and `auto`
is an **error**, by design.)

### Flags

| Flag | Default | What it does |
| --- | --- | --- |
| `-i` / `--input-column` | *(sentinel: none)* | The structure to route through. No honest default — see Traps. |
| `--from-chain` | `auto` | Chain donating the `--from-end` terminus. `auto` requires **exactly two protein chains** and takes the first; any other chain count is an **error for that design, never a guess**. |
| `--to-chain` | `auto` | Chain accepting the `--to-end` terminus. `auto` takes the **second** protein chain, same rule. |
| `--from-res` | *(empty)* | Residue **number** to anchor on instead of the chain's own terminus. Accepts the `{expr}` mini-language (`'{binder_len}'`, `'{binder_end - 2}'`), resolved per design up the lineage at **submit** time. |
| `--to-res` | *(empty)* | Likewise for `--to-chain`. |
| `--from-end` | `C` | Which terminus the route leaves from: `C` uses the backbone **C** atom (falling back to CA), `N` the backbone **N**. A fusion linker grows out of a C-terminus. |
| `--to-end` | `N` | Which terminus the route arrives at. |
| `--obstacle-chains` | `all` | **The scientifically load-bearing flag.** The chains whose atoms the route must go around, comma-joined. `all` = every chain in the file. A chain named but absent is an **error**, never a smaller selection. |
| `--include-h` | *off* | Treat hydrogens as obstacles. Off because most structures here are heavy-atom only and a mixed set would make designs incomparable. |
| `--protein-only` | *off* | Ignore ligands/cofactors/ions/nucleic acids. **Off by default: a cofactor in the way is in the way**, and turning this on can only shorten the route. |
| `--spacing` | `1.0` | Grid spacing (Å) — the precision/cost lever. Voxel count and runtime go as 1/spacing³. |
| `--probe` | **`2.0`** | Probe radius (Å): how fat the thing being routed is. **Differs from the prototype's 1.4** — see below. |
| `--pad` | `15.0` | Padding (Å) around the bounding box, so the route can leave the surface and go round the outside. |
| `--carve` | `3.0` | Bubble freed around each anchor so the chain can leave its own terminus. **The effective radius is raised to at least vdW(anchor) + probe** — see the trap. |
| `--taut-rise` | `3.5` | Å per residue of a fully extended chain → `n_res_min`. An absolute floor. |
| `--relaxed-rise` | `2.1` | Å per residue of a relaxed coil at ≤60 % of contour → `n_res_relaxed`. **Design against this one.** |
| `--no-both-directions` | *(both ON)* | Skip the reverse route. Leave it on: it is a **free trust metric**, one extra A\*. |
| `--no-residue-estimate` | *(estimate ON)* | Record distances only; `n_res_min`/`n_res_relaxed` come back **NA**. |
| `--model` | `0` | Model index for a multi-model file. |
| `--designs-per-task` | `50` | Designs per task. The work is ~0.16 s per design; this amortises the ~10 s container cold start. |

## The probe: why 2.0 and not the prototype's 1.4

`linker_path.py` defaults to **1.4 Å — a water probe**. That is the right radius for
a solvent-accessible surface and the wrong one for a polypeptide. A backbone is
thicker than a water molecule, so a 1.4 Å probe lets the route thread crevices a real
chain cannot enter, which **understates the detour**.

**2.0 Å is a deliberate campaign decision**: a chain-sized channel. The trade-off is
real in both directions:

| | 1.4 Å (water) | 2.0 Å (chain) |
| --- | --- | --- |
| Route through surface grooves | allowed | blocked |
| Measured on the C2 barnase dimer | **50.45 Å** | **52.80 Å** |
| Risk | too short — a route no chain can take | too long — a groove a thin loop could use is refused |
| Risk direction | **looks like a good design** | looks like a bad one |

Pick 2.0 because its failure mode is conservative. Use `--probe 1.4` only to
reproduce the prototype, and `--probe 3.0` if you want a bound that is safe for a
side-chain-bearing chain. **`probe` is echoed into the table** — a distance is
meaningless without it, and two runs at different probes are not comparable.

## `straight_ca_dist`: the bridge to dimerfit's retired `link_dist`

`dimerfit_link_dist` was **CA–CA** (C-term of protomer A to N-term of protomer B) and
**40 rows of it already exist**. The new tool reports two straight lines:

- `straight_dist` — **anchor atom to anchor atom** (backbone C → backbone N). The
  chemically correct endpoints for a linker, and the denominator of `detour_ratio`.
- `straight_ca_dist` — **CA to CA**. This exists *only* as the bridge: it is the
  direct equivalent of the retired `dimerfit_link_dist`, so old and new rows stay
  comparable.

They differ: **41.42 Å vs 41.76 Å** on the measured dimer. Compare historical
`dimerfit_link_dist` values to `linkpath_straight_ca_dist`, **never** to
`straight_dist` and **certainly never** to `path_dist`.

## The `path_found` contract — read this carefully

If A\* finds no free route, **the design is unlinkable at that probe radius. That is
a RESULT, not a failure.** So:

| | |
| --- | --- |
| `linkpath_status` | stays **`OK`** |
| `linkpath_path_found` | `False` |
| `path_dist`, `detour_ratio`, `n_res_min`, `n_res_relaxed` | **NA** |
| `note` | `no free path at probe 2 / spacing 1` |

**Do not read `status == "OK"` as "a route exists" — filter on `path_found`.** An
`error:` status here would conflate "this design cannot be linked" with "the tool
broke", and a standard trust filter would silently delete a real finding.

Genuine failures **do** get `error:` with every metric NA, recorded as data so a
partial array still collects and the batch exits 0:

- a chain named but absent; a chain with no amino acids;
- an unreadable or missing structure file;
- `auto` on a structure that is not exactly two protein chains;
- `--from-res`/`--to-res` naming a residue the chain does not have;
- a violated `path_dist >= straight_dist` invariant (a tool bug, asserted in-tool).

## `path_asymmetry`: a free trust metric under C2

With `--both-directions` (on by default) the tool also routes **C-terminus of
`to_chain` → N-terminus of `from_chain`**. Under **exact** C2 symmetry those two are
the same route carried by the symmetry operator, so they must be equal.
`path_asymmetry = |path_dist - rev_path_dist|` therefore measures the grid's own
discretisation error, for free, per design.

**Measured on an exact C2 dimer (barnase chain A + a true 180° copy), spacing 1.0,
probe 2.0: `path_asymmetry` = 0.162 Å** on a 52.8 Å route. So `< 1.0` is a sound
threshold at the default grid.

What a large asymmetry means — in order of likelihood:

1. **The structure is not actually C2.** A `rpxdock`/`dimerfit` dock should be; a
   predicted or relaxed one may not be. This is the useful signal.
2. **The grid is too coarse.** Raise `--spacing` resolution (lower the number) and
   see if it shrinks. Measured on the same dimer: 1.5 → 1.23 Å, 1.0 → 0.16 Å,
   0.75 → 0.50 Å, 0.5 → 0.35 Å.
3. **The two termini are genuinely in different environments** — true whenever the
   obstacle set is not symmetric (e.g. a target in the file), in which case the
   metric is not a trust metric at all and the threshold is meaningless.

On a **non**-symmetric pair the number is simply two different measurements: it read
**10.39 Å** on a deliberately asymmetric toy. Do not filter on it unless the input is
supposed to be symmetric.

## Traps

- **`-i` is effectively required.** `default_input_column` is the `"not applicable"`
  sentinel (as in `dimerfit`, `ifacegeom`, `cms`, `chainsel`): the builder refuses the
  run and lists the `_path` columns the table actually has. The structure column
  differs per campaign and a wrong one must not fail silently.
- **`auto` is not a guess.** On anything other than exactly two protein chains it is
  an `error:` for that design. Feeding `dimerfit_complex_path` (dimer **plus**
  target, 3+ chains) with `auto` is the obvious way to hit it — name the chains.
- **A narrower obstacle set only ever makes the answer SHORTER**, i.e. fails in the
  direction that looks like a good design. Read `obstacle_chains` before trusting a
  number; it is never the literal `all`, always what `all` resolved to.
- **`--protein-only` and `--include-h` are both measurement changes, not cosmetics.**
  `--protein-only` deletes obstacles; neither is echoed as its own column, so use
  `-l/--dir-label` if you run both ways.
- **The effective `--carve` is larger than you asked for, on purpose.** An anchor atom
  sits at the **centre** of its own bubble and blocks a sphere of radius vdW + probe
  around itself, so a `--carve` smaller than that wraps the bubble in a complete,
  unbroken blocked shell. With the campaign defaults (carve 3.0, probe 2.0, backbone
  C at 1.70) that shell is **0.70 Å thick around every anchor**, and the route could
  only escape by a 26-connected step jumping it. **Measured, before the fix**: the
  same dimer gave a 57.29 Å route at `--spacing 1.0` and **no path at all** at
  `--spacing 0.5`. linkpath therefore grows the bubble to `vdW(anchor) + probe`
  (3.70 Å for a backbone C at probe 2.0) and logs it to the task `.out`. After the
  fix the same structure gives 52.80 / 50.84 / 54.03 Å at spacing 1.0 / 0.75 / 0.5 —
  consistent, and ~4.5 Å shorter than the artifact. **Old prototype numbers are a
  safe upper bound, never an under-estimate** — see the next section for the proof.
- **`--spacing` is not just precision.** A coarse grid can cut corners between free
  voxel centres, which makes `path_dist` **too short**. `min_clearance` is the guard
  (below), and the residual spread across spacings is real: **50.8–54.0 Å across
  spacing 1.5–0.5** on the measured dimer, i.e. about ±3 Å / 6 %. Treat `path_dist`
  as good to a few Å, not to the 3 decimals it prints.
- **`--spacing 1.5` is the lever if it is ever too slow** (0.145 s vs 0.205 s on a
  324-residue trimer) — but it widens that spread. Prefer raising
  `--designs-per-task`.
- **`path_dist` includes the two anchor hops** — from each anchor atom to its first
  free voxel — exactly as the prototype's `total` does. It is anchor-atom to
  anchor-atom, directly comparable to `straight_dist`.
- **`direct_clear` is NOT comparable to the prototype's `direct path clear` line — do
  not reconcile the two.** They are different measurements with the same name, and
  the prototype's is useless: its grid-sampled `line_is_clear` reads
  `NO - obstructed` for **every** structure, including two residues alone in empty
  space (verified by running `linker_path.py` on exactly that), because the segment's
  first sample is the anchor atom's own voxel, which the anchor atom blocks.
  linkpath computes it **analytically** against every obstacle atom *except the two
  anchor residues' own*, which the segment necessarily starts and ends inside.
  `linkpath_direct_clear == True` ⟺ `detour_ratio` ≈ 1; a prototype "obstructed"
  means nothing at all.
- **Residue labels are `chain:resnum[icode]`** — `A:111`, matching `ifacegeom` and
  `dimerfit`. `--from-res`/`--to-res` take a bare **integer** (or `{expr}`), not a
  label.
- **`{expr}` resolves at SUBMIT time**, so a bad column name raises before anything is
  queued, naming the design. It cannot work on a root create (no lineage).
- **`-l/--dir-label` for variants.** Another structure column, another obstacle set,
  another probe → the columns collide otherwise. A `-l withtarget` run writes
  `linkpath_withtarget_*`.
- **A worker crash exits the task non-zero** (with fallback `error: worker crashed`
  TSVs); a *per-design* failure does not — the batch still exits 0 and the table
  carries the error.
- **Not yet ported to vib.** There is no `SAPIA_ACTIVATE_LINKPATH` anywhere. On Modal
  the image is `debian_slim + gemmi + numpy`, which builds in seconds.

## Are old `linker_path.py` numbers usable? Yes, as an upper bound

Because the sealed-anchor-bubble bug is a *finding about the prototype*, the question
matters for any number already written down from it: is the inflation **always** an
over-estimate, or can the blocked shell force a route out through a different gap and
come back **shorter** than the truth?

**It is always an over-estimate.** Both argued and measured.

*The argument.* Correcting the bug only ever **frees** voxels (`carve_eff ≥ carve`,
and nothing becomes blocked), so the corrected free set is a **superset** of the
buggy one. Every voxel path that is feasible on the buggy grid is therefore also
feasible on the corrected grid, with identical step costs. The start and goal voxels
are identical too — both are the rounded anchor voxel, which is inside the bubble and
free either way — so the two anchor hops are the same number. The corrected optimum
is a minimum over a superset of the same paths, hence **corrected ≤ buggy**. The
shell cannot "push" a route into a longer gap that the corrected run then has to use:
the corrected run can always still take the buggy run's route.

*The measurement.* **72 comparisons** — 6 structures × 4 spacings (1.5/1.0/0.75/0.5)
× 3 probes (1.4/2.0/2.6), same grid, carve 3.0 vs grown:

| Outcome | Count |
| --- | --- |
| Prototype inflated (corrected shorter) | 25 — up to **4.76 Å** |
| Identical | 20 |
| **Prototype false negative** ("no path" where a route exists) | 15 |
| Both report no path (a genuinely buried terminus) | 12 |
| **Corrected route LONGER than the buggy one** | **0** |

So, for a number already recorded from `linker_path.py`:

- **A reported length is a safe upper bound.** It is never shorter than what linkpath
  gives at the same probe and spacing. If the old number already cleared your linker
  budget, it still clears it.
- **A reported "No free path found" means nothing.** That is where the bug really
  bites: **15 of 72** configurations reported no route where one exists, and the
  effect gets *worse* as the grid gets finer or the probe larger (every
  `--spacing 0.5 --probe 2.0` case in the sweep). Never conclude a design is
  unlinkable from a prototype run — re-run it through the tool.
- **The inflation is not a constant you can subtract.** It ranged 0.00–4.76 Å and
  depends on the local geometry at the anchor.

## Columns collected (`linkpath_` prefix)

`<leaf>_status` is `OK`, `warn: …`, `error: …` or `missing`.
`<leaf>_path` is the route as a PDB of connected pseudo-atoms — **chain X = forward,
chain Y = reverse**. Load it alongside the structure in PyMOL and the detour is
visible directly. It is for looking at; nothing designs against it.

### Provenance / trust

A distance is meaningless without the parameters and the endpoints it was measured
between, so they are columns, not run notes.

| Column | Meaning |
| --- | --- |
| `from_res`, `to_res` | The anchor residues as `A:110` / `B:3`. **Check these against the construct's real termini** — the anchors are the last/first *resolved* amino acid. |
| `from_atom`, `to_atom` | The backbone atom actually used. Expect `C` and `N`. A **`CA`** means that residue had no C/N, so the distance is ~1.3 Å off the chemically right endpoint. |
| `from_chain`, `to_chain` | Resolved — `auto` never reaches the table. |
| `obstacle_chains` | **The auditable column**: the RESOLVED chain list the route had to go around. Never `all`. This is the answer to "was the target in there?". |
| `n_obstacle_atoms` | How many atoms that was. Halves when you drop a protomer (verified: 1728 → 864). |
| `probe`, `spacing` | Echoed, because they set the answer. Two runs at different values are not comparable. |
| `grid_shape` | `nx x ny x nz` voxels, e.g. `86x68x69`. |
| `direct_clear` | Was the straight anchor-to-anchor vector unobstructed? Analytic, excluding the two anchor residues' own atoms. `True` ⟺ `detour_ratio` ≈ 1. |

### Geometry

| Column | Meaning |
| --- | --- |
| `straight_dist` | Anchor atom to anchor atom (backbone C → backbone N). **A LOWER BOUND**, and the denominator of `detour_ratio`. |
| `straight_ca_dist` | CA to CA — **the bridge to the retired `dimerfit_link_dist`**. Compare old rows to this, not to `straight_dist`. |
| `path_dist` | The A\* route length, anchor to anchor, including the two short hops from each anchor atom to its first free voxel. **NA when no route exists.** |
| `detour_ratio` | `path_dist / straight_dist`. 1.0 = the straight line was free. Measured 1.00 (free space) to 1.28 (around a 110-aa body). |
| `path_found` | `False` = **no free route at this probe**. See the contract above — this is the column to filter on, not the status. |
| `min_clearance` | Closest approach of the route to an atom **surface** (Å), sampled continuously along the route, outside the anchor bubbles. A route that honestly stays in solvent-accessible space never comes closer than `probe`. Well below it means the grid cut a corner, which makes `path_dist` too **short**; below `probe - 0.25` the status becomes `warn:`. Measured 1.99–2.06 against probe 2.0. |
| `rev_path_dist` | The reverse route: C-terminus of `to_chain` → N-terminus of `from_chain`, each chain's **own** terminus (never `--from-res`/`--to-res`). |
| `path_asymmetry` | `|path_dist - rev_path_dist|`. **A trust metric under C2 only** — see above. |
| `n_res_min` | `ceil(path_dist / taut_rise)`. A taut, strained linker at contour length — an absolute floor, not a design. |
| `taut_rise` | The Å/residue constant that produced `n_res_min` (default 3.5). |
| `n_res_relaxed` | `ceil(path_dist / relaxed_rise)`. A relaxed coil at ≤60 % of contour. **Design against this one.** |
| `relaxed_rise` | The Å/residue constant that produced `n_res_relaxed` (default 2.1). |
| `note` | Free text: why no path was found, mostly. Empty on a clean row. |
| `seconds` | Wall time of the design's measurement. |

NA (empty) rather than 0 marks a field that did not apply. **A blank `path_dist` is
not "zero distance", it is "no route".**

**The two rise constants are columns for the same reason `probe` and `spacing` are:
a residue count cannot be re-derived from the table without them.** And the rise is
the *more* arbitrary of the two — 2.1 Å/res is a convention about how taut a coil may
be, not a physical constant, so a later reader will want to re-derive the count at a
different one. All four (`n_res_min`, `taut_rise`, `n_res_relaxed`, `relaxed_rise`)
are NA together under `--no-residue-estimate` and when no route was found: recording
a rise that was not applied to anything would be worse than recording nothing.

## Reading the numbers

- **`n_res_relaxed` is the answer.** `n_res_min` is a floor nobody should design to:
  a linker at contour length is a taut string, and it will pull the two protomers
  together rather than let them sit where the dock put them.
- **Both are lower bounds anyway.** The route hugs the surface; a real chain cannot.
  Budget above `n_res_relaxed`, not at it.
- **A `GGGGS` repeat is the usual realisation**: `ceil(n_res_relaxed / 5)` repeats.
  The tool does not emit that — it is one division and the repeat unit is a design
  choice.
- **`detour_ratio` tells you whether the question mattered.** At ≈ 1.0 the straight
  line was already right and `dimerfit`'s retired `link_dist` would have been fine.
  Well above 1 is the tool earning its place. **Context**: the retired straight-line
  `dimerfit_link_dist` measured 10.4–50.3 Å over 40 real docks and at a ~20 aa budget
  (≈60–70 Å fully extended) **never once bound** — and a straight line is only a
  lower bound, so the real routes are longer still. Expect `n_res_relaxed` to be the
  binding constraint, and size the budget before the pool, not after.
- **Compare designs only at identical `probe`, `spacing` and `obstacle_chains`.** All
  three are columns for exactly this reason.

### The filter the caller can now write

`-f` takes a **module** with `apply_filter(df) -> df`:

```python
def apply_filter(df):
    return df.query(
        'linkpath_status == "OK" '
        'and linkpath_path_found '
        'and linkpath_path_asymmetry < 1.0 '
        'and linkpath_n_res_relaxed <= 20'
    )
```

Clause by clause: `status == "OK"` drops failures **and** `warn:` rows whose route
clipped the protein; `path_found` drops the unlinkable designs that are deliberately
still `OK`; `path_asymmetry < 1.0` is the C2 trust gate (sound at the default grid,
measured 0.162 Å on an exact dimer) and must be dropped if the obstacle set is not
symmetric; `n_res_relaxed <= 20` is the science and the threshold is the campaign's
call. **Verified to run against a collected table.**

## Verification done

All of the below was run **locally** (repo `.venv`, gemmi 0.7.5 / numpy 2.5.3). The
tool has **not** yet been submitted through `sapia run` on Modal or vib; that first
run is the orchestrator's `-l verification`.

**Test structures** (built from PDB **1A2P**, barnase — 108 resolved aa per chain,
real coordinates):

- `c2_dimer.pdb` — chain A plus an **exact** 180° copy (true C2, 216 aa, 1728 heavy
  atoms, protomers in contact: min inter-chain distance 3.69 Å, 17 contacts < 5 Å).
- `opposite.pdb` — chain A with two single-residue chains planted on **opposite
  faces** along its principal axis, so the direct vector goes straight through.
- `free_pair.pdb` — two single-residue chains 30 Å apart in **empty space**.
- `buried.pdb` — chain A with a single-residue chain planted at its **centre of
  mass**: a genuinely buried terminus.

**1. Known answer, free space.** `path_dist` = `straight_dist` = **27.100 Å**
exactly (well inside one voxel diagonal, 1.732 Å), `detour_ratio` = 1.000,
`direct_clear` = True.

**2. Deliberate obstruction.** `detour_ratio` = **1.265**, `direct_clear` = False.
**Asserted numerically**: every written path vertex outside the anchor bubbles
(46 of 57) clears vdW + probe — worst vertex clearance **2.0061 Å** against probe
2.0. The continuous `min_clearance` column agrees at 1.993 Å (the ~0.007 Å dip is the
chord-vs-arc effect between two free voxel centres, well inside the 0.25 Å slack).

**3. The invariant.** `path_dist >= straight_dist` held on **45 runs** — 5 structures
× 3 spacings (1.5/1.0/0.75) × 3 probes (1.4/2.0/3.0). It is also asserted *inside*
the tool, and a violation is recorded as an `error:` row.

**4. Exact C2.** `path_asymmetry` = **0.162 Å** on a 52.80 Å route (forward 52.796,
reverse 52.634) at the defaults — comfortably under the 1.0 Å filter threshold.

**5. Deliberately bad inputs.** A batch of four (3-protein-chain file with `auto`;
nonexistent file; `--from-res 9999`; one good design): **the batch exited 0**, each
bad design got an `error:` status with **every metric NA**, and the good design still
collected `OK`.

**6. Biopython → gemmi port.** Ran `linker_path.py` (Biopython 1.88) and
`linkpath_worker.py` (gemmi) on the same files at the same settings, comparing the
route length at full precision. **7 of 8 comparisons are IDENTICAL to < 1.1e-6 Å**,
including both real-coordinate structures at both probes (c2_dimer: 50.449669 /
57.285618; opposite: 73.891729 / 74.577284). The grid shapes and straight-line
distances match exactly too. The single divergence is the synthetic `free_pair.pdb`,
whose atoms sit on **exact integer coordinates**: Biopython stores coordinates as
**float32** (1.45 → 1.45000005), which flips one knife-edge `d2 <= rr*rr` voxel test.
Nudging that structure 0.1 Å off the lattice makes the two agree to 1e-6. No real
structure is lattice-aligned.

**7. Two prototype defects found and fixed** (both verified by running
`linker_path.py` itself):

- Its `direct path clear` prints `NO - obstructed` even for `free_pair.pdb` — two
  residues alone in empty space. linkpath's analytic `direct_clear` reads `True`.
- Its free-space route comes out **27.93 Å** (probe 1.4) / **29.39 Å** (probe 2.0)
  for a problem whose exact answer is **27.10 Å**, because the anchor's own blocked
  shell seals its carve bubble. linkpath's grown carve gives 27.10 Å exactly.

**8. End to end through the CLI.** `sapia new_run` → a 4-row table → `sapia run
linkpath -t table0 -i struct_path` (manifest + sub-manifest written; `sbatch` is
absent locally so the submit itself raises, as expected) → `linkpath.sh` executed by
hand with the real prelude, `SAPIA_TASK_ID=1` → **task rc 0** → `sapia collect
linkpath`. **All 28 columns** (26 metrics + `linkpath_status` + `linkpath_path`)
landed with the right dtypes (booleans as real `True`/`False`, not strings), the
error row carried its message with every metric NA, the file-missing row collected as
`missing`, and the caller's filter expression ran against the result. Re-run after
the rise columns were added, and **both residue counts were re-derived from the
collected table alone** — `ceil(path_dist / taut_rise)` and
`ceil(path_dist / relaxed_rise)` reproduce `n_res_min` and `n_res_relaxed` on every
row, which is the entire point of carrying the constants.

**8b. The worker/collector column contract, checked mechanically.**
`tests/test_linkpath.py` asserts that the worker's `METRIC_COLUMNS` and the
collector's `RESULT_COLUMNS` are the same set, that the collector's four type
buckets partition them exactly once each, that the worker's TSV header is
`["name", "status", "path", *METRIC_COLUMNS]`, and that each rise sits immediately
after the count it produced. The two lists live in different files, so drift is
otherwise **silent** — a column added to the worker and not to the collector simply
never reaches the table, with no error anywhere. **The check was proven non-vacuous**
by deleting `relaxed_rise` from the collector: it failed with
`worker-only (written, never collected): ['relaxed_rise']`, and passed again on
restore. 8 tests, all passing.

**9. `{expr}`.** `--from-res '{binder_end - 2}'` with `binder_end = 112` wrote `110`
into the task file. `--from-res '{no_such_column}'` raised at **submit** time:
`ValueError: design 'c2': column 'no_such_column' not found (or empty) in lineage`.

**10. Timing, measured not guessed.**

| Structure | Atoms | Grid | spacing | Time/design |
| --- | --- | --- | --- | --- |
| 216-aa C2 dimer | 1728 | 86×68×69 | 1.0 | **0.161 s** (mean of 5; 0.153–0.172) |
| 324-aa trimer (1A2P, A→C) | 2627 | 101×88×96 | 1.0 | 0.205 s |
| 324-aa trimer | 2627 | 68×59×64 | 1.5 | 0.145 s |
| 324-aa trimer | 2627 | 134×117×127 | 0.75 | 0.345 s |

Both directions included. `--designs-per-task 50` is therefore ~8 s of compute per
task against a ~10 s cold start.

**11. The carve-bug direction.** 72 paired runs (6 structures × 4 spacings × 3
probes), buggy carve vs grown carve on the identical grid: the corrected route was
**never once longer**. See "Are old `linker_path.py` numbers usable?" above.

## Files

| File | Role |
| --- | --- |
| `tools/linkpath/spec.py` | The `Tool` descriptor (`action: update`, sentinel input column). |
| `tools/linkpath/run_linkpath.py` | Manifest builder: validates the flags, resolves `{expr}`, writes one sub-manifest per task, forces `gpus_per_task = 0`. |
| `tools/linkpath/linkpath.sh` | Per-task script, with the manifest-line guard copied from `dimerfit.sh`. |
| `tools/linkpath/linkpath_worker.py` | The measurement (gemmi + numpy). Runs standalone on a task file. |
| `tools/linkpath/collect_linkpath.py` | Per-design collector. |
| `tools/linkpath/modal_image.py` | `debian_slim + gemmi + numpy`, CPU-only. |
| `tests/test_linkpath.py` | The worker/collector column contract, the rise-column rules and the error-row contract. `.venv/bin/pytest tests/test_linkpath.py`. **Run it after touching either column list.** |
| `linker_path.py` (repo root) | The original Biopython prototype. Kept for reference and for the port comparison above; **not** used by the tool. |
