# pH-responsive EGFR binder via dimer protection — Phase 2 campaign plan

Companion to `dimer_binder.md`, which holds the strategy. This holds the execution shape: what
gets built, what runs in what order, and what each number is for.

History — what changed, when, and why — lives in
`campaigns/20261002_egfr_ph_dimer_binder.md`, not here. This document states the current plan.

---

## Context

A pH-responsive EGFR binder built by **dimer protection**: at pH 7.4 two copies self-associate
through a histidine-containing interface that occludes the EGFR-binding face; at low pH the
histidines protonate, the dimer opens, and the binder engages EGFR. The halves are finally
joined by a ~20 aa linker so the protection is intramolecular. Budget ≤250 aa total, ≤115 aa per
monomer. The trick that makes it tractable is **pH-agnostic design** — design only the neutral
state and let protonation break it.

**Phase 1 is skipped.** The campaign starts from BindCraft2 binders already designed against
both orthologs, in `sapia-runs-toon/inputs/bc2_output_for_anthony/`.

- **Backend:** Modal, volume `sapia-runs-toon` (`.env:28`).
- **Run_dir:** `outputs/20261002_143419_dimer_phase2`.
- **Tooling branch:** `worktree-ifacegeom` — `ifacegeom` exists only there, and the Modal
  workstation bakes `tools/` from the local working directory, so every run must be launched
  from that worktree until it is merged.

---

## 1. Status

| Step | Tool | State |
| --- | --- | --- |
| 0 | seed `table0` | **done** — 37 rows, hEGFR only |
| 1 | `ifacegeom` | **done** — 37/37 `OK` |
| 2 | independent revalidation | **skipped by decision** (§8) |
| 3 | `chainsel` | **done** — 37/37 `OK` |
| 4 | `rpxdock` C2 | **running** |
| 5 | `dimerfit` | **built, verified locally** — 40 docks, not yet submitted through `sapia run` |
| 6 | `linkpath` | **being built** |
| 7 | `hbdesigner` | not started |
| 8 | `graft` | **not built** |
| 9 | `atomium` / `proteinmpnn` | not started |
| 10 | `boltz` ×3 forks | not started |

`dimerfit`, `linkpath` and `graft` are built **just in time** — each is specced against real
output from the step before it, not against assumptions about what that output will look like.

---

## 2. Governing principles

### 2.1 Record, don't filter

Several of the strategy document's requirements are **desired outcomes, not admission
criteria**, and are treated as such until we have seen the distributions:

| Criterion | Status |
| --- | --- |
| Histidine-free epitope | Record `n_target_his`. **Unsatisfiable on this pool** — see §3. |
| C-term two-sided gate | Record `cterm_proj`, `cterm_iface_min_dist`. |
| Symmetric His H-bond network | Record network composition (§2.5). |
| ≤115 aa per monomer | Record `binder_len`. Hard only at the final linked construct. |
| Dimer-interface ⟷ epitope geometry | Record all of §4.1's columns. |
| Linker reach between protomers | Record `linkpath`'s distances and residue estimates (§4.2). The ≤20 aa linker is a hard requirement, but it binds the final construct, not a rigid dock. |

Inventing a threshold before seeing a distribution is how a campaign discards its only good
designs. Every one of these numbers lands in the table, carries through lineage, and can be
re-filtered at any cutoff later without recomputation. **Report distributions, not verdicts;
the thresholds are the user's call.**

### 2.2 The success criterion is never measured

The design must land in a *window* — the dimer must beat EGFR at pH 7.4 and lose to it at pH 6.
Nothing estimates either ΔG or its pH dependence, and no structure predictor is pH-aware, so
step 10 cannot test the switch either. Deferred by decision; recorded, not ignored.

The linker helps without being measured: intramolecular effective concentration is ~mM, so a
*weak* symmetric interface can hold the closed state.

### 2.3 Inheriting binders inverts the first job

We did not choose the epitope; we had to discover it. The upstream campaign had no reason to
respect the histidine, C-term or length requirements, so step 1 was **characterization of a
pool we did not design** — a picture of what we have, not a pass/fail.

### 2.4 "Close but not overlapping" and "strong clash" are different gates

The strategy wants the dimer interface near the epitope without full overlap, *and* strong
clashes with EGFR. Resolution: the dimer interface sits **adjacent** to the epitope while the
partner's **body** sterically covers it. Occlusion from the body high, residue overlap low —
separate columns, opposite-signed thresholds.

### 2.5 Symmetric His networks are constrained

His-N···H-donor paired with His-H···carboxylate is an *asymmetric* pair; under C2 each monomer
must supply both a histidine and a carboxylate. We record what hbdesigner finds rather than
demanding that motif up front — including networks with a single histidine, or none, which
still tell us what the interface can support.

### 2.6 Expect low yield, and rpxdock is biased against us

The target geometry is constrained at once by interface adjacency, occlusion, linker reach and
low residue overlap. Worse, rpxdock preferentially docks the C2 interface onto the **most
hydrophobic face — usually the EGFR face we are preserving**, so its top-scoring docks are
systematically the ones we want least.

Mitigation: carry the pool wide, `--nout-top 20` rather than the default 10, and treat rpxdock
`score` as **secondary to geometry**. With Phase 1 skipped the pool is fixed and finite — this
is the main risk to the campaign.

---

## 3. The pool

37 designs, hEGFR only (mouse deferred). Each design is shipped twice (`*_hEGFR.cif` +
`*_mEGFR.cif`), mmCIF only, with no accompanying metadata file.

| Fact | Value |
| --- | --- |
| **chain A** | the **target** — human EGFR 311–503 renumbered **1–193**, so `local = human − 310` |
| **chain B** | the **binder**, 60–118 aa, numbered from 1 |
| campaigns | v1 (5 designs), v2 (32) |
| inherited scores | embedded in each CIF as ~213 `_bindcraft.*` lines; both files of a pair carry both targets' scores |
| `renum/` | **ignored** — differs only in that the binder chain continues the target's numbering |

**`table0` is a hand-written input manifest.** There is no `sapia` import or seed verb (the
verbs are `new_run / init / fork-tool / modal-shell / run / collect`), so
`scripts/seed_table0_from_bc2.py` globs the top-level `*_hEGFR.cif` and writes `name`,
`input_path`, `source_file`, `ortholog`, `binder_len_src`. `binder_len_src` exists purely so
`ifacegeom_binder_len` has an independent number to be checked against.

### What step 1 measured

- **All five strategy hotspots (L325, P349, F412, V417, I467) are contacted by 37/37.** The
  pool hits exactly the patch it was aimed at.
- **The histidine-free epitope requirement is unsatisfiable.** His409 is contacted by 37/37 and
  sits **two residues from the F412 hotspot**; His346 by 29/37; His359 by 2. `n_target_his` is
  never 0. **Decision: accept and proceed pH-agnostic**, recording the column and excluding
  nothing. A histidine-free gate here would have emptied the pool.
- **Binder-side histidines in 34/37.** The binders carry histidines in their *own* epitope
  residues, which protonate at low pH in the same direction as His409. This bites at step 9:
  fixing the epitope positions during sequence fill **preserves** them. Not fixing them is a
  lever we hold and have not used.
- **The epitope footprint is large** — `n_binder_res` 15–35 (median 24) against binder lengths
  of 60–118, i.e. up to **46 % of the binder frozen** before atomium runs.
- **27/37 have `cterm_proj < 0`**, the desired side. The C-terminal tag constraint is far less
  binding on this pool than assumed and need not drive selection.

---

## 4. The tools this phase needs

Every per-design number becomes a **column**, never a side script. If the verdicts live in
files, re-exploring a threshold means recomputing everything — and no later reader can tell why
a row was carried forward.

Column names below are bare; the driver leaf-prefixes them and adds `<leaf>_status`.

### 4.1 `dimerfit` — built, action `update`

> *Premise: place a C2 dock back into the binder–target frame and measure whether the partner
> protomer occludes the target binding site.* Covers the epitope-adjacency and occlusion halves
> of the strategy's steps 3 and 4 in one tool, because both need the same superposition.

The linker question is **not** in this tool's scope and is no longer measured here — it is
`linkpath`'s job (§4.2). A straight terminal distance is only a lower bound, and the obstacle
set that makes the real number meaningful is the dimer alone, not the dimer plus target.

**Invoked** at step 5, on the dock child table.

| | |
| --- | --- |
| `-t` | the dock table |
| `-i` | `rpxdock_path` — the C2 assembly, dumped with `--use-orig-coords` |
| flags | `--ref-column` (the original complex) and `--epitope-column ifacegeom_binder_res`, both resolved from `table0` through lineage; `--binder-chain-in-ref`, `--target-chains-in-ref`; `--clash-cutoff` |

**What it does.** Kabsch-superposes dock chain A onto the reference binder, applies that one
transform to the **whole dimer**, writes the moved structure.

**Output** — two files plus columns:

- files: `path` (transformed dimer, no target) and `complex_path` (dimer + target, for
  inspection).
- mapping trust: `align_rmsd`, `seq_match_frac`, `resnum_offset` — these catch rpxdock's
  1..N renumbering silently shifting every later residue index.
- geometry: `dimer_iface_res`, `n_dimer_iface_res`, `epitope_com_dist`, `n_overlap_res`,
  `frac_overlap`.
- occlusion: `n_clash`, `clash_frac`, `occluded_frac`, `min_dist_b_target`.

`n_clash` / `clash_frac` are **heavy-atom distance counts, not `fa_rep`**. There is no Rosetta
in this image, and the strategy document's "fa_rep clash score" is not what gets computed.
They rank docks against each other; they are not an energy.

Built from `ringfit` as a template (same Kabsch → `gemmi.Transform` →
`transform_pos_and_adp` → write PDB shape, same `seq_match_frac` / `resnum_offset` sanity),
but a separate tool because the premise differs. **Verified locally on 40 real docks; not yet
submitted through `sapia run`** — that first Modal run is still the orchestrator's
`-l verification`.

### 4.2 `linkpath` — NEW, action `update`

> *Premise: measure the shortest route between two chain termini that stays out of the
> protein — the length a flexible linker would actually have to span.*

A standalone tool, not a dimerfit column, for three reasons. The premise is different
(routing through solvent, not occlusion of a site). The measurement is reusable — any
fusion-linker question in any later campaign asks it. And its answer depends entirely on
**which chains are treated as obstacles**, which is a decision the caller must make
explicitly rather than inherit from whatever file another tool happened to write.

**Invoked** at step 6, on the same dock table dimerfit annotated.

| | |
| --- | --- |
| `-t` | the dock table |
| `-i` | `dimerfit_path` — the **transformed dimer with no target** |
| flags | `--from-chain` / `--to-chain` (`auto` = chain 1's C-term → chain 2's N-term; any count other than 2 is an error, never a guess), `--obstacle-chains`, `--probe`, `--spacing`, `--pad`, `--carve`, `--no-both-directions`, `--taut-rise` / `--relaxed-rise`, `--no-residue-estimate`. CPU-only and **~0.16 s/design** — the builder forces `-g 0`, you don't pass it |

**What it does.** Builds a voxel grid over the structure, marks a voxel blocked if its centre
lies within (vdW radius + probe) of any obstacle atom, and runs A\* over the free voxels from
the C-terminal anchor to the N-terminal anchor with true Euclidean edge costs.

**The obstacle set is the whole scientific content, and here it is the dimer alone.** Not
because the target would corrupt the arithmetic — a straight distance is rigid-body invariant
and `dimerfit` never had the target in the file it measured — but because of *which state the
linker has to span*. At pH 7.4 the construct is closed and EGFR is **not bound**. Routing
around EGFR would invent an obstruction that does not exist at the moment the linker matters.
`linkpath` records the resolved chain list in `obstacle_chains` so this is auditable per row
rather than inferred from the invocation.

**Output** — one file plus columns:

- file: `path` — the route as PDB pseudo-atoms, to load beside the structure in PyMOL.
- provenance: `from_res`, `to_res`, `from_atom`, `to_atom`, `obstacle_chains`,
  `n_obstacle_atoms`, `probe`, `spacing`, `grid_shape`. A distance is meaningless without
  the probe radius it was measured at, so the parameters ride in the table beside it.
- geometry: `straight_dist` (backbone C → N anchors), `straight_ca_dist`, `path_dist`,
  `detour_ratio`, `direct_clear`, `path_found`.
- route trust: `min_clearance`, the closest approach of the route to an atom *surface*. Read
  it against `probe`. A route that cuts a corner between free voxel centres makes `path_dist`
  too **short** — wrong in the direction that flatters a design — and this is the only column
  that catches it. Below `probe − 0.25` the row is downgraded to `warn:`, which keeps every
  metric but drops it out of any `status == "OK"` selection. **Report the `warn:` count; do
  not let it vanish silently.**
- symmetry trust: `rev_path_dist` and `path_asymmetry`. Under exact C2 the two linkages are
  equivalent, so a non-zero asymmetry means the grid, the dock or the symmetry is wrong.
  Measured floor on an exact C2 dimer: **0.162 Å on a 52.8 Å route** at spacing 1.0 — that is
  grid discretization, so a threshold near 1 Å is the right order.
- budget: `n_res_min` (`path_dist`/`taut_rise`, taut and strained) and `n_res_relaxed`
  (`path_dist`/`relaxed_rise`, a relaxed coil at ≤60 % of contour length), each followed by
  the rise constant that produced it — `taut_rise` (3.5) and `relaxed_rise` (2.1) are
  columns, so both counts are re-derivable at a different convention from the table alone.
  All four are suppressed together by `--no-residue-estimate`.

`straight_ca_dist` exists for continuity: it is the direct equivalent of the retired
`dimerfit_link_dist`, so the 40 rows already measured stay comparable to anything new.

**"No path found" is a result, not a failure.** `status` stays `OK`, `path_found` is `False`,
the distance columns are NA and `note` says why. An `error:` status there would conflate *this
design cannot be linked* with *the tool broke*, and the standard `status == "OK"` trust clause
would then silently delete a real finding.

**The root-level prototype `linker_path.py` is superseded and was wrong.** Not stylistically —
numerically, in two places, both found while building this tool and both confirmed against a
known answer. Do not run it. Numbers already taken from it are salvageable **asymmetrically**:
a recorded **length is a safe upper bound** (proven, and 0 counterexamples in 72 paired runs —
so if it cleared the budget it still clears it), but a recorded **"no path found" is
worthless** (15 false negatives in those same 72 runs, and *every* fine-grid/large-probe case
among them). Details and signatures in `campaigns/20261002_egfr_ph_dimer_binder.md` §5.

**Expectation, stated before the numbers.** `dimerfit` measured the straight line at
10.4–50.3 Å over 40 docks — every one inside a 20 aa budget, so the constraint never bound.
I expect obstruction to change that: that a substantial minority of docks have a blocked
direct vector, and that among those with a straight line >30 Å the detour ratio exceeds ~1.3,
pushing some past 20 aa. **If instead `detour_ratio` comes back ≈1.0 on nearly every row, the
cheap metric was right all along and this tool should not be run at scale again.** Both
answers are cheap, and the second is worth knowing.

### 4.3 `graft` — NEW, action `update`

> *Premise: transplant named side chains, with their exact input rotamers, onto the
> structurally-equivalent positions of another structure — for both C2 copies.*

**Invoked** at step 8, on the hbdesigner child table — i.e. onto a poly-glycine dimer that
already carries a designed network. It runs **after** hbdesigner, not before: hbdesigner writes
`to_pdb(unk_to_gly=True)` and has no per-residue "do not design here" flag, so rotamers handed
to it would simply be erased.

| | |
| --- | --- |
| `-t` | the hbdesigner table |
| `-i` | `hbdesigner_path` |
| flags | `--ref-column` (via lineage), `--residues-column ifacegeom_binder_res`, `--network-column hbdesigner_network`, `--apply-symmetry`, `--on-collision network` |

**Output** — one file plus columns:

- file: `path` — the dimer with epitope rotamers restored on both copies.
- collision accounting: `n_network_res`, `n_epitope_res`, `n_overlap`, `n_grafted`, `n_skipped`.
- trust: `graft_bb_rmsd_max` — **the metric that decides whether the graft meant anything**,
  since a rotamer is only transferable where the backbone is near-identical; plus
  `n_clash_after`.
- hand-off: `fixed_positions` and `fix1..fixK`, the **union** of network and epitope positions,
  in the 1-based-per-chain convention `proteinmpnn --fixed-positions` expects.

Network residues win collisions and are counted in `n_skipped`.

### 4.4 `binder-campaign` skill — missing

Referenced by `README.md:112`, `.claude/agents/thinker.md:117` and `all-tools/SKILL.md` — and
does not exist. The thinker's instructions say to load it before reading any interface number.
Reconstruct from `campaigns/20260928_7ojg_slyb_binder.md`.

---

## 5. The chain

| Step | Tool | Reads | Writes to |
| --- | --- | --- | --- |
| 1 | `ifacegeom` | `table0` `input_path` | `table0` columns |
| 3 | `chainsel` | `table0` `input_path` | `table0` columns + monomer PDBs |
| 4 | `rpxdock` | `chainsel_path` | **child** dock table |
| 5 | `dimerfit` | `rpxdock_path` + `table0` via lineage | dock-table columns + PDBs |
| 6 | `linkpath` | `dimerfit_path` | dock-table columns + route PDBs |
| 7 | `hbdesigner` | `dimerfit_path` | **child** network table |
| 8 | `graft` | `hbdesigner_path` + lineage | network-table columns + PDBs |
| 9 | `atomium` / `proteinmpnn` | `graft_path`, `--fixed-positions` from `graft_fix*` | **child** sequence table |
| 10 | `mkcomplex` → `boltz` | sequence table, 3 `--table-label` forks | forked columns |

Steps 5 and 6 both annotate the **same** dock table — both are `update` tools, and `linkpath`
reads the file `dimerfit` wrote. Step 7 takes `dimerfit_path`, not `linkpath_path`:
`linkpath`'s output is a route of pseudo-atoms for inspection, never a design input.

**Step 1 — characterize.**
`ifacegeom -t table0 -i input_path --binder-chains B --target-chains A`.
Note the chain order: the binder is **B**. With the tool's defaults every number is wrong, and
the signature is `binder_len` coming back as 193.

**Step 3 — extract the monomer.**
`chainsel -t table0 -i input_path --chains B --rename-to A`, no renumber flag. Chain B is
already numbered 1..N, so this is an identity map `B:n → A:n` with nothing to record.
Invariants: `chainsel_n_res == ifacegeom_binder_len`, `chainsel_n_chains == 1`.

**Step 4 — dock.**
`rpxdock -t table0 -i chainsel_path --architecture C2 --nout-top 20 --hscore-files afilmv_ehl
--use-orig-coords --recenter-input --mem 64G`, no dir-label. All 37 at once.

**Step 5 — dock geometry.** `dimerfit` → report the joint distribution of `epitope_com_dist`,
`occluded_frac`, `frac_overlap`; choose cutoffs with the user. Measured on 40 docks,
`occluded_frac` and `frac_overlap` are **nearly collinear (r = 0.97)** and 0 of the 23 docks
with `occluded_frac > 0.5` had `frac_overlap ≤ 0.3` — so this is a trade-off to pick
deliberately, not two independent gates to AND together.

**Step 6 — linker reach.** `linkpath -i dimerfit_path -g 0` → report `path_dist`,
`detour_ratio`, `n_res_relaxed` and how many rows have `direct_clear == False`. Gate on
`path_asymmetry` first: it must be near zero on a true C2 dock, and is not if the geometry or
the grid is wrong.

**Step 7 — H-bond network.** `hbdesigner` on selected dimers, `--symm-chains A,B`, steered away
from the epitope with `--guide-res` / `--guide-radius`. Record network composition including
histidine content.

**Step 8 — transplant.** `graft` → `n_overlap` recorded, combined `fixed_positions` emitted.

**Step 9 — fill the sequence.** `atomium` / `proteinmpnn` with network ∪ epitope positions
fixed. (atomium has no `score` column, so its sequences cannot be ranked the way MPNN's can.)
**Open decision:** whether to also fix the binder-side epitope histidines (§3) or let them be
redesigned.

**Step 10 — validate.** `boltz` on three forked tables (`--table-label`): dimer alone,
EGFR+monomer, EGFR+dimer. Tests foldability and assembly, **not** the pH switch.

---

## 6. Hard tool constraints

Not scientific preferences — the ones that silently produce wrong numbers.

| Constraint | Consequence |
| --- | --- |
| rpxdock dumps **backbone-only** (N/CA/C/O/CB + `CEN`) by default | `--use-orig-coords` is **mandatory** — fa_rep clash, rotamer transplant and sequence remapping are all meaningless without real atoms. ~83 distinct atom names in a dump means it took; ~6 means it did not |
| rpxdock renumbers the dump **1..N per chain** | `dimerfit` re-derives the mapping by sequence, which is why real atoms are required |
| rpxdock needs an **origin-centred, single-chain** input | `chainsel` first; a cyclic dock **hard-errors** on >1 chain at submit time. Verify `n_chains_in == 1` and `input_com_dist` |
| rpxdock `--hscore-files afilmv_ehl` is a standing campaign decision | the default `ilv_h` is helix-only and returns **plausible, wrongly-scaled** numbers with no error. `rpxdock_hscore` on every row is the only post-hoc proof. A run that forgot the flag is discarded, not rescaled |
| `afilmv_ehl` needs **`--mem 64G`** | the 16G image default is not enough |
| rpxdock loads motif tables for **~200 s before docking starts** | a healthy task looks exactly like a hang for its first few minutes. The real hang is the same picture past ~10 min with a stall at `willutil/pdb/pdbfile.py:295` |
| rpxdock is a `create` tool | a scaffold that failed mints **no row at all**. Audit by comparing the child table's distinct parent names against the parent table |
| `linkpath`'s answer is decided entirely by `--obstacle-chains` | pointed at `dimerfit_complex_path` instead of `dimerfit_path` it routes around EGFR and silently inflates every distance, with no error. Read the resolved list back from `linkpath_obstacle_chains`; never infer it from the command you think you ran |
| a `linkpath` distance is meaningless without its probe radius | 1.4 Å (water) threads crevices a polypeptide backbone cannot physically enter and understates the detour; the campaign default is **2.0 Å**. `linkpath_probe` and `linkpath_spacing` ride in the table for exactly this reason |
| hbdesigner writes `to_pdb(unk_to_gly=True)`, no per-residue "don't design here" flag | hence graft runs after it; `--guide-res`/`--guide-radius` steer but do not forbid |
| hbdesigner `--symm-chains` is experimental; symmetrization failure writes **no PDB and exits 0**, recording `_hb0` | gate on `hbdesigner_status == 'OK'`, reject `_hb0` rows |
| bindcraft2's scores come from the **same AF2 that designed the binder** | inherited metrics are a floor, never evidence |
| mkcomplex `path` is always empty | sequence tool only |
| `modal-shell --cmd` **always exits 0** | never test `$?`. The authoritative submit evidence is `<out_dir>/<script>_logs/<script>_modal.json` holding `{app_id, n_tasks}` |

---

## 7. Verification

Each new tool is verified on **real data in the run_dir**, under a `verification` dir-label,
before use:

1. Run on a handful of designs; show the collected columns.
2. Check an invariant against a number the tool did **not** compute — a length from the parent
   table, a residue identity at a known position, a count of raw result files.
3. Confirm the trust metrics say the mapping was right (`seq_match_frac`, `resnum_offset`,
   `graft_bb_rmsd_max`, `align_rmsd`).
4. Confirm a deliberately bad input yields an error status, not a plausible number.

Specifically: `dimerfit`'s transform cross-checked by `usalign` of its output against the
reference; `graft` cross-checked by residue identity at grafted positions.

`linkpath` has two invariants that need no second tool. **`path_dist >= straight_dist` on
every row** — a violation is a broken grid, not a short route, and the tool errors rather
than reporting it. And **`path_asymmetry ≈ 0` on a true C2 dock**, because the two linkages
are symmetry-equivalent; a non-zero value means the dock, the chain assignment or the grid is
wrong. Beyond that, check the written route numerically: every pseudo-atom on the path must
sit at least (vdW + probe) from every obstacle atom.

---

## 8. Known gaps

- **No pH-responsiveness measurement anywhere** — the campaign's actual success criterion
  (§2.2). Closest future fix: PyRosetta rescoring with network histidines neutral vs.
  doubly-protonated.
- **The histidine-free epitope requirement is unsatisfiable on this pool** (§3). Accepted,
  unmeasured.
- **Binder-side histidines in 34/37** (§3). Unaddressed; the lever is at step 9.
- **The epitope footprint is up to 46 % of the binder** (§3), and nothing gates on that ratio.
- **Step 2 revalidation skipped.** Nothing independently confirms the inherited poses;
  BindCraft2's own scores come from the AF2 that designed them. If step 10 disappoints, this is
  the first place to look.
- **Phase 1 deferred**, so the pool is fixed and finite — the main yield risk (§2.6). If Phase 2
  exhausts it, Phase 1 is the unblock, and multi-target bindcraft2 would need a wrapper edit
  (`run_bindcraft2.py:438-446` hard-codes a single-element `targets` list).
- **Thresholds deliberately unset** (§2.1). They must be chosen and recorded before any claim
  about yield.
- **The linker is only ever measured in the closed state.** `linkpath` routes through the
  assembled dimer, which is the right obstacle set for pH 7.4 (§4.2) — but the open,
  EGFR-bound state has the protomers apart, and no structure for it exists until step 10.
  A linker that suffices closed may be strained open. The gap closes cheaply once step 10
  produces an EGFR+dimer prediction: rerun `linkpath` on that, with `--obstacle-chains` set
  to the binder chains only.
- **`linkpath` treats the structure as rigid.** Termini are the most mobile part of a
  protein and a real linker is flexible, so every number is a geometric floor, not a free
  energy. It says a route exists and how long it is, never that the chain will take it.
- **Mouse EGFR entirely unmeasured.** Whether His409 is conserved in mouse is a cheap sequence
  check that has not been done, and it interacts with the dual-species requirement.
- Provenance of the inherited binders is outside this repo.
- `CLAUDE.md`'s env table says `sapia-runs`; `.env` says `sapia-runs-toon`. Stale.
- rpxdock cage/dihedral architectures are out of scope; C2 is verified.
- No AF2 tool exists (relevant only if Phase 1 is revived).
