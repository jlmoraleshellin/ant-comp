# Campaign: pH-responsive EGFR binder via dimer protection

**Date:** 2026-10-02 · **Status: in progress — Phase 2, step 5 (`dimerfit` being built)**
**Backend:** Modal, volume `sapia-runs-toon` · **Run_dir:** `outputs/20261002_143419_dimer_phase2`
**Branch:** `worktree-ifacegeom`
**Strategy doc:** `dimer_binder.md` · **Execution plan:** `dimer_binder_plan.md`

Steps done: 0 seed (37 rows) → 1 `ifacegeom` (37/37 OK) → 3 `chainsel` (37/37 OK) →
4 `rpxdock` C2 (37/37 OK, **740 docks** in `table1`). Step 2 skipped by decision.

**Brief:** design a binder that engages EGFR only at low pH. Two copies self-associate at
pH 7.4 through a histidine-containing interface that occludes the EGFR-binding face; at low
pH the histidines protonate, the dimer opens, and the binder engages EGFR. The two halves
are finally joined by a ~20 aa linker so the protection is intramolecular. Budget ≤250 aa
total, ≤115 aa per monomer. Must bind both human and mouse EGFR.

The trick that makes it tractable is **pH-agnostic design**: design only the neutral state
and let protonation break it.

**Phase 1 is skipped by decision.** The campaign starts from BindCraft2 binders already
designed against both orthologs, at `inputs/bc2_output_for_anthony/`. This log covers
**Phase 2 only**.

---

## 1. The inherited pool

Established by read-only inspection on `sapia-runs-toon`, 2026-10-02.

| quantity | value |
| --- | --- |
| designs | **37**, each shipped twice (`_hEGFR` + `_mEGFR`) = 74 top-level files |
| format | mmCIF only — no PDB, and **no shipped TSV/CSV/JSON** |
| chains | exactly A and B in every file; no HETATM records anywhere |
| **chain A** | **the target**, 193 res (human) / 194 res (mouse), numbered 1–193 |
| **chain B** | **the binder**, 60–118 res, numbered from 1 |
| campaigns | v1 (5 designs, hotspots A72/A105/A128/A107), v2 (32 designs, + A36/A99) |
| naming | `egfr_human_mouse_<v1\|v2>_multitarget_l<binder_len>_<16hex>_seq0_<hEGFR\|mEGFR>.cif` |

**The target sequence is constant.** All 37 chain-A sequences in `_hEGFR` files share one
md5; likewise the `_mEGFR` files share a different one. Each design's binder sequence is
identical between its h and m file — one multitarget design written out twice, target
swapped.

**Chain A is human EGFR 311–503 renumbered 1–193**, i.e. `local resnum = human − 310`.
Verified by residue identity at all five design-document hotspots:

| design doc (human) | local | residue |
| --- | --- | --- |
| L325 | A:15 | L |
| P349 | A:39 | P |
| F412 | A:102 | F |
| V417 | A:107 | V |
| I467 | A:157 | I |

**`binder_length_range` is a campaign setting, not a per-design value**, and the observed
lengths (60–118) are narrower than the ≤115 budget implies. Only the single `l118` design
exceeds it, so **length is very nearly a non-issue on this pool** — contrary to the plan's
assumption that it would be a live constraint.

### Ignore `renum/`

The `renum/` subdirectory holds 37 files, all `_hEGFR`, same names as the top-level ones.
The **only** difference is that the binder chain continues the target's numbering
(`B:194…` rather than `B:1…`); target numbering is identical. Per-chain numbering is what
`chainsel` and `rpxdock` expect downstream, so **the top-level files are the ones to use**.

### Inherited scores live inside the CIFs

There is no shipped metadata file. BindCraft2's scores are embedded as ~213 `_bindcraft.*`
key-value lines per file, and **both** files of a pair carry **both** targets' scores
(`i_pTM_hEGFR` and `i_pTM_mEGFR` alike). Keys appear per-target (`X_hEGFR` / `X_mEGFR`) and
per-validation-model (`X_<target>_model_1_ptm`). There is **no `rank` and no
`Binder_Length` key** — binder length exists only in the filename and the chain-B count.

These are a **floor, never evidence**: they come from the same AF2 that designed the binders.

---

## 2. Decisions taken, and why

| Decision | Choice | Reasoning |
| --- | --- | --- |
| Backend | Modal, `sapia-runs-toon` | `.env:28`. `CLAUDE.md`'s env table still says `sapia-runs` — stale. |
| Ortholog order | **human first**, mouse later | User's call. Mouse is a second `ifacegeom` run (own rows, or `-l` on a second table). |
| `renum/` | **excluded** | See above. |
| Governing principle | **record, don't filter** | Plan §2.1. The pool is fixed and finite; inventing a threshold before seeing a distribution is how a campaign discards its only good designs. |
| Histidine in epitope | **accept, record, proceed pH-agnostic** | See §3 — the requirement is unsatisfiable on this pool. |
| Campaign run_dir | **clean mint**, label `dimer_phase2` | The earlier run_dir was a verification run and was deleted before this session. ifacegeom is CPU-only and batched, so re-running is nearly free and buys clean lineage from row zero. |
| Seeding `table0` | hand-written manifest via `scripts/seed_table0_from_bc2.py` | There is **no `sapia` import/seed verb** — confirmed, the verbs are `new_run/init/fork-tool/modal-shell/run/collect`. Plan §3 asks for this to be done once, explicitly, and recorded. This is the record. |
| Linker metric | **obstruction-aware path, not a straight line** | A terminal-Cα chord is only a lower bound on what a linker must span. Measured, `dimerfit_link_dist` was 10.4–50.3 Å over 40 docks — always inside a 20 aa budget, i.e. a gate that never fired. Replaced by an A\* route through solvent-accessible space. |
| Where that metric lives | **new standalone tool `linkpath`**; `link_dist` and `cterm_to_nterm_res` **removed from `dimerfit`** | Different premise (routing vs. occlusion), reusable for any future fusion-linker question, and its answer is set by the obstacle set — a decision the caller must make explicitly. User's call to remove rather than keep the chord as a pre-screen. |
| `linkpath` obstacle set | **the dimer alone** (`dimerfit_path`, not `dimerfit_complex_path`) | Not because the target corrupted the old number — it never did; `link_dist` was computed on a dimer-only model and a chord is rigid-body invariant anyway. Because of *which state the linker spans*: at pH 7.4 the construct is closed and EGFR is not bound, so routing around EGFR would invent an absent obstruction. |
| `linkpath` probe radius | **2.0 Å**, not 1.4 | A polypeptide backbone is thicker than water; a water probe threads crevices a real chain cannot enter and understates the detour. Recorded per row as `linkpath_probe`. |

---

## 3. The histidine problem — the finding that reshapes the strategy

`dimer_binder.md` asks for an EGFR epitope **without histidines**, because their protonation
at low pH could ruin the very interface we need at low pH.

**On this pool that requirement is unsatisfiable.** Confirmed in `dimer_phase2`, 37/37 `OK`
(an earlier verification run had found the same numbers and was deleted before this session;
these are the re-measured ones):

- **His409 is contacted by 37/37 binders.** It sits **two residues from the F412 hotspot** —
  the intended hydrophobic patch has a histidine built into it.
- **His346 is contacted by 29/37**; His359 by 2.
- **0/37 have a histidine-free target epitope.** `n_target_his` is 1 in 8 designs, 2 in 27,
  3 in 2 — never 0.

This is the "record, don't filter" principle earning its keep: a histidine-free gate applied
at step 1 would have **emptied the pool**.

**Decision: accept it and proceed pH-agnostic** (user, 2026-10-02). `n_target_his` stays a
recorded column; nothing is excluded. The risk is explicit and currently unmeasured: the
EGFR interface may *also* weaken at low pH, partly cancelling the switch.

Options considered and **not** taken, recorded so they need not be re-derived:

- **Exploit it** — place a carboxylate/acceptor near His409 so protonation at low pH
  *strengthens* binding, widening the switch window. Costs the pH-agnostic simplification on
  the EGFR face and needs pH-aware scoring we do not have (§6).
- **Prefer low-His designs** — prioritise the ~8 contacting only His409, not His346.
  Still available as a cheap way to manage rpxdock cost by sampling.
- **Re-epitope via Phase 1** — the only real unblock; expensive, and needs the bindcraft2
  multi-target wrapper edit (`run_bindcraft2.py:438-446` hard-codes a single-element
  `targets` list).

**Open and worth one cheap query:** is His409 a histidine in **mouse** EGFR too? If the
orthologs differ there, the dual-species requirement and the pH behaviour interact.

---

## 4. The pipeline as planned

`dimer_binder_plan.md` §3.6 holds the full table. Status as of this log:

| Step | Tool | State |
| --- | --- | --- |
| 0 | seed `table0` | script written; running |
| 1 | **`ifacegeom`** (new) | **built, verified, running on the campaign run_dir** |
| 2 | `cms`, `pyrosetta`, `boltz` | not started — independent revalidation of inherited binders |
| 3 | **`chainsel`** (needs edit: record `resnum_map`) | not started |
| 4 | `rpxdock` C2 | not started — **needs a single-design timing probe first** |
| 5 | **`dimerfit`** (new) | not built |
| 6 | `hbdesigner` | not started |
| 7 | **`graft`** (new) | not built |
| 8 | `atomium` / `proteinmpnn` | not started |
| 9 | `boltz` ×3 forks | not started |

Build order is **just-in-time**: `dimerfit` is specced once we have seen a real `rpxdock`
dump, `graft` once we have seen real hbdesigner output. Specifying them against assumptions
is how they get the chain mapping wrong.

---

## 4b. Results so far

### Step 1 — `ifacegeom` on 37 hEGFR complexes

Invariant that validates everything else: **`ifacegeom_binder_len` matched the filename's
`l<N>` on all 37**, a number the tool never reads. Had the chain assignment been wrong it
would have returned 193 (the target) on every row.

| | min | Q1 | median | Q3 | max |
| --- | --- | --- | --- | --- | --- |
| `binder_len` | 60 | 64 | 74 | 86 | 118 |
| `n_binder_res` (epitope footprint) | 15 | 21 | 24 | 28 | 35 |
| `n_target_res` | 22 | 28 | 30 | 33 | 45 |
| `cterm_proj` | −18.1 | −7.0 | −4.7 | 0.1 | 16.4 |
| `cterm_iface_min_dist` | 0.0 | 5.2 | 8.7 | 12.2 | 21.9 |

- **All five hotspots contacted by 37/37** (L325, P349, F412, V417, I467), residue identity
  confirmed independently from the raw CIFs. The pool hits the patch it was aimed at — but this
  is uniform, so it carries **no discriminating information**.
- **27/37 have `cterm_proj < 0`** (C-terminus on the far side from the epitope, the desired
  side). **A prediction was recorded that most would be positive, and it was wrong** — the
  C-terminal tag constraint is far less binding on this pool than the plan assumed and need not
  drive selection.
- **`n_binder_his > 0` in 34/37** — binder-side epitope histidines, anticipated nowhere in
  either planning document. Only 3 designs have none.
- 36/37 are ≤115 aa; only `v1_l118` exceeds the budget.

Cross-tabulating the two recorded criteria, 5 designs carry **only** His409 *and* have the
C-terminus on the far side — `v2_l73_e306ad4b` (the standout: 15-residue footprint,
`cterm_iface_min_dist` 21.9), `v2_l73_48b9cdd4`, `v2_l74_f6ae3fc9`, `v2_l61_d94f18a9`,
`v2_l86_adb0cabb`. **Recorded, not used as a filter.**

### Steps 3–4 — `chainsel` → `rpxdock` C2

`chainsel --chains B --rename-to A`, no renumber flag: 37/37 `OK`,
`chainsel_n_res == ifacegeom_binder_len` on all 37, `chainsel_n_chains == 1` on all 37.

`rpxdock` C2, `--nout-top 20 --hscore-files afilmv_ehl --use-orig-coords --recenter-input
--mem 64G`: 37/37 `OK`, **740 rows** in `table1`, 740 dock PDBs on disk.

| | min | Q1 | median | Q3 | max |
| --- | --- | --- | --- | --- | --- |
| `score` | 63.6 | 100.3 | 115.8 | 131.3 | 192.8 |
| `rpx` | 63.1 | 99.3 | 114.9 | 129.9 | 191.6 |
| `ncontact` | 24 | 78 | 95 | 118 | 286 |
| `n_docks` | 55 | 74 | 81 | 99 | 121 |
| `hscore_seconds` | 167.0 | 194.3 | 197.1 | 207.8 | 250.9 |

Trust: `hscore` = `afilmv_ehl` on 740/740; `n_chains_in` = 1 on 740/740; `recentered` True on
all; `input_com_dist` 85.5–96.6 Å, so **`--recenter-input` was essential, not decorative**.
`--use-orig-coords` took: 8.25 atoms/residue, **0 `CEN`**. Lineage audit: 37 distinct parents,
exact match, no scaffold dropped.

> **`--nout-top 20` bound on ALL 37 scaffolds.** `n_docks` is 55–121, so **~3,000 poses were
> discarded against 740 kept.** This matters because rpxdock is expected to favour docking onto
> the hydrophobic face we are preserving — i.e. the useful geometry may sit *below* the score
> cut. Re-dumping is not cheap: `_Result.txz` cannot be reopened in this image (no PyRosetta),
> so more poses means paying the ~200 s table load again per scaffold.
> **Decision deferred to `dimerfit`:** measure whether the wanted geometry correlates with
> `score` rank before buying more poses. **Prediction on record: it will correlate negatively
> or not at all.**

---

## 5. Traps found, with their signatures

- **`linker_path.py` (the root-level prototype) seals its own anchors inside a blocked shell.**
  The anchor atom sits at the centre of the bubble `--carve` frees, but the atom still blocks a
  sphere of `vdW + probe` around itself. With the defaults in play (`carve` 3.0, `probe` 2.0,
  backbone C `vdW` 1.70) the freed bubble is wrapped in a **complete 0.70 Å-thick blocked
  shell**, so the route can only escape through a grid artifact. **Signature:** the answer moves
  with `--spacing` and gets *worse* as the grid gets finer — measured on one C2 dimer, **57.29 Å
  at spacing 1.0 and NO PATH AT ALL at spacing 0.5**. A path-length that is not stable under
  refinement is the tell. **Proof it is a bug and not a modelling choice:** on a free-space pair
  whose answer is known exactly to be 27.10 Å, the prototype returns **27.93 Å (probe 1.4) /
  29.39 Å (probe 2.0)**; `linkpath`, which grows the bubble to `vdW(anchor) + probe`, returns
  **27.10 Å** exactly. On the real dimer the fix shortens the route by ~4.5 Å and makes it
  stable across spacings (52.80 / 50.84 / 54.03 Å at 1.0 / 0.75 / 0.5).
  **Are old prototype numbers salvageable? Lengths yes, verdicts no.** Correcting the defect
  only ever *frees* voxels, so the corrected free set is a strict superset of the buggy one
  and the start/goal voxels are identical — the corrected optimum is a minimum over a
  superset of the same paths, hence **corrected ≤ buggy necessarily**. Measured over 72
  paired runs (6 structures × 4 spacings × 3 probes): 25 inflated by up to **4.76 Å**, 20
  identical, **0 cases where the corrected route came back longer**. So a recorded length is
  a safe **upper bound** — if it cleared the linker budget it still clears it; if it failed,
  it may be a false reject by up to ~5 Å. The inflation is **not a constant you can
  subtract** (0.00–4.76 Å, set by local geometry at the anchor). But a recorded **"No free
  path found" is worthless**: 15 of those 72 runs reported no route where one exists, and it
  worsens as the grid gets finer or the probe larger — **every** `--spacing 0.5 --probe 2.0`
  case in the sweep was a false negative. Never conclude a design is unlinkable from a
  prototype run.
- **`linker_path.py`'s `direct path clear` line is always "obstructed".** Its grid-sampled
  `line_is_clear` walks the straight segment including the two anchor residues' own atoms,
  which the segment necessarily starts and ends inside. **Signature:** it reports obstructed
  for *two residues alone in empty space* — verified. `linkpath_direct_clear` is computed
  analytically and excludes the anchor residues' own atoms, so **the two are not comparable
  and must never be reconciled.**
- **A voxel route can cut corners between free voxel centres**, which makes the path come back
  too **short** — wrong in the direction that flatters a design, and therefore the dangerous
  direction. **Signature:** none visible in `path_dist` itself; this is why `linkpath` carries
  `min_clearance` (closest approach of the continuously-sampled route to an atom surface) and
  downgrades a row to `warn:` below `probe − 0.25`. Measured 1.99–2.06 Å against probe 2.0 on
  healthy rows.

- **The binder is chain B, the target chain A** — the *opposite* of `ifacegeom`'s default and
  of the plan's own §3.1 example. **Signature if missed:** `ifacegeom_binder_len` comes back
  as 193 (the target) on every row. The independent check is the filename's `l<N>`, which the
  tool never reads.
- **"`--contact-cutoff` 8 Å to match InterGroupInterfaceByVector" is wrong**, and was in the
  plan. 8 Å is that selector's **CB–CB** cutoff, not an atom–atom one; 8 Å of heavy-atom
  separation is not a contact. Setting it selects most of the binder. Default is **5.5 Å**
  (Rosetta's `nearby_atom_cut`) with the CB–CB cutoffs exposed separately.
- **A tool's `SKILL.md` can cite a run_dir that no longer exists.** `ifacegeom`'s "Verified
  on real data (Modal)" section cites `outputs/20261002_140341_ifacegeom_hegfr_test`, which
  was deleted before this session. **Signature:** the orchestrator reports "No such file or
  directory" and `outputs/` is empty on every volume. The verification was genuine; the
  evidence is not re-readable. **Lesson: a commit message is not a measurement.** Numbers
  inherited from prose must be re-established before anything is decided on them — and a
  campaign decision was very nearly taken on these.
- **`modal-shell --cmd` always exits 0.** Never test `$?`. The authoritative submit evidence
  is `<out_dir>/<script>_logs/<script>_modal.json` holding `{app_id, n_tasks}`; no file means
  nothing was queued.
- **An empty manifest line makes a task run on nothing, silently.** The single most expensive
  failure of this session: **18 of 37 rpxdock tasks** died with
  `IsADirectoryError: [Errno 21] Is a directory: '.'`. Cause: the prelude does
  `SAPIA_LINE="$(sed -n "${SAPIA_TASK_ID}p" "$MANIFEST")"`, and **`sed -n Np` on a file that
  does not yet show N lines exits 0 with empty output** — `set -e` does not trip. Three empty
  fields reached the worker and `--config ''` became `Path('') == '.'`.
  **Signature:** a `.out` reading `task N: rpxdock on  ()` with the name *and* path blank;
  optionally `warning: command substitution: ignored null byte in input` from the prelude.
  **What it was NOT:** the manifest was provably intact afterwards (37 lines, 3 fields, zero
  NUL bytes, byte-identical in structure between the 18 that failed and the 19 that passed);
  versions matched across the boundary (sha256-identical prelude/`modal.py`/`sapia_modal_task`);
  off-by-one was ruled out. It is a **read-side visibility failure in the task container**.
  Note prosapia already defends against this — `publish_manifest` uses `batch_upload` and its
  docstring names this exact failure — **and it was not sufficient.** One task re-read ~2
  minutes later and still saw nothing, so the stale view can outlive a short backoff; warm
  container reuse carrying an older Volume mount is the leading hypothesis, **unproven**.
  **Fix:** commit `30566b5` hardens `tools/rpxdock/rpxdock.sh` — retry 5× with backoff logging
  the line count the container sees, then fail loudly; check the **field count** before the
  empty check (because `cut -f2` on a line with no tabs returns the *whole line*, so a
  one-field line would pass as `name == input == config`); check each path is a non-empty file.
  **The hazard is workspace-wide** — every tool's `.sh` reads the manifest this way.
  `ifacegeom` and `chainsel` got lucky, not immune. The real fix belongs upstream in
  `sapia_task_prelude.sh`.
- **The rpxdock skill's "~83 distinct atom names proves `--use-orig-coords` took" is wrong.**
  83 implies ~15.8 atoms/residue, which only works if hydrogens were counted; ~35 is the
  correct heavy-atom name count and is **not diagnostic**. The real test is
  **~8.4 atoms per residue and zero `CEN` pseudo-atoms**. Backbone-only is ~5/residue + `CEN`.
- **A deleted run_dir leaves a tool's `SKILL.md` citing evidence nobody can re-read.** Harmless
  here (the run was deliberately cleaned up and the numbers reproduced exactly), but a campaign
  decision was very nearly taken on figures that could not be verified. **A commit message is
  not a measurement.**
- **Custom tools are baked from the local working directory.** `ifacegeom` lives only on
  `worktree-ifacegeom`, so the workstation must be launched from that worktree. From the main
  checkout the tool is simply absent from `sapia run --help` — it looks like it was never
  written rather than like a path problem.
- **Residue-label conventions differ between tools.** `ifacegeom` emits `A:12`;
  `ringfit --hotspots` wants `A47`. Do not paste one into the other.

---

## 6. Known gaps

- **Nothing measures pH responsiveness** — the campaign's actual success criterion. The design
  must land in a *window*: the dimer beats EGFR at pH 7.4 and loses to it at pH 6. Nothing
  estimates either ΔG or its pH dependence, and **no structure predictor is pH-aware**, so the
  step-9 Boltz validation tests foldability and assembly, not the switch. Closest future fix:
  PyRosetta rescoring with network histidines neutral vs. doubly-protonated. *Deferred by
  decision, not by oversight.*
- **His409 is in every epitope** (§3). Accepted, unmeasured.
- **Phase 1 deferred**, so the pool is fixed, finite and the main yield risk. rpxdock is
  additionally **biased against us**: it preferentially docks the C2 interface onto the most
  hydrophobic face — usually the EGFR face we are preserving — so its top-scoring docks are
  systematically the ones we want least.
- **Thresholds deliberately unset.** They must be chosen and recorded before any claim about
  yield, or the table cannot say why a row was carried forward.
- Provenance of the inherited binders is outside this repo; what we know of them is only what
  step 1 measures.
- `binder-campaign` skill is referenced by `README.md:112`, `.claude/agents/thinker.md:117`
  and `all-tools/SKILL.md` — and **does not exist**.
- hbdesigner `--symm-chains` is experimental upstream; symmetrization failure writes **no PDB
  and exits 0**, recording a `_hb0` row.
- Mouse EGFR entirely unmeasured so far.

---

## 7. Next

1. **Verify `dimerfit`** on real data under a `verification` dir-label, then run it on all 740
   docks. Checks that must pass before any number is read: `seq_match_frac` ≈ 1.0,
   `resnum_offset` == 0, per-protomer `n_res` == the parent's `ifacegeom_binder_len`, and a bad
   input giving `error:` with NA rather than a plausible number.
2. **Settle the `--nout-top` question** (see §4b): does the wanted geometry correlate with
   `rpxdock_score` rank? If the good docks cluster at ranks 15–20 we lost poses and should
   re-dock wider; if they are at ranks 1–5, the cut cost nothing.
3. **Choose thresholds with the user** for `occluded_frac` (want high), `frac_overlap` (want
   low), `epitope_com_dist`, and `linkpath`'s `n_res_relaxed` against the ~20 aa budget.
   Record them here — without them the table cannot say why a row was carried forward.
   Note `occluded_frac` and `frac_overlap` measured r = 0.97, so those two are one
   trade-off dial, not two independent gates.
4. **Decide step 8's histidine policy:** fix the binder-side epitope histidines (34/37) or let
   atomium redesign them.
5. **Cheap and not done:** is His409 a histidine in mouse EGFR? One sequence vs one sequence.

### Housekeeping debt

- **`uv.lock` on this branch pins prosapia `post18`**, whose Modal executor raises
  `commit() can only be called on a mounted volume inside a container` and **queues nothing
  while printing `Submitting N designs`**. The worktree venv was repaired in place to `post19`
  / `4782f10`; the lockfile was not. The next person to build this branch gets the broken pin.
- `ifacegeom/SKILL.md` cites the deleted run_dir; repoint it at
  `outputs/20261002_143419_dimer_phase2`.
- `rpxdock/SKILL.md` carries the wrong atom-name diagnostic (see §5).
- `worktree-ifacegeom` is unmerged — `dimer_tools_dev` has no `ifacegeom` and no `dimerfit`.
  A second stale worktree, `worktree-tool-commits`, is also unmerged.
- `binder-campaign` skill is still missing while three files reference it.
