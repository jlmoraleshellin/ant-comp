# Campaign: pH-responsive EGFR binder via dimer protection

**Date:** 2026-10-02 · **Status: in progress — Phase 2, step 1**
**Backend:** Modal, volume `sapia-runs-toon` · **Branch:** `worktree-ifacegeom`
**Strategy doc:** `dimer_binder.md` · **Execution plan:** `dimer_binder_plan.md`

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
| Governing principle | **record, don't filter** | Plan §1.1. The pool is fixed and finite; inventing a threshold before seeing a distribution is how a campaign discards its only good designs. |
| Histidine in epitope | **accept, record, proceed pH-agnostic** | See §3 — the requirement is unsatisfiable on this pool. |
| Campaign run_dir | **clean mint**, label `dimer_phase2` | The earlier run_dir was a verification run and was deleted before this session. ifacegeom is CPU-only and batched, so re-running is nearly free and buys clean lineage from row zero. |
| Seeding `table0` | hand-written manifest via `scripts/seed_table0_from_bc2.py` | There is **no `sapia` import/seed verb** — confirmed, the verbs are `new_run/init/fork-tool/modal-shell/run/collect`. Plan §2.2 asks for this to be done once, explicitly, and recorded. This is the record. |

---

## 3. The histidine problem — the finding that reshapes the strategy

`dimer_binder.md` asks for an EGFR epitope **without histidines**, because their protonation
at low pH could ruin the very interface we need at low pH.

**On this pool that requirement is unsatisfiable.** From the (now deleted) verification run,
pending re-confirmation in `dimer_phase2`:

- **His409 is contacted by 37/37 binders.** It sits **two residues from the F412 hotspot** —
  the intended hydrophobic patch has a histidine built into it.
- **His346 is contacted by ~29/37.**
- **0/37 have a histidine-free target epitope.**

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

## 5. Traps found, with their signatures

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

1. Land the `dimer_phase2` ifacegeom run; confirm the `l<N>` invariant on all 37.
2. Report distributions of `cterm_proj`, `cterm_iface_min_dist`, `n_binder_res`,
   `n_target_his`; choose thresholds **with the user**, and record them here.
3. Update `ifacegeom/SKILL.md` to cite the live run_dir.
4. Decide the scope of step 2 revalidation before building further.
