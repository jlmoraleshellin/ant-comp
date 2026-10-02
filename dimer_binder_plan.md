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
| 2 | independent revalidation | **skipped by decision** (§7) |
| 3 | `chainsel` | **done** — 37/37 `OK` |
| 4 | `rpxdock` C2 | **running** |
| 5 | `dimerfit` | **not built** |
| 6 | `hbdesigner` | not started |
| 7 | `graft` | **not built** |
| 8 | `atomium` / `proteinmpnn` | not started |
| 9 | `boltz` ×3 forks | not started |

`dimerfit` and `graft` are built **just in time** — each is specced against real output from the
step before it, not against assumptions about what that output will look like.

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

Inventing a threshold before seeing a distribution is how a campaign discards its only good
designs. Every one of these numbers lands in the table, carries through lineage, and can be
re-filtered at any cutoff later without recomputation. **Report distributions, not verdicts;
the thresholds are the user's call.**

### 2.2 The success criterion is never measured

The design must land in a *window* — the dimer must beat EGFR at pH 7.4 and lose to it at pH 6.
Nothing estimates either ΔG or its pH dependence, and no structure predictor is pH-aware, so
step 9 cannot test the switch either. Deferred by decision; recorded, not ignored.

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
  residues, which protonate at low pH in the same direction as His409. This bites at step 8:
  fixing the epitope positions during sequence fill **preserves** them. Not fixing them is a
  lever we hold and have not used.
- **The epitope footprint is large** — `n_binder_res` 15–35 (median 24) against binder lengths
  of 60–118, i.e. up to **46 % of the binder frozen** before atomium runs.
- **27/37 have `cterm_proj < 0`**, the desired side. The C-terminal tag constraint is far less
  binding on this pool than assumed and need not drive selection.

---

## 4. Tools still to build

Every per-design number becomes a **column**, never a side script. If the verdicts live in
files, re-exploring a threshold means recomputing everything — and no later reader can tell why
a row was carried forward.

Column names below are bare; the driver leaf-prefixes them and adds `<leaf>_status`.

### 4.1 `dimerfit` — NEW, action `update`

> *Premise: place a C2 dock back into the binder–target frame and measure whether it occludes
> the binding site and can be linked.* Covers the strategy's steps 3 and 4 in one tool, because
> both need the same superposition.

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
- geometry: `dimer_iface_res`, `n_dimer_iface_res`, `iface_com_dist`, `n_overlap_res`,
  `frac_overlap`, `link_dist` (Cα C-term A → N-term B).
- occlusion: `n_clash`, `clash_fa_rep`, `occluded_frac`.

*Build from `ringfit` as a template* — `ringfit_worker.py` already does Kabsch →
`gemmi.Transform` → `transform_pos_and_adp` on a whole model → write PDB, with
`seq_match_frac` / `resnum_offset` sanity. The premise differs, so this is a **new tool**, not
a fork.

### 4.2 `graft` — NEW, action `update`

> *Premise: transplant named side chains, with their exact input rotamers, onto the
> structurally-equivalent positions of another structure — for both C2 copies.*

**Invoked** at step 7, on the hbdesigner child table — i.e. onto a poly-glycine dimer that
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

### 4.3 `binder-campaign` skill — missing

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
| 6 | `hbdesigner` | `dimerfit_path` | **child** network table |
| 7 | `graft` | `hbdesigner_path` + lineage | network-table columns + PDBs |
| 8 | `atomium` / `proteinmpnn` | `graft_path`, `--fixed-positions` from `graft_fix*` | **child** sequence table |
| 9 | `mkcomplex` → `boltz` | sequence table, 3 `--table-label` forks | forked columns |

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

**Step 5 — dock geometry.** `dimerfit` → report the joint distribution of `iface_com_dist`,
`link_dist`, `occluded_frac`, `frac_overlap`; choose cutoffs with the user.

**Step 6 — H-bond network.** `hbdesigner` on selected dimers, `--symm-chains A,B`, steered away
from the epitope with `--guide-res` / `--guide-radius`. Record network composition including
histidine content.

**Step 7 — transplant.** `graft` → `n_overlap` recorded, combined `fixed_positions` emitted.

**Step 8 — fill the sequence.** `atomium` / `proteinmpnn` with network ∪ epitope positions
fixed. (atomium has no `score` column, so its sequences cannot be ranked the way MPNN's can.)
**Open decision:** whether to also fix the binder-side epitope histidines (§3) or let them be
redesigned.

**Step 9 — validate.** `boltz` on three forked tables (`--table-label`): dimer alone,
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

---

## 8. Known gaps

- **No pH-responsiveness measurement anywhere** — the campaign's actual success criterion
  (§2.2). Closest future fix: PyRosetta rescoring with network histidines neutral vs.
  doubly-protonated.
- **The histidine-free epitope requirement is unsatisfiable on this pool** (§3). Accepted,
  unmeasured.
- **Binder-side histidines in 34/37** (§3). Unaddressed; the lever is at step 8.
- **The epitope footprint is up to 46 % of the binder** (§3), and nothing gates on that ratio.
- **Step 2 revalidation skipped.** Nothing independently confirms the inherited poses;
  BindCraft2's own scores come from the AF2 that designed them. If step 9 disappoints, this is
  the first place to look.
- **Phase 1 deferred**, so the pool is fixed and finite — the main yield risk (§2.6). If Phase 2
  exhausts it, Phase 1 is the unblock, and multi-target bindcraft2 would need a wrapper edit
  (`run_bindcraft2.py:438-446` hard-codes a single-element `targets` list).
- **Thresholds deliberately unset** (§2.1). They must be chosen and recorded before any claim
  about yield.
- **Mouse EGFR entirely unmeasured.** Whether His409 is conserved in mouse is a cheap sequence
  check that has not been done, and it interacts with the dual-species requirement.
- Provenance of the inherited binders is outside this repo.
- `CLAUDE.md`'s env table says `sapia-runs`; `.env` says `sapia-runs-toon`. Stale.
- rpxdock cage/dihedral architectures are out of scope; C2 is verified.
- No AF2 tool exists (relevant only if Phase 1 is revived).
