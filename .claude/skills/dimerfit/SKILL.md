---
name: dimerfit
description: How to run the custom dimerfit tool — placing a C2 dock back into the binder/target frame through one protomer, then measuring whether the partner protomer occludes the target binding site (we want it to) and whether the two protomers can be linked. Covers the premise and where its numbers are meaningless, the lineage-resolved reference and epitope columns, the chain-convention trap (binder on B, target on A), the mapping-trust columns that catch a silently wrong superposition, the two opposite-signed requirements that must not be collapsed, and every column it collects. Load before composing a dimerfit run or interpreting its columns.
---

# dimerfit

**Custom tool** (lives in `tools/dimerfit/`, not bundled with prosapia).
**`action: update`** — annotates the dock table it reads, in place. `-t` is required.
**CPU-only**; the manifest builder forces `gpus_per_task = 0`, so you do not pass `-g 0`.
**No default input column** — `-i` is effectively required (see Traps).
**Fast**: ~90 ms per design measured, so one batched task covers a whole dock table.

> **Premise.** *Place a C2 dock back into the binder–target frame through protomer A,
> and measure whether the partner protomer occludes the target binding site and
> whether the two protomers can be linked.*

Per design it:

1. reads the C2 dock — **exactly two protomers of the same binder**;
2. Kabsch-superposes dock protomer **A** onto the reference complex's **binder** chain;
3. applies that **one** transform to the **whole dimer**, so protomer B lands wherever
   the C2 operator put it relative to the target;
4. measures occlusion of the target binding site by protomer B, where the C2 interface
   sits relative to the epitope, and the C-term(A) → N-term(B) distance;
5. writes the transformed **dimer** (no target) and the dimer **plus** target.

## Why it exists

The pH-switch strategy needs a dimer that *blocks* its own binding site: the two
protomers are linked, and in the "off" state the partner protomer sits where EGFR
would. Deciding that from a dock requires putting the dock back into the frame of the
binder–target complex — which means superposing on **one chain** and moving **both**.

Nothing in the workspace could do that:

- **`pyrosetta`** (columns: `total_score`, `score_per_res`, the weighted terms,
  `sasa`, `packstat`, `buried_unsat`, `relax_ca_rmsd`, `if_dG`, `if_dSASA`,
  `if_hbonds`, `if_delta_unsat`) has **no superposition and no rigid-body machinery**
  at all — its only RMSD is a relaxed pose against its own start.
- **`ringfit`** (columns: `align_rmsd`, `seq_match_frac`, `resnum_offset`, `n_clash`,
  `n_clash_res`, `min_dist_ring`, `lipid_clash`, `min_dist_lipid`, `bsa_t1/t2/total`,
  `bridge_ratio`, `n_contact_res_t1/_t2`, `hotspot_recall`, `hotspot_hits`,
  `binder_len`, `binder_chains`) has the same *machinery* — Kabsch, trust metrics,
  write the moved structure — but a different **premise**: a binder straddling two
  adjacent protomers of a larger oligomer. Its `bridge_ratio` is undefined here, and
  it has no `link_dist`, no `occluded_frac`, no epitope-overlap columns. Same tools,
  different question → new tool, not a fork.
- **`usalign`** (`TM1`, `TM2`, `RMSD`, `ID*`, `L*`, `sup_path`) superposes two
  structures but cannot apply one chain's transform to a second chain, and emits no
  occlusion, linker or overlap geometry.
- **`cms`** (`target`, `binder`, `sc`, `sc_area`, `sc_median_dist`, `n_atoms_*`, plus
  a per-residue file) and **`ifacegeom`** (`binder_res`, `target_res`, the COMs, the
  terminus projections) both measure **one file in its own frame**. Neither can see
  the partner protomer against the target.

And it is a **tool**, not a script, because the verdicts have to be filterable
columns that ride lineage into the hbdesigner table — and because `dimerfit_path` is
the file the next step consumes.

## Scope — where these numbers are meaningless

- **Exactly two protomers.** A monomer or a C3+ assembly is an `error:` row: "the
  partner" and `link_dist` are undefined with any other count.
- **The dock must be the same protein as the reference binder.** The tool does no
  sequence alignment: it pairs residues ordinally (or by residue number) and
  *reports* whether that pairing holds (`seq_match_frac`, `resnum_offset`,
  `align_rmsd`). It cannot fix a wrong pairing, only make it visible.
- **Backbone-only docks understate the clash/occlusion counts.** Dump the dock with
  real side chains (`rpxdock --use-orig-coords`); check `rpxdock_path` is not the
  N/CA/C/O/CB+CEN form.
- **The epitope lists are taken as given.** Nothing here re-detects an interface. If
  `--epitope-column` was measured on a different structure, every epitope-derived
  number is wrong and no column will say so.
- **`n_clash`/`clash_frac` are heavy-atom distance counts, not an energy.** They are
  **not** `fa_rep`; there is no Rosetta in this image. Compare them between docks,
  never to a Rosetta number.
- **Nothing is relaxed or repacked.** Rigid-body geometry of the inputs as given.
- **Blind to non-amino-acid content.** Ligands, glycans and HETATM on either side are
  ignored everywhere (selection, clashes, occlusion). A target whose binding site is
  partly glycan will read as less occluded than it is.
- **It does not score the C2 interface.** That is `rpxdock_score` / `rpxdock_rpx`
  upstream; `dimerfit` only says *where* that interface sits.

## Invocation

```bash
sapia run dimerfit <run_dir> \
    -t table1 \
    -i rpxdock_path \
    --ref-column input_path \
    --epitope-column ifacegeom_binder_res \
    --target-epitope-column ifacegeom_target_res \
    --binder-chain-in-ref B --target-chains-in-ref A

sapia collect dimerfit <run_dir> -t table1
```

`-t` is the **dock table** (the child table rpxdock minted). `--ref-column`,
`--epitope-column` and `--target-epitope-column` are resolved **up the lineage**, so
they may live on the parent complex table — which is where `ifacegeom` wrote them.

### Flags

| Flag | Default | What it does |
| --- | --- | --- |
| `-i` / `--input-column` | *(sentinel: none)* | The C2 dock column, normally `rpxdock_path`. No honest default — see Traps. |
| `--ref-column` | `input_path` | The reference binder/target **complex**: the frame the dock is placed back into. Lineage-resolved. |
| `--epitope-column` | `ifacegeom_binder_res` | **Binder-side** epitope, `B:12,B:15,…` in the reference binder's numbering. Mapped onto the dock through the residue match. Drives `epitope_com_dist`, `n_overlap_res`, `frac_overlap`. Lineage-resolved. |
| `--target-epitope-column` | `ifacegeom_target_res` | **Target-side** binding site, `A:12,…` in the reference target's own numbering — the residues protomer B has to cover. Drives `occluded_frac`, `n_occluded_res`. Read straight off the reference, no mapping. Lineage-resolved. |
| `--binder-chain-in-ref` | `B` | The binder's chain in the reference complex. **The inherited hEGFR complexes have binder = B, target = A** — the opposite of the usual convention. |
| `--target-chains-in-ref` | `A` | The target's chains, comma-joined, or `auto` = every protein chain that is not the binder chain. |
| `--dock-chains` | `auto` | The two protomer chains, **ordered**: first = protomer A (superposed), second = protomer B (the partner). `auto` uses the dock file's own chain order. The dock must hold exactly two protein chains either way. |
| `--resnum-match` | `auto` | `ordinal` (i-th CA to i-th CA — right for rpxdock's 1..N dump), `resnum` (equal residue numbers), or `auto` (ordinal when the CA counts match, resnum otherwise). The choice is reported as `resnum_match`. |
| `--clash-cutoff` | `2.5` | Heavy-atom distance (Å) below which a protomer-B atom counts as clashing with the target. Below any real contact (C–C vdW ≈ 3.4 Å, shortest heavy-atom H-bond ≈ 2.6 Å), so a non-zero count means genuine interpenetration, not a tight interface. Same default as `ringfit`. |
| `--occlusion-cutoff` | `5.0` | Heavy-atom distance (Å) within which a **target** epitope residue counts as occluded by protomer B. First-shell contact distance, so "occluded" means "B is in contact with it". |
| `--dimer-contact-cutoff` | `5.0` | Heavy-atom distance (Å) within which a protomer-A residue counts as part of the C2 interface. Deliberately equal to `--occlusion-cutoff` so the two residue sets are comparable. |
| `--designs-per-task` | `100` | Designs per task. The work is ~90 ms per design; this exists to amortise the container cold start. |

## Traps

- **`-i` is effectively required.** `default_input_column` is the `"not applicable"`
  sentinel (as in `ifacegeom`, `cms`, `chainsel`): the builder refuses the run and
  lists the `_path` columns the table actually has.
- **Binder = chain B, target = chain A in the inherited complexes.** Backwards from
  the usual convention, and it has already caused one bug in this campaign. Get it
  wrong and the dock is superposed onto **EGFR**: `seq_match_frac` collapses and the
  row goes `warn:`, so it is caught — but only if you read that column.
- **`--dock-chains` is ordered, and the order is a *decision*.** The first chain is
  the one that keeps the binder's pose; the second is the one whose occlusion is
  measured. For a true C2 dimer of identical protomers the metrics are symmetric
  (verified: swapping gives identical numbers), but for anything that is not exactly
  C2 they are not.
- **Residue labels are `chain:resnum[icode]`** — `B:12`, never `B12`. `B12` is a named
  parse error, not a guess. (`ringfit --hotspots` uses the other convention.)
- **A chain named but absent is an `error:` for that row**, never a smaller
  selection: a silently narrower target would *understate* occlusion, i.e. fail in
  the direction that looks like a bad design rather than a bad run.
- **An epitope residue that cannot be mapped is an `error:`**, not a quietly smaller
  denominator — `occluded_frac` and `frac_overlap` are only meaningful on the exact
  residue sets.
- **A missing `--ref-column` / epitope value raises at SUBMIT time**, naming the
  design and the column. It does not silently drop the row, because a smaller
  `Submitting N designs` reads like a filter problem rather than a missing column.
- **`resnum_offset != 0` is a `warn:`, not an error.** The epitope is mapped through
  the residue *pairs*, so a renumbered reference still gives correct geometry
  (verified: a +17 shift reproduces every number). The warning exists because in this
  campaign the mapping is supposed to be the identity — a non-zero offset means
  something renumbered that nobody accounted for.
- **A blank list column is not NA-the-metric.** `dimer_iface_res`, `overlap_res` and
  `occluded_res` are empty strings when the corresponding count is genuinely 0 — read
  the `n_*` column, not the list.
- **`-l/--dir-label` for variants.** Another dock column, another reference, other
  cutoffs → the columns collide otherwise. A `-l verification` run writes
  `dimerfit_verification_*`.
- **A worker crash exits the task non-zero** (with fallback `error: worker crashed`
  TSVs); a *per-design* failure does not — the batch still exits 0 and the table
  carries the error.

## Columns collected (`dimerfit_` prefix)

`<leaf>_status` is `OK`, `warn: …`, `error: …` or `missing`.
`<leaf>_path` is the **transformed dimer, no target** — the file the next step
(`hbdesigner`) consumes.

### Mapping trust — the point of the tool, not an extra

A wrong superposition produces perfectly plausible geometry, so read these first.

| Column | Meaning |
| --- | --- |
| `align_rmsd` | CA RMSD (Å) of dock protomer A onto the reference binder. **Expect ~0** when the dock is the same scaffold rigidly re-placed. |
| `n_align_atoms` | How many CA pairs the fit used. |
| `seq_match_frac` | Fraction of those pairs whose residue identity agrees. **Expect 1.0.** Below 0.9 → `warn:` (metrics still reported). |
| `resnum_offset` | Modal `ref_resnum − dock_resnum`. **Expect 0.** Non-zero → `warn:`. |
| `resnum_match` | `ordinal` or `resnum` — how the pairs were formed. |
| `dock_chains` | What was actually used; `auto` resolved. **First = protomer A, second = protomer B.** |
| `ref_binder_chain`, `ref_target_chains` | Likewise for the reference. |
| `complex_target_chains` | The chain IDs the target got inside `complex_path` (renamed only on collision — the dimer keeps its IDs). |
| `n_protomers` | 2. Any other count is an `error:` row whose message carries the real number. |
| `n_res_a`, `n_res_b` | Amino acids per protomer. **Check against the parent's `ifacegeom_binder_len`** — an invariant this tool does not compute, so it is independent evidence. |
| `ref_binder_len` | Amino acids in the reference binder chain. |
| `n_atoms_b` | Protomer B heavy atoms = the denominator of `clash_frac`. |

### C2 interface geometry, relative to the epitope — **LOW is what we want**

| Column | Meaning |
| --- | --- |
| `dimer_iface_res`, `n_dimer_iface_res` | The C2 interface residues on the **protomer-A** side, `A:12,A:15,…` (by C2 symmetry the same residue numbers apply to B). |
| `dimer_iface_com`, `epitope_com` | Mass-weighted centres as `x,y,z`, **in the reference frame** (so they are comparable across docks of the same parent). |
| `epitope_com_dist` | `|dimer_iface_com − epitope_com|` (Å). |
| `n_epitope_res` | Epitope residues mapped onto the dock = the denominator of `frac_overlap`. |
| `overlap_res`, `n_overlap_res`, `frac_overlap` | Residues that are **both** C2 interface and epitope. **We want this LOW**: the dimer interface should sit *adjacent* to the epitope, not on top of it (a dimer interface built from the binding face cannot bind at all once dissociated). |
| `link_dist` | CA(C-term of A) → CA(N-term of B), Å. A ~20 aa linker reaches roughly 60–70 Å fully extended. **Recorded, never gated on.** |
| `cterm_to_nterm_res` | e.g. `A:111->B:1` — checkable by hand against `dimerfit_path`. |

### Occlusion of the target binding site — **HIGH is what we want**

| Column | Meaning |
| --- | --- |
| `n_clash`, `clash_frac` | Protomer-B heavy atoms within `--clash-cutoff` of the target, absolute and as a fraction of `n_atoms_b`. A geometric count, **not** an energy. |
| `min_dist_b_target` | Closest protomer-B-to-target heavy-atom distance (Å). `n_clash == 0` ⟺ this is above the cutoff — use it as the continuous version. |
| `n_target_epitope_res` | Target binding-site residues read from `--target-epitope-column` = the denominator below. |
| `occluded_res`, `n_occluded_res`, `occluded_frac` | Target epitope residues within `--occlusion-cutoff` of protomer B. **HIGH = the partner protomer blocks EGFR** — the "off" state. |
| `complex_path` | The transformed dimer **plus** target, for inspection only. |
| `seconds` | Wall time of the design's measurement. |

NA (empty) rather than 0 marks a field that did not apply: the design failed, no TSV
was written, or the geometry was undefined. **A blank `occluded_frac` is not "nothing
occluded".**

## Reading the numbers

**`frac_overlap` low and `occluded_frac` high are two different, opposite-signed
requirements, and the tool deliberately does not combine them into a score.**

- `frac_overlap` is about *which residues build the C2 interface*. If the dimer
  interface is made of the epitope itself, the protomer that dissociates has nothing
  left to bind with.
- `occluded_frac` / `n_clash` are about *where the partner's body ends up*. The
  partner must physically cover the target binding site for the dimer to be the "off"
  state.

The ideal dock has the C2 interface **beside** the epitope while the partner's bulk
**covers** it. **Measured on 40 real docks (two parents) these two are nearly
collinear — Pearson r = 0.97 between `occluded_frac` and `frac_overlap`, and 0 of the
23 docks with `occluded_frac > 0.5` had `frac_overlap ≤ 0.3`.** That is a finding
about the dock pool, not about the metrics: pick the trade-off deliberately (or widen
the pool) rather than hoping one threshold satisfies both.

`link_dist` was 10–50 Å across those 40 docks, i.e. **the linker never binds** at a
~20 aa budget. Treat it as a record, not a gate, exactly as specified.

`epitope_com_dist` (0.8–13.7 Å measured) is the continuous version of the overlap
question and does not depend on a contact cutoff — useful for ranking inside a tied
`frac_overlap` band.

### The filter the caller can now write

`-f` takes a **module** with `apply_filter(df) -> df`:

```python
def apply_filter(df):
    return df.query(
        'dimerfit_status == "OK" '
        'and dimerfit_align_rmsd < 0.5 '
        'and dimerfit_seq_match_frac > 0.99 '
        'and dimerfit_resnum_offset == 0 '
        'and dimerfit_occluded_frac >= 0.6 '
        'and dimerfit_frac_overlap <= 0.35 '
        'and dimerfit_link_dist <= 60'
    )
```

The first four clauses are the trust gate; drop them only if you enjoy ranking
nonsense. The last three are the science, and the thresholds are the campaign's call.

## Verification done

All of the below was run **locally** (worktree `.venv`, gemmi 0.7.5 / numpy 2.5.3) on
files pulled from the `sapia-runs-toon` Volume — real docks from
`outputs/20261002_143419_dimer_phase2`, real BindCraft2 references from
`inputs/bc2_output_for_anthony/`. The tool has **not** yet been submitted through
`sapia run` on Modal; that first run is the orchestrator's `-l verification`.

**1. Synthetic case with a known answer.** Built a C2 dimer as `[G·binder, G·S·binder]`
from the reference binder, with `S` a true 180° operator and `G` an arbitrary
"junk" transform (37° rotation + a 300 Å translation), all computed with plain numpy
outside the tool:

- `align_rmsd` = **0.000**, `seq_match_frac` = 1.000, `resnum_offset` = 0.
- `link_dist` = **34.262 Å** = the independently computed CA(A:111)→CA(B:1) distance,
  to the digit; `cterm_to_nterm_res` = `A:111->B:1`.
- Protomer B in the output lands within **0.0011 Å** (PDB rounding) of the
  independently computed `S·binder` — i.e. the **single** transform really was applied
  to the whole dimer.

**2. Real data, 40 docks** (parents `v1_l111_c947eb80`, `v1_l90_30b6751c`), 3.7 s total:

- 40/40 `OK`, empty stderr.
- `align_rmsd` 0.001 Å, `seq_match_frac` 1.000, `resnum_offset` 0 on **every** row —
  the identity map holds, as the campaign assumed.
- **`n_res_a` == `n_res_b` == the parent's `ifacegeom_binder_len` on 40/40** — the
  invariant the tool does not compute.
- **`epitope_com` reproduces the parent's `ifacegeom_binder_iface_com` to within
  0.05 Å on 40/40.** ifacegeom computed that COM on the reference complex, dimerfit
  computes it on the dock after the superposition: they agree, so protomer A really
  does land where the reference binder sat.
- Independent re-checks on the written files (gemmi's own `superpose_positions`, not
  the tool's Kabsch): protomer A vs the reference binder = **0.0007 Å** direct CA RMSD
  with no refitting; all pairwise distances inside the dimer preserved to
  **1.4e-3 Å** (rigid); A→B is a **180.00°** rotation (the C2 survived the move);
  `link_dist` identical before and after the transform; the target inside
  `complex_path` is bit-identical (0.0005 Å) to the untransformed reference.
- Distributions: `occluded_frac` 0–1, `frac_overlap` 0–0.97, `n_clash` 0–710,
  `link_dist` 10.4–50.3 Å, `n_dimer_iface_res` 15–51.

**3. Deliberately bad inputs** → `error:` status, **every metric NA**, batch still
exits 0 (verified with an `OK` control row in the same batch):

| Input | Status |
| --- | --- |
| monomer instead of a dimer | `error: … holds 1 protein chain(s) (A), expected exactly 2 protomers of a C2 dock: n_protomers=1` |
| missing dock file | `error: structure missing: nope.pdb` |
| `--binder-chain-in-ref Z` (absent) | `error: --binder-chain-in-ref 'Z' absent from ref.cif (protein chains present: A,B) …` |
| epitope written `B12,B15` | `error: cannot parse --epitope-column residue 'B12': expected 'chain:resnum' …` |
| epitope on the wrong chain (`A:12`) | `error: --epitope-column lists chain(s) ['A'] but --binder-chain-in-ref is 'B' …` |
| epitope residue `B:999` | `error: epitope residue(s) B:999 are not among the 111 matched dock/reference residues …` |
| target epitope residue `A:9999` | `error: --target-epitope-column residue(s) A:9999 are not in the reference target chain(s) A …` |

**4. Warn paths.** Reference binder renumbered **+17** → `warn: resnum_offset 17`
with every metric *identical* to the unshifted run (the epitope goes through the
pairs, not the numbers). A dock scored against the **wrong parent's** reference →
`warn: seq_match_frac 0.06`, `align_rmsd` 11.8 Å, plus a stderr warning: the row
stays out of any `== "OK"` selection but keeps its numbers for diagnosis.

**5. Submit-time guards** (fake `ManifestCtx`): sentinel/unknown `-i` refused with the
available `_path` columns listed; a missing lineage column or reference file raises
naming the design; `--dock-chains A,B,C` / `A,A`, a binder chain also listed as a
target chain, a non-positive cutoff and `--designs-per-task 0` all refused;
`gpus_per_task` forced to 0; batching correct (5 designs at 2/task → 3 tasks of
2/2/1), 9 manifest fields, 5 sub-manifest fields.

**6. Collector**: 34 columns; `error:` rows come back with 0/34 non-NA; `warn:` rows
keep all of them; a design with no TSV collects as `missing`.

**7. The task script** (`dimerfit.sh`, run under a stub prelude): a 1-field manifest
line, an invisible task file, a line number beyond the manifest (the measured stale-
Volume case → 5 retries, then FATAL), and a **task file holding fewer designs than
the manifest says** all fail loudly with rc=1 and a named message; the fallback
`error: worker crashed` TSV is written for the affected designs. The good path runs
end to end through the script.

### Not verified

- **Never run through `sapia run` / `sapia collect` on Modal or SLURM.** The Modal
  image (`debian_slim` + gemmi + numpy, same as ringfit/ifacegeom) has not been built.
- `--resnum-match resnum` on a reference with genuine numbering **gaps** (only the
  dense case was exercised).
- Docks with ligands or HETATM content, and multi-chain targets (`--target-chains-in-ref`
  with more than one chain) — the code path exists, no real input was available.
- Whether any threshold on `occluded_frac` / `frac_overlap` predicts a real pH switch.
  That is a wet-lab question; this tool only supplies the numbers.
