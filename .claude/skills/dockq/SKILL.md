---
name: dockq
description: >-
  DockQ v2 on predicted binders: did the prediction keep the reference binding mode over every native binder-target interface? Writes Fnat, iRMSD, LRMSD, DockQ and CAPRI class. Load before composing a dockq run or interpreting its numbers. Covers the explicit binder mapping, the two residue-correspondence modes and the trust columns.
---

# dockq

**Custom tool** (`tools/dockq/`, not bundled with prosapia). **`action: update`**: it annotates the table it ran against and mints no child. **`-t` is required**, and so is `-i/--input-column` (no honest default model column). CPU only: the manifest builder forces `gpus_per_task = 0`, so you never need `-g 0`. A `modal_image.py` ships (slim Debian + `build-essential` + `python3-dev` + `DockQ==2.1.3` and pinned deps). **There is no `SAPIA_ACTIVATE_DOCKQ` anywhere, so on vib every task would die at activation**; it is Modal-only until someone writes one.

It answers one question: **per native binder-target interface, how well does the model reproduce the reference's binding mode, judged by DockQ v2?** It wraps the package and only selects, maps, reports and cross-checks: no part of Fnat, iRMSD, LRMSD or DockQ is computed by this tool.

**No column here has come back from a run on Modal yet** — treat every number below as unverified until one has, and run that one under a `verification` label (see `verifying-a-tool`). Unknown until then: whether the image builds (`build-essential` plus the sdist compile), whether `structure_path` resolves per row from the real table, and what the parents' `nat` counts are as that table sees them.

## Premise, and when its numbers are meaningless

Per design, two structures with the same layout:

- a **model** (`-i`, e.g. `boltz_cof_path`): one binder chain (`--binder-chain`) plus N target chains (`--target-chains`);
- a **reference** (`--ref-col`, default `structure_path`, resolved up the lineage, so on a child table it reads the root table's parent complex): one binder chain (`--ref-binder-chain`) plus the same N target chains (`--ref-target-chains`).

The target may be N identical copies of one protomer (the TNF-alpha C3 trimer: chains A, B, C, 152 res each). The model binder must be a **substitution variant** of the reference binder: same length, same residue numbering, at least `--min-binder-identity` of positions identical.

**In this workspace the reference is itself a prediction** (the campaign's input complex, a Boltz/AF model), not an experimental structure. A DockQ of 0.95 means "the same pose as the parent", never "the correct pose". A native-sequence control that comes back incorrect means the prediction recipe cannot reproduce the parent pose. It does not mean DockQ is wrong, and it does not mean the parent pose is right.

| Situation | What happens |
| --- | --- |
| model binder unrelated to the reference binder (another parent's binder, wrong chain letter) | **ERROR.** In `number` mode a different length errors on numbering; at equal length, identity below `--min-binder-identity` (0.5) errors on identity. |
| binder with insertions/deletions vs the reference | **ERROR** in `number` mode (numbering differs). In `sequence` mode it runs, and DockQ compares only the identical aligned positions. |
| target chains that are not the reference target | **ERROR**: DockQ's `group_chains` finds no match at `--allowed-mismatches` (default 0), or every permutation exceeds the allowance. |
| target permutations that disagree on which native binder interfaces exist | **ERROR** naming the differing key sets: the homo-oligomer symmetry assumption does not hold. |
| a named chain absent from either file | **ERROR**, naming the chains that are present. |
| the reference binder touches none of the target chains (no residue pair within 5 A) | **ERROR**: nothing to reproduce. `n_native_binder_interfaces = 0` is kept on the row as a trust fact. |
| residue numbers that DockQ's `no_align` encoding cannot represent (`number` mode) | DockQ encodes each residue number as `chr(resn + offset)` with offset >= 45, so a number above about **1,114,000** errors the row. A number encoding to `-` (45) would read as a gap — harmless only because of the identical-numbering guard below. |
| more than one binder, an antibody's two chains, a heteromeric target with repeated chains of different roles | **Out of scope.** The mapping search assumes one homo-oligomeric target plus one binder. A heteromer would trip the identical-interface-set check or map each chain only onto its sequence twin; never checked. |
| the best interface is an INCIDENTAL one (1-2 native contacts) | DockQ counts an interface as native with a single residue pair within 5 A, and there Fnat can only be 0, 0.5 or 1, so a wrong pose can reach 0.23 on it. `--min-native-contacts N` excludes such interfaces from `dockq_best`, `dockq_mean_natw`, `n_interfaces_acceptable` and `dockq_best_full_objective` (the threshold gates DockQ's own-mapping best too, not only the kept mapping's; excluded interfaces stay in `per_interface` as `used=0`); **default 1 = DockQ's own rule, i.e. no exclusion**. Read `nat_best` beside `dockq_best`. A row with no interface at or above N **ERRORS**. |
| model and reference are the same file, or `-i` equals `--ref-col` | The row errors (same file); the submit is refused (same column). |

**What DockQ measures, as defined by its authors** (thresholds are DockQ's, not ours): a native contact is a binder-target residue pair with any heavy atoms within **5 A** in the reference; the interface for iRMSD is every residue within **10 A** of the partner; RMSDs use backbone **N, CA, C, O**. **LRMSD** superposes the model onto the reference on the **receptor**, the larger chain of the pair (here one target protomer), and then takes the binder backbone RMSD without refitting — exactly the target-fixed binder RMSD. It is per protomer, so a binder spanning two protomers gets two LRMSDs that differ slightly with the protomer's own prediction error.

## Trust columns: read these before the science

| Column | Expected | How to read it |
| --- | --- | --- |
| `binder_identity` | **1.0** for a native-sequence control; **1 - n_mut/len** for a design | Fraction of binder positions identical to the reference. **How many substitutions a design carries — and so what this column reads — is a property of YOUR batch, not of the tool**: it varies between design arms of one campaign (a capped arm and an uncapped one differ severalfold in mutation count). Take `n_mut` from the designer's own column for the rows you are scoring, or read the `dockq_binder_identity` distribution off your own run; never inherit a number from another campaign, here or elsewhere. Near 1 - n_mut/len means the right binder was compared; far below it means the wrong reference or the wrong chain. |
| `binder_n_res` / `ref_binder_n_res` | equal | Unequal means the wrong reference. In `number` mode it is already an error. |
| `binder_n_res_compared` | `binder_n_res` in `number` mode; **the identical positions only** in `sequence` mode | How many binder residues DockQ actually used. |
| `target_max_mismatches` | **0** | Positions differing in the worst target pair **of the kept mapping**. A permutation that exceeds `--allowed-mismatches` is never a candidate. |
| `target_max_mismatches_any_mapping` | **0** | The same over every permutation DockQ produced, including ones skipped for exceeding the allowance. Non-zero with the kept value 0 means some permutations were excluded; `per_mapping` shows which. |
| `mapping` | a cyclic permutation for a C3 trimer, e.g. `B>A;C>B;A>C;D>D` | **model>reference**, in reference order. The letter in `per_interface` and `best_ref_chain` is the REFERENCE chain; the model chain is `model=` in each entry. Boltz assigns A/B/C to the copies arbitrarily, so a non-identity mapping is normal. |
| `mapping_full_objective` / `mapping_agrees` | usually True | What DockQ's OWN objective (the sum over ALL native interfaces, target-target included) would have picked, and whether it is the same mapping. **False is a finding, not an error:** the target-target interfaces outvoted the binder. Read `dockq_best_full_objective` beside it. |
| `n_mappings_tried` | **6** for a trimer (3!) | The target permutations DockQ's `group_chains`/`get_all_chain_maps` produced. 1 means the targets did not cluster as identical copies. |
| `per_mapping` | one entry per permutation, every one `n_binder=` the same | `<map, ','-joined>:target_mm=:n_binder=:binder_sum=:full_sum=`, or `:skipped`, `;`-joined. **Empty** = the row failed before the first permutation was checked; **shorter** on an error row = the permutations computed before the failure. On a C3 trimer the three cyclic permutations score a full sum near the maximum (2.8-2.9 measured) and the three reflections far lower (~0.2): the symmetry, visible. |
| `min_native_contacts` / `n_interfaces_used` | 1 / = `n_native_binder_interfaces` at the default | The `--min-native-contacts` echo, and how many native binder interfaces reach it. Fewer used than native = some interfaces were excluded as incidental. |
| `n_native_binder_interfaces` | 1-3 | Target chains the reference binder contacts. The number of `per_interface` entries, and the number of terms in every candidate's binder sum. |
| `model_chain_ids` / `ref_chain_ids` | e.g. `A,B,C,D` | Every protein chain in each file (auth_asym_id, HETATM excluded). |
| `ref_path` | the parent complex | The reference this row was scored against. Count the distinct values: one per parent. |
| `dockq_version` | `2.1.3` | `unknown` if the installed metadata could not be read; the numbers are unaffected and the row stays OK. |
| `residue_match`, `allowed_mismatches`, `min_binder_identity`, the four chain echoes | as passed | Rows are comparable only at equal values of these. |

### The chain mapping, and why it is not DockQ's default

DockQ v2 maps model to reference chains **by sequence**: two chains may correspond only if their alignment has at most `allowed_mismatches` mismatches (default **0**). At that default **a single substitution is enough to match nothing**, whatever the design arm's mutation count happens to be, and DockQ exits with `For chains ['D'] no identical corresponding chain was found`. Raising `allowed_mismatches` far enough to admit the substitutions your own batch carries (~18 for a 30 % variant of a 60-mer) would work, but it applies to every chain pair, so it is also an invitation for the binder to cluster with a target chain — and the number needed is a property of the batch, to be read off the designer's mutation counts, not a constant of this tool. So the tool does what DockQ's `--mapping D*:D*` does:

1. **The binder is mapped explicitly** (reference binder <- model binder). Its relatedness is then checked by `--min-binder-identity` and, in `number` mode, by identical numbering.
2. **Only the target chains go through DockQ's own `group_chains`**, at `--allowed-mismatches`. Every permutation `get_all_chain_maps` yields is checked: one whose target pairs differ at more than the allowance is skipped, each remaining one is scored with `run_on_all_native_interfaces`.
3. **Every scored permutation must yield the SAME set of native binder interfaces**, or the row errors. This is the homo-oligomer symmetry assumption made checkable, and it guarantees that the binder sums compared below have the same number of terms, so no permutation can win on count.
4. **The permutation kept maximises the summed DockQ of the BINDER interfaces only**, where DockQ's CLI maximises the sum over ALL native interfaces, including the target-target ones. In a symmetric trimer those are nearly equal across the cyclic permutations, so they mostly do not move the choice; when the binder interfaces are all poor, small target-target differences can decide it. DockQ has no option to restrict its objective to chosen interfaces (`--mapping :ABD` restricts chains, not interfaces), which is why the tool does the selection itself. Ties go to the higher full sum, then to DockQ's enumeration order.

**The kept mapping maximises the binder-interface SUM, not the best single interface — so the summed binder DockQ can only flatter the design, while `dockq_best` can come out BELOW DockQ's own choice.** Example: kept 0.5/0.5 (sum 1.0, best 0.5) vs DockQ's 0.7/0.2 (sum 0.9, best 0.7). `dockq_best_full_objective` is the `dockq_best` DockQ's own choice would have given, so it can sit **above or below** `dockq_best`, and the difference is on every row. **Near the 0.23 gate, a row with `mapping_agrees = False` must be read with `dockq_best_full_objective`.** If the two straddle the bar, the pass depends on the mapping rule, not on a clear pose.

### Residue correspondence, and why the default differs from DockQ's

Inside a mapped chain pair, **DockQ's default aligns by sequence and then keeps only IDENTICAL positions** (`get_aligned_residues` keeps `|` columns only). For a design, every substituted residue is silently dropped from iRMSD and LRMSD and can never count as a reproduced contact, while Fnat's denominator is all native contacts of the **full** reference (`calc_DockQ` computes `nat_total` on the untouched native). So **a design cannot reach Fnat 1.0 even in a perfect pose**, and the cap grows with the number of mutated interface residues — which for a designed-His binder sit exactly at the interface. Comparing designs with different mutation counts would confound pose with mutation count.

**`--residue-match number` (the default here, not DockQ's) passes DockQ's own `no_align=True`:** residues correspond **by residue number** and every position is compared. iRMSD and LRMSD are backbone only (N, CA, C, O), so a substitution does not change what they compare. **Fnat does NOT escape the mutations even here:** a contact is any HEAVY-ATOM pair within 5 A (`get_residue_distances(all_atom=True)`), side chains included, so a design that shrinks an interface side chain can lose a native contact with a perfect backbone, and one that grows a side chain can gain contacts. In `number` mode Fnat therefore still mixes pose with mutation, though less than in `sequence` mode, where the substituted residue cannot count at all. **The identical-numbering guard is what makes `no_align` exact:** the tool requires the residue-number list of every mapped chain pair to be identical, or the row errors. With identical lists, DockQ's encoded "sequences" are identical and `get_aligned_residues` returns every residue in order; without the guard, DockQ would align the encoded numbers with gaps and silently drop residues. Boltz numbers every chain from 1. If a reference keeps its original numbering (e.g. 77-233) while the prediction is renumbered from 1, this mode errors rather than guessing; use `sequence` mode or renumber the reference.

## The question it answers, and the bar

**Did the predicted binder keep the parent's binding mode?** DockQ's classes are the CAPRI docking conventions: **< 0.23 incorrect, >= 0.23 acceptable, >= 0.49 medium, >= 0.80 high** (DockQ v2 README). The campaign bar is **acceptable or better on the best interface**: `dockq_dockq_best >= 0.23`. These are conventions for docking-pose quality, not for binding: an acceptable pose can still have a different contact set, and a high one says nothing about affinity.

"Best" is the native binder interface with the highest DockQ, among those with at least `--min-native-contacts` native contacts (all of them at the default 1). `fnat_best`, `irmsd_best`, `lrmsd_best` and `best_ref_chain` all describe **that same interface**, never the best of each metric. An exact tie goes to the alphabetically first reference chain; **the `best_ref_chain` letter is not special**, since Boltz's chain lettering of identical copies is arbitrary. **Read `nat_best` beside it:** a best interface with 1-2 native contacts has a Fnat that can only be 0, 0.5 or 1, so its DockQ can clear 0.23 for a wrong pose. **Choosing `--min-native-contacts` is the user's decision, to be made from the per-interface `nat` counts of both parents** (measured locally on the parent files: tnf A=13, C=31; f3 A=28, B=32 — no incidental interface in either, which makes the default harmless there but is not a reason to keep it for another target). For a binder that spans two protomers, also read `dockq_mean_natw` or `n_interfaces_acceptable`: one protomer reproduced and the other lost means the binder pivoted on one contact patch.

## Invocation and flags

```bash
sapia run dockq <run_dir> -t table1_input_readme_tnf -i boltz_cof_path \
    --binder-chain D --target-chains A,B,C \
    --ref-binder-chain D --ref-target-chains A,B,C
sapia collect dockq <run_dir> -t table1_input_readme_tnf
```

(On Modal both run inside `sapia modal-shell --cmd '...'` from `/runs`; detect a refused submit from `<out_dir>/dockq_logs/dockq_modal.json`, never from the exit code. Modal only.)

| Flag | Default | Meaning |
| --- | --- | --- |
| `-i/--input-column` | **none: required** (sentinel `not applicable` raises, listing the `_path` columns) | The MODEL complex column. Resolved up the lineage. |
| `--ref-col` | `structure_path` | The REFERENCE complex column, resolved up the lineage. A column this table HAS is read from the row's own cell; it never falls back to an ancestor (the hisbond rule). Must differ from `-i` (refused at submit). |
| `--binder-chain` | **none: required** | Model binder, exactly one chain. |
| `--target-chains` | **none: required** | Model target chains, mini-language `A,B,C` or `A:C`. Must not contain the binder; no duplicates. |
| `--ref-binder-chain` | **none: required** | Reference binder, exactly one chain. |
| `--ref-target-chains` | **none: required** | Reference target chains, same count as `--target-chains`. Their order does NOT fix the mapping; it is searched. |
| `--allowed-mismatches` | `0` | DockQ's `allowed_mismatches`, applied to the **target** chains only. |
| `--residue-match` | `number` | `number` = DockQ `no_align` (by residue number, all positions, numbering must match). `sequence` = DockQ's own default (identical aligned positions only; see above). |
| `--min-binder-identity` | `0.5` | Identity floor for the explicitly mapped binder. 30 % substitutions = 0.7; the other TNF parent's binder measured 0.27-0.36. |
| `--min-native-contacts` | `1` (= DockQ's own rule, no exclusion) | Interfaces with fewer native contacts are excluded from `dockq_best`, `dockq_mean_natw`, `n_interfaces_acceptable` and `dockq_best_full_objective`, and kept in `per_interface` as `used=0`. The set excluded depends only on the reference (native contact counts), not on the mapping; the mapping objective still sums every native binder interface. No interface at or above N = error row. |
| `--designs-per-task` | `20` | About 0.4 s and ~150 MB per design measured locally (1 CPU). |

Resources: `RESOURCES = {"cpu": 1, "memory": "2G", "timeout": "00:20:00"}` (Modal) and the matching `#SBATCH` lines. Submit-time refusals (nothing queued, but a phantom registry row, as for every tool): binder in its own target set (model or reference), more than one binder chain, a duplicate chain, unequal target counts, `-i` equal to `--ref-col`, no `-i`, `--min-binder-identity` outside [0, 1], negative `--allowed-mismatches`, `--min-native-contacts` below 1, and every row resolving no model or no reference.

The builder prints `dockq: N design(s) against M distinct reference structure(s) from 'structure_path'`. **M should equal the number of parent complexes.** M = 1 on a two-parent table means the reference column did not resolve per lineage.

## Gotchas, each with its signature

- **DockQ's default refuses a design.** Running the CLI on a design without `--mapping` prints `Model chain D, native chain D: 8 mismatches ... no identical corresponding chain was found` and exits. That is DockQ working as designed, not a broken file. The tool avoids it by mapping the binder explicitly.
- **DockQ calls `sys.exit()` on several failures.** The worker converts `SystemExit` into an error row (`error: DockQ group_chains (...) exited`). Without that, one bad design would kill the whole task.
- **DockQ caches results on Biopython chain objects.** `align_chains`, `get_residue_distances`, `get_aligned_residues`, `subset_atoms`, `run_on_chains` and `list_atoms_per_residue` are `lru_cache`d, so the worker clears every cache in the module before each design and names what it cleared in its log line (`caches cleared between designs: ...`). *Signature:* a function missing from that list means the package changed.
- **`sequence` mode under-reads designs.** *Signature:* `binder_n_res_compared < binder_n_res` and a `fnat_best` that can never reach 1.0 even when `lrmsd_best` is small.
- **`mapping_agrees = False`.** DockQ's own objective would have scored a different permutation. The reported numbers are for the binder-optimal mapping; `dockq_best_full_objective` gives DockQ's, and `per_mapping` gives every permutation's sums.
- **A blank OWN model cell gives an EMPTY `dockq_status`, not `error:`.** prosapia's `filter_ready` drops the row at collect; the task did write `error: could not resolve ...`. Count non-empty statuses against the row count.
- **`-f` acts on the run only.** Rows it leaves out are never submitted, but collect still iterates every ready row and marks the unsubmitted ones `missing` under that label. *Signature:* `missing` rows exactly where the filter excluded them. Run a table that is already the subset you want, whole.
- **The package is an sdist only.** `DockQ==2.1.3` has no wheel on PyPI; pip compiles its Cython extension, so the image needs `build-essential` (`python3-dev` rides along as insurance). If the extension is missing, DockQ's fallback import fails inside the package and every row reads `error: could not import the DockQ package: ...`. It never returns different numbers.
- **Crashed worker.** Every design of that task gets `error: worker crashed` from `dockq.sh`'s fallback.
- **An OK row that fails its own arithmetic is an error row.** The collector re-parses `per_interface` (written at full repr precision, so no rounding can move a value across 0.23 or change a tie) and requires: each entry's `used` = (`nat` >= `min_native_contacts`); over the used entries, `dockq_best` = the largest DockQ, `best_ref_chain` and `nat_best` = that entry's, `dockq_mean_natw` = sum(dockq_i * nat_i)/sum(nat_i) (both to 1e-6), `n_interfaces_acceptable` = the count >= 0.23; entry count = `n_native_binder_interfaces`; used count = `n_interfaces_used`. *Signatures:* `error: dockq_best ... != max entry ...`, `error: per_interface entry '...' ...`, `error: dockq_mean_natw ... != ...`, `error: best_ref_chain ... != entry ...`. An OK row must also carry every science column plus `mapping`, `mapping_full_objective`, `mapping_agrees`, the binder identity/length facts, `target_max_mismatches`, `n_native_binder_interfaces`, `n_interfaces_used`, `min_native_contacts` and `ref_path`, or it reads `error: OK row missing ...`. Descriptive metadata (`dockq_version`, chain-id lists) is not in that gate.

## NA is never 0

An `error: ...` or `missing` row has every science column (`dockq_best`, `fnat_best`, `irmsd_best`, `lrmsd_best`, `class_best`, `best_ref_chain`, `nat_best`, `dockq_mean_natw`, `n_interfaces_acceptable`, `dockq_best_full_objective`, `per_interface`) **NA**. It keeps the echo of the flags and whatever trust facts were established before the failure (e.g. `binder_n_res` 44 vs `ref_binder_n_res` 60 on a cross-parent row; `per_mapping` when the interface sets disagreed).

**A native interface the model does not reproduce is NOT NA.** DockQ scores it (Fnat 0.0, iRMSD and LRMSD of 10-30 A, DockQ ~0.03) and that measured value is reported — DockQ's score is never exactly 0 even for a binder on the far side of the target, because its RMSD terms are 1/(1+(rms/d0)^2). So "not reproduced" reads as `fnat=0.0` with a DockQ below ~0.05, not as a zero DockQ. An interface that is not native (the reference binder does not touch that protomer) is not scored at all and does not appear in `per_interface`; a model binder that newly contacts it goes unpenalised except through its LRMSD.

## What it collects (`dockq_` prefix, or `dockq_<label>_` with `-l`)

Every column in the trust table above is collected under the same prefix. **Read those first**; the science columns below mean nothing without them.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `dockq_best` | highest DockQ over native binder interfaces | the gate; >= 0.23 acceptable |
| `fnat_best`, `irmsd_best`, `lrmsd_best` | Fnat, iRMSD (A), LRMSD (A) **of that same interface** | LRMSD = binder backbone RMSD with that protomer superposed |
| `class_best` | `incorrect` / `acceptable` / `medium` / `high` | CAPRI class of `dockq_best` |
| `best_ref_chain` | reference target chain of the best interface | a REFERENCE letter; ties alphabetical; the letter is not special |
| `nat_best` | native contacts of the best interface | 1-2 = incidental: Fnat can only be 0, 0.5 or 1 there, so do not trust a pass carried by it |
| `dockq_mean_natw` | DockQ mean over the USED native binder interfaces, **weighted by native contact count** (DockQ's `nat_total`) | weighted so that a 2-contact incidental interface (Fnat quantised in steps of 0.5) cannot count as much as the 30-contact real one; the plain mean is recomputable from `per_interface` |
| `n_interfaces_acceptable` | used native binder interfaces with DockQ >= 0.23 | equal to `n_interfaces_used` = the whole (non-incidental) footprint kept |
| `dockq_best_full_objective` | `dockq_best` under DockQ's OWN mapping choice | equals `dockq_best` when `mapping_agrees`; otherwise it can be above OR below it (the kept mapping maximises the binder sum, not the max) |
| `per_interface` | `;`-joined, reference-chain order, full precision: `A:dockq=0.9479342446364449:fnat=1.0:irmsd=0.6311189037543716:lrmsd=0.6489641488883116:class=high:nat=13:model=B:used=1` | `nat` = native contacts, `model` = model chain mapped onto that reference chain, `used=0` = excluded by `--min-native-contacts` |
| `_status` | `OK` or `error: ...` or `missing` | `OK` is the only proof of success, and an OK row has passed the arithmetic cross-checks |
| `_path` | the MODEL structure scored | |

On disk: `<run_dir>/<table>/dockq[_<label>]/<name>.tsv` per design (floats written with full precision), and `dockq_tasks/task_<n>.tsv` (name, model, reference) per task.

## A filter a caller can write

```python
# gate_dockq.py  -- for `-f gate_dockq.py` on the next tool, against the annotated table
LABEL = ""  # the -l used for the dockq run, if any
MIN_NAT = 1  # decide from the parents' per-interface nat counts; 1 = no extra condition
P = f"dockq_{LABEL}_" if LABEL else "dockq_"


def apply_filter(df):
    ok = (
        (df[P + "status"] == "OK")
        & (df[P + "target_max_mismatches"] == 0)
        & (df[P + "binder_identity"] >= 0.6)
        & (df[P + "dockq_best"] >= 0.23)
        # A pass carried by an incidental interface is not a pass. Set MIN_NAT once
        # the per-interface nat counts are known (or rerun with --min-native-contacts,
        # which does the same inside the tool).
        & (df[P + "nat_best"] >= MIN_NAT)
    )
    # A pass that exists only under the binder-optimal mapping is not a clear pass:
    # require DockQ's own mapping to agree, or its dockq_best to clear the bar too.
    robust = (df[P + "mapping_agrees"] == True) | (  # noqa: E712 - column is object dtype
        df[P + "dockq_best_full_objective"] >= 0.23
    )
    return df[ok & robust]
```

`mapping_agrees` comes back as an object column (True/False/NaN), so compare with `== True` rather than relying on truthiness. With `-l v1` every column is renamed `dockq_v1_<name>`, which is why the prefix is interpolated. A forked table (one minted with `--table-label` upstream) changes the **table name** you pass to `-t`, not the column names. For a design-set comparison, read `dockq_dockq_best` and `dockq_lrmsd_best` per set (a distribution question for `table-analyst`, not a column), always beside the native-sequence control of the same parent and recipe.

## What it is blind to

- **Side chains and chemistry.** Backbone RMSDs plus residue-level contacts. Whether the designed His reaches the target is `hisbond`; interface energy and packing are `pyrosetta` (`if_dG`) and `cms` (CMS, SC).
- **The reference's own accuracy.** The parent is a prediction. "High" means the same pose as the parent.
- **Confidence.** A pose can reproduce the parent with low pLDDT/ipSAE, or miss it confidently. Read `boltzconf` beside it.
- **A different pose on the same footprint.** A binder rotated or flipped in place that still touches the same target residues earns partial Fnat. LRMSD exposes it: measured on one design, the centroid moved only ~4 A but the binder CA RMSD is 12.6 A. Read `lrmsd_best` beside `fnat_best`.
- **New, non-native contacts** to a protomer the reference binder never touched: not scored as an interface.
- **Binding.** CAPRI classes grade docking poses; they say nothing about affinity or specificity.
- **Whole-complex agreement.** For a TM-score/RMSD over the whole complex use `usalign`. Its numbers are a different quantity (e.g. whole-complex superposition), so do not expect them to equal LRMSD.
