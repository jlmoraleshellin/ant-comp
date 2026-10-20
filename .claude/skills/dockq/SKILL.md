---
name: dockq
description: >-
  Did the predicted binder keep the reference binding mode over every native binder-target interface? DockQ v2 on predicted binders: writes Fnat, iRMSD, LRMSD, DockQ and CAPRI class per native interface, with an explicit binder mapping and two residue-correspondence modes. Load before composing a dockq run or interpreting its numbers. Traps: the reference is itself a prediction, `--ref-col` must resolve to the parent's model, and a pass can ride on an incidental 1-2 contact interface.
---

# dockq

**Custom tool** (`tools/dockq/`, not bundled with prosapia). **`action: update`**: annotates the table it ran against and mints no child. **`-t` is required**, and so is `-i/--input-column` (no honest default model column). CPU only (the manifest builder forces `gpus_per_task = 0`, no `-g 0`). On vib `SAPIA_ACTIVATE_DOCKQ` points at the lab's shared `dockq.sh` (DockQ pins numpy < 2, so it is not served by `prosapia_python.sh`); on Modal a `modal_image.py` ships (slim Debian + `build-essential` + `python3-dev` + `DockQ==2.1.3` and pinned deps).

It answers one question: **per native binder-target interface, how well does the model reproduce the reference's binding mode, judged by DockQ v2?** It wraps the package and only selects, maps, reports and cross-checks: no part of Fnat, iRMSD, LRMSD or DockQ is computed by this tool.

## Premise, and when its numbers are meaningless

Per design, two structures with the same layout:

- a **model** (`-i`, e.g. `boltz_path`): one binder chain (`--binder-chain`) plus N target chains (`--target-chains`);
- a **reference** (`--ref-col`, default `structure_path`, resolved up the lineage, so on a child table it reads the root table's parent complex): one binder chain (`--ref-binder-chain`) plus the same N target chains (`--ref-target-chains`).

The target may be N identical copies of one protomer (e.g. a C3 trimer, chains A, B, C). The model binder must be a **substitution variant** of the reference binder: same length, same residue numbering, at least `--min-binder-identity` of positions identical.

**Here the reference is itself a prediction**, not an experimental structure: a DockQ of 0.95 means "the same pose as the parent", never "the correct pose". A native-sequence control that comes back incorrect means the prediction recipe cannot reproduce the parent pose; it does not mean DockQ is wrong or the parent pose right.

| Situation | What happens |
| --- | --- |
| model binder unrelated to the reference binder (another parent's binder, wrong chain letter) | **ERROR.** In `number` mode a different length errors on numbering; at equal length, identity below `--min-binder-identity` (0.5) errors on identity |
| binder with insertions/deletions vs the reference | **ERROR** in `number` mode (numbering differs). In `sequence` mode it runs and DockQ compares only the identical aligned positions |
| target chains that are not the reference target | **ERROR**: DockQ's `group_chains` finds no match at `--allowed-mismatches` (default 0), or every permutation exceeds the allowance |
| target permutations that disagree on which native binder interfaces exist | **ERROR** naming the differing key sets: the homo-oligomer symmetry assumption does not hold |
| a named chain absent from either file | **ERROR**, naming the chains present |
| the reference binder touches none of the target chains (no residue pair within 5 Å) | **ERROR**: nothing to reproduce; `n_native_binder_interfaces = 0` is kept as a trust fact |
| residue numbers DockQ's `no_align` encoding cannot represent (`number` mode) | DockQ encodes each residue number as `chr(resn + offset)` with offset ≥ 45, so a number above about **1,114,000** errors the row. A number encoding to `-` (45) would read as a gap, harmless only because of the identical-numbering guard |
| more than one binder, an antibody's two chains, a heteromeric target with repeated chains of different roles | **Out of scope.** The mapping search assumes one homo-oligomeric target plus one binder; a heteromer would trip the identical-interface-set check or map each chain only onto its sequence twin (never checked) |
| the best interface is an INCIDENTAL one (1–2 native contacts) | DockQ counts an interface as native with a single residue pair within 5 Å, where Fnat can only be 0, 0.5 or 1, so a wrong pose can reach 0.23 on it. `--min-native-contacts N` excludes such interfaces from `dockq_best`, `dockq_mean_natw`, `n_interfaces_acceptable` and `dockq_best_full_objective` (the threshold also gates DockQ's own-mapping best; excluded interfaces stay in `per_interface` as `used=0`); **default 1 = DockQ's own rule, no exclusion**. Read `nat_best` beside `dockq_best`. A row with no interface at or above N **ERRORS** |
| model and reference are the same file, or `-i` equals `--ref-col` | The row errors (same file); the submit is refused (same column) |

**What DockQ measures** (thresholds are DockQ's, not ours): a native contact is a binder-target residue pair with any heavy atoms within **5 Å** in the reference; the interface for iRMSD is every residue within **10 Å** of the partner; RMSDs use backbone **N, CA, C, O**. **LRMSD** superposes the model onto the reference on the **receptor** (the larger chain of the pair, here one target protomer) and takes the binder backbone RMSD without refitting: the target-fixed binder RMSD. It is per protomer, so a binder spanning two protomers gets two LRMSDs that differ with the protomer's own prediction error.

## Trust columns: read these before the science

| Column | Expected | How to read it |
| --- | --- | --- |
| `binder_identity` | **1.0** for a native-sequence control; **1 − n_mut/len** for a design | Fraction of binder positions identical to the reference. How many substitutions a design carries is a property of the batch (it varies between design arms): take `n_mut` from the designer's own column, or read the `dockq_binder_identity` distribution off your run. Near 1 − n_mut/len = the right binder was compared; far below = the wrong reference or chain |
| `binder_n_res` / `ref_binder_n_res` | equal | Unequal = the wrong reference (already an error in `number` mode) |
| `binder_n_res_compared` | `binder_n_res` in `number` mode; **the identical positions only** in `sequence` mode | How many binder residues DockQ used |
| `target_max_mismatches` | **0** | Positions differing in the worst target pair **of the kept mapping**; a permutation exceeding `--allowed-mismatches` is never a candidate |
| `target_max_mismatches_any_mapping` | **0** | The same over every permutation DockQ produced, skipped ones included; non-zero with the kept value 0 means some permutations were excluded (`per_mapping` shows which) |
| `mapping` | a cyclic permutation for a C3 trimer, e.g. `B>A;C>B;A>C;D>D` | **model>reference**, in reference order. The letter in `per_interface` and `best_ref_chain` is the REFERENCE chain; the model chain is `model=` in each entry. Boltz assigns A/B/C to copies arbitrarily, so a non-identity mapping is normal |
| `mapping_full_objective` / `mapping_agrees` | usually True | What DockQ's OWN objective (the sum over ALL native interfaces, target-target included) would have picked, and whether it is the same. **False is a finding, not an error**: the target-target interfaces outvoted the binder. Read `dockq_best_full_objective` beside it |
| `n_mappings_tried` | **6** for a trimer (3!) | Target permutations from DockQ's `group_chains`/`get_all_chain_maps`; 1 = the targets did not cluster as identical copies |
| `per_mapping` | one entry per permutation, every one with the same `n_binder=` | `<map, ','-joined>:target_mm=:n_binder=:binder_sum=:full_sum=`, or `:skipped`, `;`-joined. **Empty** = failed before the first permutation; **shorter** on an error row = the permutations computed before the failure. On a C3 trimer the three cyclic permutations score a full sum near the maximum (typical 2.8–2.9) and the three reflections far lower (~0.2): the symmetry, visible |
| `min_native_contacts` / `n_interfaces_used` | 1 / = `n_native_binder_interfaces` at the default | The `--min-native-contacts` echo, and how many native binder interfaces reach it (fewer used than native = some excluded as incidental) |
| `n_native_binder_interfaces` | 1–3 | Target chains the reference binder contacts: the number of `per_interface` entries and of terms in every candidate's binder sum |
| `model_chain_ids` / `ref_chain_ids` | e.g. `A,B,C,D` | Every protein chain in each file (auth_asym_id, HETATM excluded) |
| `ref_path` | the parent complex | The reference this row was scored against; one distinct value per parent |
| `dockq_version` | `2.1.3` | `unknown` if the metadata could not be read; numbers unaffected |
| `residue_match`, `allowed_mismatches`, `min_binder_identity`, the four chain echoes | as passed | Rows are comparable only at equal values |

### The chain mapping, and why it is not DockQ's default

DockQ v2 maps model to reference chains **by sequence**: two chains may correspond only if their alignment has at most `allowed_mismatches` mismatches (default **0**). At that default **a single substitution is enough to match nothing**, and DockQ exits with `For chains ['D'] no identical corresponding chain was found`. Raising `allowed_mismatches` enough to admit a batch's substitutions (~18 for a 30% variant of a 60-mer) would work but applies to every chain pair, inviting the binder to cluster with a target chain, and the needed number is a property of the batch. So the tool does what DockQ's `--mapping D*:D*` does:

1. **The binder is mapped explicitly** (reference binder ← model binder); its relatedness is checked by `--min-binder-identity` and, in `number` mode, by identical numbering.
2. **Only the target chains go through DockQ's own `group_chains`** at `--allowed-mismatches`. Every permutation `get_all_chain_maps` yields is checked: one whose target pairs differ beyond the allowance is skipped, each remaining one is scored with `run_on_all_native_interfaces`.
3. **Every scored permutation must yield the SAME set of native binder interfaces**, or the row errors: the homo-oligomer symmetry assumption made checkable, which also guarantees the binder sums have the same number of terms, so no permutation wins on count.
4. **The permutation kept maximises the summed DockQ of the BINDER interfaces only**, where DockQ's CLI maximises the sum over ALL native interfaces including target-target ones. In a symmetric trimer those are nearly equal across cyclic permutations, so they mostly do not move the choice; when the binder interfaces are all poor, small target-target differences can decide it. DockQ cannot restrict its objective to chosen interfaces (`--mapping :ABD` restricts chains, not interfaces), hence the tool's own selection. Ties go to the higher full sum, then to DockQ's enumeration order.

**The kept mapping maximises the binder-interface SUM, not the best single interface: the summed binder DockQ can only flatter the design, while `dockq_best` can come out BELOW DockQ's own choice** (kept 0.5/0.5, sum 1.0, best 0.5, vs DockQ's 0.7/0.2, sum 0.9, best 0.7). `dockq_best_full_objective` is the `dockq_best` DockQ's own choice would give, so it can sit above or below `dockq_best`. **Near the 0.23 gate, a row with `mapping_agrees = False` must be read with `dockq_best_full_objective`:** if the two straddle the bar, the pass depends on the mapping rule, not on a clear pose.

### Residue correspondence, and why the default differs from DockQ's

Inside a mapped chain pair, **DockQ's default aligns by sequence and keeps only IDENTICAL positions** (`get_aligned_residues` keeps `|` columns only). For a design every substituted residue is silently dropped from iRMSD and LRMSD and can never count as a reproduced contact, while Fnat's denominator is all native contacts of the **full** reference (`calc_DockQ` computes `nat_total` on the untouched native), so **a design cannot reach Fnat 1.0 even in a perfect pose**, and the cap grows with the number of mutated interface residues (which for a designed-His binder sit exactly at the interface). Comparing designs with different mutation counts would confound pose with mutation count.

**`--residue-match number` (the default here, not DockQ's) passes DockQ's own `no_align=True`:** residues correspond **by residue number** and every position is compared. iRMSD and LRMSD are backbone only, so a substitution does not change what they compare. **Fnat does NOT escape the mutations even here:** a contact is any HEAVY-ATOM pair within 5 Å (`get_residue_distances(all_atom=True)`), side chains included, so a design that shrinks an interface side chain can lose a native contact with a perfect backbone, and one that grows a side chain can gain contacts; `number` mode mixes pose with mutation less than `sequence` mode (where a substituted residue cannot count at all). **The identical-numbering guard is what makes `no_align` exact:** the tool requires the residue-number list of every mapped chain pair to be identical, or the row errors (otherwise DockQ would align the encoded numbers with gaps and silently drop residues). Boltz numbers every chain from 1: if a reference keeps its original numbering (e.g. 77–233) while the prediction is renumbered, this mode errors rather than guessing; use `sequence` mode or renumber the reference.

## The question it answers, and the bar

**Did the predicted binder keep the parent's binding mode?** DockQ's classes are the CAPRI conventions: **< 0.23 incorrect, ≥ 0.23 acceptable, ≥ 0.49 medium, ≥ 0.80 high** (DockQ v2 README). The usual bar is **acceptable or better on the best interface**: `dockq_dockq_best >= 0.23`. These are conventions for docking-pose quality, not binding: an acceptable pose can still have a different contact set, and a high one says nothing about affinity.

"Best" is the native binder interface with the highest DockQ among those with at least `--min-native-contacts` native contacts. `fnat_best`, `irmsd_best`, `lrmsd_best` and `best_ref_chain` all describe **that same interface**, never the best of each metric. An exact tie goes to the alphabetically first reference chain; **the `best_ref_chain` letter is not special**, since Boltz's lettering of identical copies is arbitrary. **Read `nat_best` beside it** (a best interface with 1–2 native contacts can clear 0.23 for a wrong pose); **choose `--min-native-contacts` from the per-interface `nat` counts of the parents** (`per_interface` of a native-sequence control shows them). For a binder spanning two protomers also read `dockq_mean_natw` or `n_interfaces_acceptable`: one protomer reproduced and the other lost means the binder pivoted on one contact patch.

## Invocation and flags

```bash
sapia run dockq <run_dir> -t <table> -i boltz_path \
    --binder-chain D --target-chains A,B,C --ref-binder-chain D --ref-target-chains A,B,C
sapia collect dockq <run_dir> -t <table>
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `-i/--input-column` | **none: required** (sentinel `not applicable` raises, listing the `_path` columns) | The MODEL complex column, resolved up the lineage |
| `--ref-col` | `structure_path` | The REFERENCE complex column, resolved up the lineage. A column this table HAS is read from the row's own cell and never falls back to an ancestor (the hisbond rule). Must differ from `-i` (refused at submit) |
| `--binder-chain` | **none: required** | Model binder, exactly one chain |
| `--target-chains` | **none: required** | Model target chains, mini-language `A,B,C` or `A:C`; must not contain the binder; no duplicates |
| `--ref-binder-chain` | **none: required** | Reference binder, exactly one chain |
| `--ref-target-chains` | **none: required** | Reference target chains, same count as `--target-chains`; their order does NOT fix the mapping (it is searched) |
| `--allowed-mismatches` | `0` | DockQ's `allowed_mismatches`, applied to the **target** chains only |
| `--residue-match` | `number` | `number` = DockQ `no_align` (by residue number, all positions, numbering must match); `sequence` = DockQ's own default (identical aligned positions only) |
| `--min-binder-identity` | `0.5` | Identity floor for the explicitly mapped binder (30% substitutions = 0.7; an unrelated parent's binder lands at ~0.3) |
| `--min-native-contacts` | `1` (DockQ's own rule, no exclusion) | Interfaces with fewer native contacts are excluded from `dockq_best`, `dockq_mean_natw`, `n_interfaces_acceptable` and `dockq_best_full_objective` and kept in `per_interface` as `used=0`. The excluded set depends only on the reference, not the mapping; the mapping objective still sums every native binder interface. No interface at or above N = error row |
| `--designs-per-task` | `20` | About 0.4 s and ~150 MB per design (1 CPU) |

Resources: 1 CPU, 2 GiB, 20 min (Modal `RESOURCES` and the matching `#SBATCH` lines). Submit-time refusals (nothing queued; an empty leaf dir and `.meta.json` remain, no phantom table since it is an `update` tool): binder in its own target set (model or reference), more than one binder chain, a duplicate chain, unequal target counts, `-i` equal to `--ref-col`, no `-i`, `--min-binder-identity` outside [0, 1], negative `--allowed-mismatches`, `--min-native-contacts` below 1, and every row resolving no model or no reference.

The builder prints `dockq: N design(s) against M distinct reference structure(s) from 'structure_path'`. **M should equal the number of parent complexes**; M = 1 on a two-parent table means the reference column did not resolve per lineage.

**Reference = the PARENT's own model, on a child table that has its own model.** The natural `--ref-col alphafold3_path` is wrong twice: the child table has that column, so it is read from the child's own cell (the model scored against itself, DockQ ≈ 1.0, spectacular and meaningless), and the submit refuses it anyway (`-i` equals `--ref-col`). Use a column the **child lacks and the parent has** (e.g. `dockq_<label>_path` from an earlier dockq run on the parent holds the parent's model path), which resolves by lineage to the parent's model. Prove it on a few rows BEFORE the run (resolved reference equals the parent's model path, not the child's); after it check `dockq_ref_path != dockq_path` on every row, that no DockQ is exactly 1.0, and that the distinct reference count equals the parent count.

## Gotchas, each with its signature

- **DockQ's default refuses a design.** The CLI on a design without `--mapping` prints `Model chain D, native chain D: 8 mismatches ... no identical corresponding chain was found` and exits: DockQ working as designed. The tool avoids it by mapping the binder explicitly.
- **DockQ calls `sys.exit()` on several failures.** The worker converts `SystemExit` into an error row (`error: DockQ group_chains (...) exited`); otherwise one bad design would kill the task.
- **DockQ caches results on Biopython chain objects** (`align_chains`, `get_residue_distances`, `get_aligned_residues`, `subset_atoms`, `run_on_chains`, `list_atoms_per_residue` are `lru_cache`d), so the worker clears every cache in the module before each design and names what it cleared in its log line (`caches cleared between designs: ...`). *Signature:* a function missing from that list means the package changed.
- **`sequence` mode under-reads designs.** *Signature:* `binder_n_res_compared < binder_n_res` and a `fnat_best` that never reaches 1.0 even when `lrmsd_best` is small.
- **`mapping_agrees = False`.** DockQ's own objective would have scored a different permutation; the reported numbers are for the binder-optimal mapping, `dockq_best_full_objective` gives DockQ's and `per_mapping` every permutation's sums.
- **A blank OWN model cell gives an EMPTY `dockq_status`, not `error:`** (`filter_ready` drops the row at collect; the task did write `error: could not resolve ...`). Count non-empty statuses against the row count.
- **`-f` acts on the run only.** Rows it leaves out are never submitted, but collect still iterates every ready row and marks the unsubmitted ones `missing` under that label. *Signature:* `missing` rows exactly where the filter excluded them; run a table that is already the subset you want, whole.
- **The package is an sdist only** (`DockQ==2.1.3` has no wheel; pip compiles its Cython extension), so the Modal image needs `build-essential`. If the extension is missing, DockQ's fallback import fails and every row reads `error: could not import the DockQ package: ...`; it never returns different numbers.
- **A hand check of `lrmsd_best` must use DockQ's chain mapping and its own fit.** A fit pairing chain A with A (B, C likewise) is wrong by tens of Å when the model chains are cyclically permuted (`mapping` `C>A;A>B;B>C`); superposing only the best-interface receptor protomer on N,CA,C,O with the table's mapping and `best_ref_chain` reproduces `lrmsd_best` to 4 decimals (a whole-trimer CA fit differs by 0.04–0.30 Å, a different superposition). Check a passing design, not only failing ones.
- **A second collect on the same table re-serialises the earlier label's columns as text** (integer-valued `1` to `1.0`, floats at 15 significant digits, change ~1e-15): read DockQ columns with `pd.to_numeric` and never compare text across collects.
- **Crashed worker:** every design of that task gets `error: worker crashed` from `dockq.sh`'s fallback.
- **An OK row that fails its own arithmetic is an error row.** The collector re-parses `per_interface` (full repr precision, so rounding cannot move a value across 0.23 or change a tie) and requires: each entry's `used` = (`nat` ≥ `min_native_contacts`); over used entries `dockq_best` = the largest DockQ, `best_ref_chain` and `nat_best` = that entry's, `dockq_mean_natw` = sum(dockq_i·nat_i)/sum(nat_i) (both to 1e-6), `n_interfaces_acceptable` = the count ≥ 0.23; entry count = `n_native_binder_interfaces`; used count = `n_interfaces_used`. *Signatures:* `error: dockq_best ... != max entry ...`, `error: per_interface entry '...' ...`, `error: dockq_mean_natw ... != ...`, `error: best_ref_chain ... != entry ...`. An OK row must also carry every science column plus `mapping`, `mapping_full_objective`, `mapping_agrees`, the binder identity/length facts, `target_max_mismatches`, `n_native_binder_interfaces`, `n_interfaces_used`, `min_native_contacts` and `ref_path`, or it reads `error: OK row missing ...` (descriptive metadata such as `dockq_version` and chain-id lists is not in that gate).

## NA is never 0

An `error: ...` or `missing` row has every science column (`dockq_best`, `fnat_best`, `irmsd_best`, `lrmsd_best`, `class_best`, `best_ref_chain`, `nat_best`, `dockq_mean_natw`, `n_interfaces_acceptable`, `dockq_best_full_objective`, `per_interface`) **NA**, and keeps the echo of the flags and whatever trust facts were established before the failure (e.g. unequal `binder_n_res`/`ref_binder_n_res` on a cross-parent row; `per_mapping` when the interface sets disagreed).

**A native interface the model does not reproduce is NOT NA.** DockQ scores it (Fnat 0.0, iRMSD and LRMSD of 10–30 Å, DockQ ~0.03): its score is never exactly 0 even for a binder on the far side of the target, because its RMSD terms are 1/(1+(rms/d0)^2). "Not reproduced" reads as `fnat=0.0` with a DockQ below ~0.05, not a zero DockQ. An interface that is not native (the reference binder does not touch that protomer) is not scored and does not appear in `per_interface`; a model binder that newly contacts it goes unpenalised except through its LRMSD.

## What it collects (`dockq_` prefix, or `dockq_<label>_` with `-l`)

Every column of the trust table is collected under the same prefix; read those first.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `dockq_best` | highest DockQ over native binder interfaces | the gate; ≥ 0.23 acceptable |
| `fnat_best`, `irmsd_best`, `lrmsd_best` | Fnat, iRMSD (Å), LRMSD (Å) **of that same interface** | LRMSD = binder backbone RMSD with that protomer superposed |
| `class_best` | `incorrect` / `acceptable` / `medium` / `high` | CAPRI class of `dockq_best` |
| `best_ref_chain` | reference target chain of the best interface | a REFERENCE letter; ties alphabetical; not special |
| `nat_best` | native contacts of the best interface | 1–2 = incidental: Fnat can only be 0, 0.5 or 1 there, so do not trust a pass carried by it |
| `dockq_mean_natw` | DockQ mean over the USED native binder interfaces, **weighted by native contact count** | a 2-contact incidental interface cannot count as much as a 30-contact real one; the plain mean is recomputable from `per_interface` |
| `n_interfaces_acceptable` | used native binder interfaces with DockQ ≥ 0.23 | equal to `n_interfaces_used` = the whole (non-incidental) footprint kept |
| `dockq_best_full_objective` | `dockq_best` under DockQ's OWN mapping choice | equals `dockq_best` when `mapping_agrees`; otherwise above OR below it |
| `per_interface` | `;`-joined, reference-chain order, full precision: `A:dockq=0.9479342446364449:fnat=1.0:irmsd=0.6311189037543716:lrmsd=0.6489641488883116:class=high:nat=13:model=B:used=1` | `nat` = native contacts, `model` = model chain mapped onto that reference chain, `used=0` = excluded by `--min-native-contacts` |
| `_status` | `OK` / `error: ...` / `missing` | `OK` is the only proof of success, and an OK row has passed the arithmetic cross-checks |
| `_path` | the MODEL structure scored | |

On disk: `<run_dir>/<table>/dockq[_<label>]/<name>.tsv` per design (full precision) and `dockq_tasks/task_<n>.tsv` (name, model, reference) per task.

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
        & (df[P + "nat_best"] >= MIN_NAT)   # a pass carried by an incidental interface is not a pass
    )
    # A pass that exists only under the binder-optimal mapping is not a clear pass:
    # require DockQ's own mapping to agree, or its dockq_best to clear the bar too.
    robust = (df[P + "mapping_agrees"] == True) | (  # noqa: E712 - column is object dtype
        df[P + "dockq_best_full_objective"] >= 0.23
    )
    return df[ok & robust]
```

`mapping_agrees` comes back as an object column (True/False/NaN): compare with `== True`. With `-l v1` every column is `dockq_v1_<name>`, hence the interpolated prefix. A forked table (`--table-label` upstream) changes the **table name** you pass to `-t`, not the column names. For a design-set comparison read `dockq_dockq_best` and `dockq_lrmsd_best` per set (a distribution question for `table-analyst`), always beside the native-sequence control of the same parent and recipe.

## What it is blind to

- **Side chains and chemistry.** Backbone RMSDs plus residue-level contacts. Whether the designed His reaches the target: `hisbond`; interface energy and packing: `pyrosetta` (`if_dG`) and `cms` (CMS, SC).
- **The reference's own accuracy.** The parent is a prediction; "high" means the same pose as the parent.
- **Confidence.** A pose can reproduce the parent with low pLDDT/ipSAE, or miss it confidently: read `boltzconf` / `af3conf` beside it.
- **A different pose on the same footprint.** A binder rotated or flipped in place that still touches the same target residues earns partial Fnat while its centroid moves only a few Å and its binder CA RMSD exceeds 10 Å: read `lrmsd_best` beside `fnat_best`.
- **New, non-native contacts** to a protomer the reference binder never touched: not scored as an interface.
- **Binding.** CAPRI classes grade docking poses, not affinity or specificity.
- **Whole-complex agreement.** For a TM-score/RMSD over the whole complex use `usalign`; its numbers are a different quantity and need not equal LRMSD.
