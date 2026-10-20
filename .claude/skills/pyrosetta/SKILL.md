---
name: pyrosetta
description: >-
  How does a structure score under the Rosetta energy function after a FastRelax, and for a binder/target complex how strong and well packed is the interface? PyRosetta scoring on Modal and vib. Load before composing a pyrosetta run or reading any `if_*` number. A project-local FORK in `tools/pyrosetta/` shadows the bundled tool. Traps: the `--interface` spec must be explicit and is checked against the pose, the interface step runs AFTER the relax, the choice of structure column matters, and failures exit 0.
---

# pyrosetta

Scores one structure per design with [PyRosetta](https://www.pyrosetta.org/) (ref2015 by default), after an optional FastRelax, and writes Rosetta energy metrics as columns. **`action: update`**, annotating the table it reads in place; **`-t` is required**. CPU only: the builder forces `gpus_per_task = 0`, so `-g 0` is not needed. Runs on both executors: on vib PyRosetta 2025.37 (py3.13, conda env `pyrosetta`) is activated by `$SAPIA_VIB_SHARED_ACTIVATION/pyrosetta.sh` via `SAPIA_ACTIVATE_PYROSETTA` in `.env`. Source: **`tools/pyrosetta/`** (the fork below, which is what runs); the bundled original is `prosapia/tools/pyrosetta/` in the installed package. `sapia run pyrosetta --help` is authoritative for flags; the run loop is in `running-a-step`.

## A project-local FORK shadows the bundled tool

`tools/pyrosetta/` is a fork of the bundled tool and **wins on tool discovery** (same `Tool.name`, so columns and output leaf are unchanged: `pyrosetta_*`, `<table>/pyrosetta/`). It exists because the bundled `interface_metrics` goes through `rosetta.core.pose.DockingPartners`, which **does not exist** in PyRosetta 2025.37 (the vib conda build): every `--interface` run raised `module 'pyrosetta.rosetta.core.pose' has no attribute 'DockingPartners'` after the relax had already run. One functional change plus three guards, all in `interface_metrics` / `_interface_spec`; relax, `score_metrics`, the columns and error-as-data behaviour are the bundled ones.

- **The analyzer is built from the dock-chains STRING**, `InterfaceAnalyzerMover("<chains>_<chains>")`.
- **`--interface auto` is refused for more than two chains** (it meant `A_BCD` on a four-chain complex: a believable number for an interface nobody asked for).
- **An explicit spec is validated against the pose**: exactly one underscore, both sides non-empty, no repeated letter, every letter a chain that exists in the structure, and the two sides together covering **every** chain of it. Otherwise the task raises naming the offending chains.
- **`pyrosetta_if_partners`** records the spec that produced the numbers.

## Premise, and when its numbers are meaningless

It reads **`-i/--input-column`** (default `boltz_path`): a full-atom structure with a designed sequence. mmCIF is fine (converted to PDB at manifest-build time, cached under `<run_dir>/.cif_to_pdb/`). Unlike usalign, `-i` is **not** resolved up the lineage: the column must be on the `-t` table.

- **Score predictions, not rfd3 backbones** (placeholder residues, so the energy means nothing): the usual target is the Boltz or AF3 prediction of a ProteinMPNN sequence, on the sequence table.
- **Unrelaxed predictor output measures the predictor**: it has never seen the Rosetta energy function, so unrelaxed scores are dominated by `fa_rep` clashes and `fa_dun` rotamer penalties. One FastRelax cycle removes most of that, hence the default.
- **Interface energies need a relax.** InterfaceAnalyzer repacks the *separated* partners but not the complex, so on an unrelaxed structure `if_dG` comes out **positive** from the complex's own clashes (barnase-barstar as-is: +210 REU; after one cycle: -58.6). Keep `--relax-cycles >= 1` whenever `--interface` is set.
- **Non-canonical residues and ligands are dropped on load** (`-ignore_unrecognized_res`): a ligand-bound model is scored as apo protein; say so if it applies.

## Canonical runs

```bash
# 1 FastRelax cycle, then score (default)
sapia run     pyrosetta <run_dir> -t <table>
sapia collect pyrosetta <run_dir> -t <table>

# Score as-is into its own column family (pyrosetta_raw_*)
sapia run     pyrosetta <run_dir> -t <table> --relax-cycles 0 -l raw

# A binder / complex: spec ALWAYS explicit, covering every chain (3-chain target A,B,C vs binder D)
sapia run     pyrosetta <run_dir> -t <table> --interface ABC_D

# Test a spec in seconds without the relax (the interface step runs after it)
sapia run     pyrosetta <run_dir> -t <table> --interface ABC_D --relax-cycles 0 -l spectest
```

**Collect takes only `-t` and `-l`** (the metric set comes from what the worker wrote); with `-l <label>` on the run, collect with the same `-l` or it reads the wrong dir.

## Flags

| Flag | Default | Note |
| --- | --- | --- |
| `-i` | `boltz_path` | Structure column to score. |
| `--relax-cycles` | `1` | FastRelax repeats before scoring; `0` scores as-is; `5` is Rosetta's standard relax (about 5x slower). |
| `--constrain-relax` | off | Coordinate-constrain the relax to the input: fixes clashes and rotamers without moving the backbone. Use it when the structure itself is what you are judging. |
| `--interface` | none | InterfaceAnalyzer spec, **explicit**: `A_B`, `ABC_D`. Must name only chains of the scored structure and cover **all** of them. `auto` is accepted only for a two-chain structure. Skipped when empty. |
| `--scorefxn` | `ref2015` | Any Rosetta score function name; the term columns follow it. |

Resources: **1 CPU, 4 GiB, 1 h timeout, no GPU** (the `.sh` matches); Rosetta is single-threaded, so `-c` above 1 buys nothing. One task per design. Typical (1 CPU): ~10 s for a 76-aa protein with one cycle; ~3 min for a 588-residue complex on Modal; ~6.5-7 min and ~1 GB for a 535-residue 4-chain complex with one constrained cycle on vib. Roughly linear in size x cycles: `--relax-cycles 5` on a large complex can approach the timeout, so raise `-T`. **(Modal only)** the image installs a ~1.5 GB wheel at build, inside the submit call on first use or after any `modal_image.py` change: give that call 15+ minutes.

## Traps, each with its signature

- **A failed design still exits 0.** The worker catches every error (missing file, unreadable structure, bad `--interface` chains) and writes `error: <reason>` (truncated to 300 characters; pybind errors list every overload); if the worker dies the `.sh` writes `error: worker crashed`. Zero task exit codes prove only that the tasks *ran*: `pyrosetta_status` after collect proves they worked. Report the error strings.
- **`missing` at collect** means no TSV on disk: the task never ran or was killed (timeout on a large complex). Rerun with `-T`.
- **Reruns skip `OK` rows** unless `--force`: a relaxed and an unrelaxed pass share a leaf unless separated with `-l`, and without it the second run submits nothing (`0 designs`).
- **The interface step runs AFTER the relax, and a failure there costs the whole row.** A bad spec (a chain not in the structure, an uncovered chain, `auto` on more than two chains) raises once the relax is done: the design becomes an `error:` row with **no metrics and no scored PDB**, and the relax time is spent. **Test a spec with `--relax-cycles 0 -l spectest` first** (seconds instead of minutes).
- **`--interface A_BCD` passes the guards and is NOT the same measurement as `ABC_D`** (tens of REU apart in `if_dG`); for a binder you want target chains vs binder. **Always read `pyrosetta_if_partners` before comparing `if_*` across rows or runs.**

## What it collects

Columns are `pyrosetta_<field>` (`pyrosetta_<label>_<field>` with `-l`). On disk: `<run_dir>/<table>/pyrosetta/<name>.tsv` (one row) and `<name>.pdb`.

| Column | Meaning |
| --- | --- |
| `_status` | `OK`, `error: <reason>`, or `missing`. The only success signal. |
| `_path` | The scored (relaxed) PDB, run-relative like every other `_path`. |
| `_total_score`, `_score_per_res`, `_score_raw` | Total energy of the scored pose (REU; scales with size, so compare `score_per_res` = `total_score / n_res` across lengths) and the total **before** relax. |
| `_relax_ca_rmsd` | CA RMSD between input and relaxed pose, Å (`0.0` with `--relax-cycles 0`). |
| `_n_res`, `_n_chains` | Pose size. |
| `_fa_atr`, `_fa_rep`, `_fa_sol`, `_fa_elec`, `_hbond_*`, `_rama_prepro`, `_fa_dun`, `_p_aa_pp`, `_ref`, ... | Every nonzero-weighted term of the score function, **already weighted** (they sum to `total_score`). |
| `_sasa`, `_sasa_hydrophobic` | Total and hydrophobic SASA (Å²). |
| `_packstat` | RosettaHoles packing score, 0-1, higher is better packed; slightly stochastic. |
| `_buried_unsat`, `_dssp` | Buried unsatisfied polar atoms (BuriedUnsatHbondFilter; scales with size); DSSP string (H/E/L per residue). |
| `_if_partners` | Only with `--interface`: the dock-chains spec these `if_*` numbers describe (`ABC_D`). **Read it before comparing interfaces.** |
| `_if_dG`, `_if_dSASA`, `_if_dG_per_dSASA`, `_if_hbonds`, `_if_delta_unsat`, `_if_n_res`, `_if_sc`, `_if_packstat` | Only with `--interface`. `if_dG` is the separated-minus-bound energy (REU, **more negative = stronger**); `if_dSASA` is the buried area summed over **both** sides; `if_dG_per_dSASA` is Rosetta's dG/dSASA ratio times 100 (`if_dG / if_dSASA * 100`, no other factor); `if_n_res` uses Rosetta's own interface definition, **wider than a 5 Å heavy-atom contact count**, so an independent contact count comes out smaller and the two are not interchangeable; `if_sc` is shape complementarity (0-1). |

## Invariants a worker can check

- **Real:** the term columns sum to `total_score` (separate energy lookups); `_n_res` and `_n_chains` match the structure; with `--interface`, `if_partners` equals the spec passed and covers every chain letter; `if_dG_per_dSASA` equals `if_dG / if_dSASA * 100`.
- **Tautologies:** `score_raw == total_score` and `relax_ca_rmsd == 0.0` when `--relax-cycles 0`.
- **Break-it control** (own `-l`): a spec with an absent chain letter, an uncovered chain, or `auto` on more than two chains must give `error:` rows with no metrics and no PDB.

## Reading the result

Rough ref2015 conventions for de-novo designs: rank within a batch, not hard cutoffs.

- Typical: a natural well-packed 76-aa protein (ubiquitin) scores about **-3.4 REU/res** after one cycle (+0.43 unrelaxed: a large relax gap is normal for a non-Rosetta structure), packstat 0.67, 10 buried unsats. `score_per_res` **<= -2.0** is typical of a well-packed designed monomer; near 0 or positive means clashes or a poor core; a large positive `fa_rep` after relax is a clash the relax could not resolve.
- `relax_ca_rmsd` above ~1.5-2 Å with one cycle: the structure was not near a Rosetta minimum, distrust its energies. `packstat` >= ~0.6 is well packed.
- Interfaces: `if_dG` <= ~-30 REU with `if_sc` >= ~0.6 is the usual binder bar (barnase-barstar: -58.6 REU, `if_sc` 0.76, dG/dSASA x100 -3.5); use `if_dG_per_dSASA` to compare interfaces of different size.
- Typical for a designed binder/target complex (AF3 model, 535 residues, 4 chains, one **constrained** cycle, explicit `ABC_D`): `relax_ca_rmsd` 0.5-0.6 Å, `score_per_res` about -3.44, `if_dG` -70...-78 REU, `if_sc` 0.62-0.65, `if_hbonds` 15-17, `if_dSASA` 2470-2540 Å², `if_n_res` 90-91, `buried_unsat` 61-67 (size-dependent, not comparable with a monomer's ~10), `packstat` 0.56-0.59. A design far outside such ranges on the same kind of input is the one to look at.

## What it is blind to

One relaxed conformation on the Rosetta scale (REU, not kcal/mol); `if_dG` is not an affinity. Nothing here says whether the structure is confident or folds back to its designed backbone: report energies **with** confidence (`boltzconf`) and self-consistency (`usalign`), since a design that folds back *and* scores well is the one to keep. Whether a designed histidine reaches the target is `hisbond`.
