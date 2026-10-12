---
name: pyrosetta
description: >-
  PyRosetta energy scoring on Modal: FastRelax by default, plus interface metrics for complexes. Load before composing a pyrosetta run. Traps: choose the structure column to score deliberately, know when to turn relax off, and failures exit 0.
---

# pyrosetta

Scores one structure per design with [PyRosetta](https://www.pyrosetta.org/) (ref2015 by
default), after an optional FastRelax, and writes Rosetta energy metrics as columns.
**`action: update`** — annotates the table it reads, in place. `-t` is required.

Source: `prosapia/tools/pyrosetta/` in the installed package
(`.venv/lib/python<version>/site-packages/`), or `src/prosapia/tools/pyrosetta/` checkout
 — `run_pyrosetta.py`, `pyrosetta_worker.py`, `collect_pyrosetta.py`. 
`sapia run pyrosetta --help` is authoritative for flags.

> **Validated end to end on Modal.** `outputs/20260928_113443_pyrosetta_test` on
> `sapia-runs`: rfd3 → ProteinMPNN → Boltz → pyrosetta (defaults) scored 2 designs, both
> `OK`, submit ~8 s, tasks done in ~45 s. The worker was also checked on ubiquitin (1UBQ)
> and barnase–barstar (1BRS), including `--interface A_D` / `auto`, constrained relax and
> the missing-file path. `--interface` has not yet been run through `sapia run` on a
> designed complex.

## What to score

It reads **`-i/--input-column`**, default **`boltz_path`** — a full-atom predicted
structure. mmCIF is fine: structures are converted to PDB at manifest-build time, in the
workstation, cached under `<run_dir>/.cif_to_pdb/`.

- **Score predictions, not rfd3 backbones.** An rfdiffusion3 backbone carries no designed
  sequence (only placeholder residues), so its energy means nothing. The usual target
  is the Boltz (or AF3) prediction of a ProteinMPNN sequence, on the sequence table.
- Unlike usalign, `-i` is **not** resolved up the lineage — the column must be on the
  table you pass with `-t`.

## Canonical runs

```bash
# After boltz on table1: 1 FastRelax cycle, then score (default)
sapia run     pyrosetta <run_dir> -t table1
sapia collect pyrosetta <run_dir> -t table1

# Score as-is, into its own column family (pyrosetta_raw_*) alongside the relaxed one
sapia run     pyrosetta <run_dir> -t table1 --relax-cycles 0 -l raw
sapia collect pyrosetta <run_dir> -t table1 -l raw

# A binder / complex: add interface metrics (chain A vs chain B)
sapia run     pyrosetta <run_dir> -t table2 --interface A_B
sapia collect pyrosetta <run_dir> -t table2
```

**Collect takes only `-t` and `-l`** — the metrics set comes from what the worker wrote.
If you ran with `-l <label>`, collect with the same `-l`, or it reads the wrong dir.

## Flags

| Flag | Default | Note |
| --- | --- | --- |
| `-i` | `boltz_path` | Structure column to score. |
| `--relax-cycles` | `1` | FastRelax repeats before scoring. `0` = score as-is. `5` is Rosetta's standard relax (≈5× slower). |
| `--constrain-relax` | off | Coordinate-constrain relax to the input, so it fixes clashes/rotamers without moving the backbone. Use it when the structure itself is what you are judging. |
| `--interface` | none | InterfaceAnalyzer spec, `A_B` / `AB_C`, or `auto` (first chain vs the rest). Skipped when empty. |
| `--scorefxn` | `ref2015` | Any Rosetta score function name. The term columns follow it. |

Why relax by default: a Boltz/AF3 structure has never seen the Rosetta energy function, so
unrelaxed scores are dominated by `fa_rep` clashes and `fa_dun` rotamer penalties — they
measure the predictor, not the design. One cycle removes most of that.

Resources: **1 CPU, 4 GiB, 1 h timeout, no GPU.** The builder forces `gpus_per_task = 0`,
so **`-g 0` is not needed.** Rosetta is single-threaded; `-c` above 1 buys nothing.
One task per design.

**Image**: `pyrosetta-installer` downloads a ~1.5 GB wheel from RosettaCommons at build
time. It is **already built and cached** on this Modal workspace (the install step took
~2 min; one earlier cold build took far longer). If `modal_image.py` changes it rebuilds
inside the submit call — give that Bash call a 15+ min timeout before assuming it is stuck.

**Timings measured (1 CPU)**: ~5 s to score a 76-aa protein as-is, ~10 s with 1 relax
cycle; ~3 min for a 588-residue complex with 1 cycle. PyRosetta init is ~3–4 s of that.
Budget roughly linear in size × cycles; `--relax-cycles 5` on a large complex can approach
the 1 h timeout — raise `-T` for those.

## Gotchas

- **A failed design still exits 0.** The worker catches every error (missing file,
  unreadable structure, bad `--interface` chains) and writes `error: <reason>` into the
  design's TSV; if the worker itself dies, the `.sh` writes `error: worker crashed`. So all
  `.exit` = `0` only proves the tasks **ran** — `pyrosetta_status` after collect is what
  proves they worked. Check it and report the error strings.
- **`missing` at collect** = no TSV on disk: the task never ran or was killed (timeout on a
  large complex with `--relax-cycles 5`) — check the app per the worker's `.exit`
  loop, and rerun with `-T`.
- **Reruns skip `pyrosetta_status == OK` rows** unless `--force`. A relaxed and an
  unrelaxed pass share a leaf unless you separate them with `-l`; without it, the second
  run submits nothing ("0 designs").
- **Don't trust interface energies without relax.** InterfaceAnalyzer repacks the
  *separated* partners but not the complex, so on an unrelaxed structure `if_dG` comes out
  **positive** from the complex's own clashes (1BRS as-is: +210 REU; after 1 cycle: −58.6).
  Keep `--relax-cycles ≥ 1` whenever `--interface` is set.
- **Error strings are truncated to 300 characters** (PyRosetta's pybind errors list every
  overload). Enough to diagnose; re-run the worker by hand if you need the rest.
- **Non-canonical residues / ligands** are dropped on load (`-ignore_unrecognized_res`),
  so a ligand-bound Boltz model is scored as apo protein. Say so if it applies.

## What it collects

Columns are `pyrosetta_<field>` (with `-l <label>`, `pyrosetta_<label>_<field>`).

| Column | Meaning |
| --- | --- |
| `_status` | `OK`, `error: <reason>`, or `missing`. The only success signal. |
| `_path` | The scored (relaxed) PDB: `<run_dir>/<table>/pyrosetta/<name>.pdb`, run-relative like every other `_path`. |
| `_total_score` | Total energy of the scored pose (REU). Scales with size — compare `score_per_res` across lengths. |
| `_score_per_res` | `total_score / n_res`. |
| `_score_raw` | Total energy **before** relax (equals `total_score` with `--relax-cycles 0`). |
| `_relax_ca_rmsd` | CA RMSD between input and relaxed pose (Å). |
| `_n_res`, `_n_chains` | Pose size. |
| `_fa_atr`, `_fa_rep`, `_fa_sol`, `_fa_elec`, `_hbond_*`, `_rama_prepro`, `_fa_dun`, `_p_aa_pp`, `_ref`, … | Every nonzero-weighted term of the score function, **already weighted** (they sum to `total_score`). |
| `_sasa` | Total SASA (Å²). |
| `_sasa_hydrophobic` | Hydrophobic SASA (Å²). |
| `_packstat` | RosettaHoles packing score, 0–1, higher = better packed. Slightly stochastic. |
| `_buried_unsat` | Buried unsatisfied polar atoms (BuriedUnsatHbondFilter). |
| `_dssp` | DSSP string (H/E/L per residue). |
| `_if_dG`, `_if_dSASA`, `_if_dG_per_dSASA`, `_if_hbonds`, `_if_delta_unsat`, `_if_n_res`, `_if_sc`, `_if_packstat` | Only with `--interface`. `if_dG` is the separated-minus-bound energy (REU, more negative = stronger); `if_dG_per_dSASA` is ×100; `if_sc` is shape complementarity (0–1). |

On disk: `<run_dir>/<table>/pyrosetta/<name>.tsv` (one row) and `<name>.pdb`.

## Reading the result (for the report back)

Rough ref2015 conventions for de-novo designs — rank within a batch, don't treat as hard
cutoffs:

- Reference points measured with this tool: ubiquitin (a natural, well-packed 76-aa
  protein) scores **−3.4 REU/res**, packstat 0.67, 10 buried unsats after 1 cycle;
  unrelaxed it is +0.43 REU/res — the relax gap is normal for a non-Rosetta structure.
- `score_per_res` **≤ −2.0** is typical of a well-packed, well-designed monomer; near 0 or
  positive means clashes or a poor core.
- `fa_rep` large and positive after relax → a clash the relax could not resolve.
- `relax_ca_rmsd` > ~1.5–2 Å with one cycle → the predicted structure was not near a
  Rosetta minimum; distrust the energies of that design.
- `packstat` ≳ 0.6 is well packed; `buried_unsat` should be small (a few at most for ~100
  aa).
- Interfaces: `if_dG` ≲ −30 REU with `if_sc` ≳ 0.6 is the usual binder bar (barnase–
  barstar, a femtomolar complex, gives −58.6 REU, `if_sc` 0.76, dG/dSASA ×100 −3.5); use
  `if_dG_per_dSASA` to compare interfaces of different size.

Report energies **together with** confidence (Boltz pLDDT) and self-consistency (usalign
RMSD): a design that folds back to its backbone *and* scores well is the one to keep.
