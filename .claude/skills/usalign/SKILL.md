---
name: usalign
description: >-
  Structure-vs-structure comparison per design with USalign: TM-score, RMSD, alignment lengths. Load before composing a usalign run or reading any TM or RMSD it wrote. Traps: which column is honest depends on `--mm` and on monomer vs complex; on a complex the RMSD inverts, so low RMSD marks designs that drifted; it takes a `--col-a`/`--col-b` pair; failures exit 0.
---

# usalign

Compares two structures per design with [USalign](https://github.com/pylelab/USalign) and writes
TM-score / RMSD columns. **`action: update`** — annotates the table it reads, in place. `-t` is
required.

There is no `docs/tools/usalign.md`; the source is the reference
(`.venv/lib/python3.13/site-packages/prosapia/tools/usalign/`, or
`src/prosapia/tools/usalign/` in a `../prosapia` checkout), and `sapia run usalign --help` is
authoritative for flags. Everything mechanical below was read from that source.

## Premise, and when its numbers are meaningless

usalign superposes two structures **and then measures the residues that survived the
superposition.** Both halves matter, and the second is where campaigns have lost weeks.

**It always returns a number.** There is no input for which it reports "this comparison is
meaningless" — so *which of its columns is honest* is decided by `--mm`, by chain composition, and by
nothing the tool will tell you.

Three regimes. Know which one you are in **before** composing the flags:

| what you are comparing | mode | the honest column | what lies |
| --- | --- | --- | --- |
| one chain vs one chain (binder vs its designed backbone) | `--mm 0` | `_RMSD`, `_TM1`/`_TM2` | — |
| a complex vs the same complex (did the binder stay in its pose?) | `--mm 1` | **`_Lali − L_target`** | `_RMSD` **and** `_TM*` |
| a complex, read as if it were a fold score | `--mm 1` | **none** | all of them |

### The default is `--mm 1`, and on a complex the RMSD inverts

`--mm` defaults to **1** (multimer). Run it on a predicted complex against a reference complex and
two things happen at once:

1. **The large fixed target dominates the superposition**, so `_TM1`/`_TM2` read ~1.0 regardless of
   what the binder did.
2. **`--mm 1` keeps only the residues that structurally align.** A binder that drifted off its
   designed pose *drops out of the alignment*, leaving a **target-only** superposition — which
   reports a beautiful RMSD around **0.4 Å**.

So on a complex, **a low `_RMSD` is evidence the binder is gone.**

*Measured on `binder_nohis`, 818 rows:* `Spearman(_RMSD, Lali − L_target) = +0.945` — low RMSD
coincides with the drifted, target-only rows, so **ranking by low RMSD selects the worst designs in
the batch.** Over the same 818 rows `Lali − L_target` had **median 1 and max 67 against 62–88
residue binders**: the binder was off its designed pose for essentially the whole batch, while the
RMSD column looked excellent throughout.

**`Lali − L_target` is the honest column.** `_Lali` counts aligned residues over all chains; subtract
the fixed target length and what remains is **binder residues that aligned in the target frame**:

```
Lali − L_target  ≈ binder length   the binder stayed in its designed pose
Lali − L_target  ≈ 0               the binder is gone; any RMSD you read is target-only
```

Cross-check against `epitope`-on-the-prediction hotspot recall, which it tracked at
`Spearman = +0.671` on the same rows.

### `--ter` will not rescue you

`--ter 2` isolates the **first** chain, and in these complexes the binder is usually **last**. Use
`chainsel` to extract the binder and compare with `--mm 0` — that is a different gate (below), not a
fix for this one.

## The two gates it serves, and which threshold belongs to each

Fold and pose are independent. *Measured:* the best self-consistency in a campaign (TM 0.971,
RMSD 0.58 Å) was also its most impossible pose.

**Gate "did it fold" — binder chain only, `--mm 0`, against the designed backbone.** Extract the
binder with `chainsel` first. Here the conventional bar applies: **RMSD < 2 Å**, with TM above ~0.9
saying the predicted fold matches the designed one. A design with confident Boltz metrics but a poor
RMSD back to its own backbone is a failed **design**, not a failed prediction — report both numbers
together.

**Gate "did it stay" — the whole predicted complex vs the rfd3 backbone complex, `--mm 1`, read as
`Lali − L_target`.** Both hold the *same* target, so the superposition locks onto it and the binder's
displacement shows up in that frame.

```bash
sapia run usalign <run_dir> -t table1_biasC \
    --col-a boltz_cofold_path --col-b rfdiffusion3_path --mm 1 --output-prefix pose
```

**Do not carry the `RMSD < 2 Å` bar from the fold gate to the pose gate.** It is correct for the
first and inverted for the second, and that single substitution is the most expensive mistake
available with this tool. There is no dedicated pose tool — `framefit` was built and removed on
2026-10-02 — so this is the route. See `binder-campaign` for where both gates sit in the order.

## Trust columns — read these before any TM or RMSD

usalign aligns sequences as well as structures, and the identity columns tell you it aligned the
thing you meant.

| column | read it for | expect |
| --- | --- | --- |
| `_ID1`, `_ID2` | sequence identity normalised by A and by B | ~1.0 comparing a design to its own backbone; well below means you compared two different molecules |
| `_IDali` | identity over the aligned region only | ~1.0; high `_IDali` with low `_ID1` means only a fragment aligned |
| `_L1`, `_L2` | lengths of A and B | must match your expectation from the parent table; a surprise here is a wrong column, not a bad design |
| `_Lali` | residues that aligned | the pose-gate numerator; compare against `_L1`/`_L2` |

Generators renumber: **rfd3 renumbers every output chain from 1**, so a target numbered 18–155 comes
back as 1–138. The identity columns are how you notice the comparison still lined up.

## It does not use `--input-column`

Its `default_input_column` is the literal string `"not applicable"`. The manifest builder selects
rows on **`--col-a`** instead (rows where that column is non-empty), so `-i` is meaningless here.
Pass:

- **`--col-a`** — the column holding structure A. Required, and it also chooses which rows run.
- **`--col-b`** *or* **`--ref`** — structure B. Exactly one; both, or neither, raises at submit time.

**`--col-b` is resolved up the lineage.** For each row it takes that row's own value if set,
otherwise the nearest ancestor's (via `lookup`). That is what makes self-consistency possible: a
child row's prediction compared against the backbone its parent holds, with no column copied
forward.

`--ref` is a single fixed structure for every comparison, resolved to an absolute `/runs` path, so
**the file must already be on the Volume**.

```bash
# self-consistency: table1 holds the sequences and their predictions, table0 the parent backbones
sapia run     usalign <run_dir> -t table1 --col-a boltz_path --col-b rfdiffusion3_path --mm 0
sapia collect usalign <run_dir> -t table1 --col-a boltz_path --col-b rfdiffusion3_path
```

`boltz_path` lives on `table1`; `rfdiffusion3_path` does not — it is found on `table0` through the
parent link. Columns land as `usalign_boltz_vs_rfdiffusion3_*`.

## The prefix is the contract between run and collect

Every column and the results dir are keyed by a **prefix** derived from the flags:
`<col_a minus "_path">_vs_<col_b minus "_path">` (with `--ref`, the ref file's stem), or whatever
`--output-prefix` you pass.

**Collect recomputes it from its own flags**, so it must be given the same ones — which is why
`--col-a`/`--col-b` reappear on the collect line above. If they disagree, collect fails with
`Results dir not found: .../<prefix>`. **Passing `--output-prefix <name>` on both phases is the
safest form** whenever the comparison is anything non-obvious, and it is what lets one table hold
several comparisons side by side (`pose`, `fold`, `boltz_vs_openfold3`, …) without `-l/--dir-label`.

## Flags

| Flag | Default | Note |
| --- | --- | --- |
| `--col-a` | — | **Required.** Column for structure A; also selects which rows run. |
| `--col-b` | none | Column for structure B, resolved up the lineage. Mutually exclusive with `--ref`. |
| `--ref` | none | One fixed structure for all comparisons. Must be on the Volume. |
| `--output-prefix` | `<a>_vs_<b>` | Names the results dir and every column. **Must match at collect.** |
| `--mm` | **`1`** | USalign `-mm`: 0 monomer, 1 multimer, 2 chain-to-complex, 3 circular permutation, 4 >2 structures, 5 fully non-sequential, 6 semi-non-sequential. **For a single-chain comparison you must pass `--mm 0`** — the default is multimer; see the premise section for what that does to the numbers. |
| `--ter` | `0` | USalign `-ter`: 0 all chains in all models, 1 all chains of the first model, 2 first chain only, 3 first chain split at TER. **`--ter 2` takes the first chain; the binder is usually last.** |

Resources: **1 CPU, 4 GiB, 30 min timeout, no GPU.** The manifest builder forces
`gpus_per_task = 0` itself, so **you do not need `-g 0`** — passing it is harmless but redundant.

The image compiles USalign from source (`git clone` + `g++`), so the *first* submit builds it inside
the call — give that Bash call a long timeout like any other first run. After that it is cached.

Structures are converted CIF→PDB at **manifest-build time**, in the workstation, cached under
`<run_dir>/.cif_to_pdb/`. A large table therefore makes the submit call itself slow; the tasks are
trivial.

## Gotchas

- **A failed comparison still exits 0.** The task script catches a missing input, a USalign crash and
  unparsable output, writes `ERROR: <reason>` into that design's TSV, and exits `0`. So **`n_tasks`
  `.exit` files all `0` does not mean the step succeeded** — for this tool the `.exit` loop only
  proves the tasks ran. The truth is the `usalign_<prefix>_status` column after collect.
- **Rows can be dropped from the manifest silently.** If `--col-b` resolves to nothing on a row *and*
  all its ancestors, that design is skipped. `Submitting N designs` may be well under the table's row
  count — compare the two and say so if they differ.
- **`missing` at collect means no TSV on disk**, usually one of those skipped rows: collect iterates
  every row in the table, since the `"not applicable"` input column filters nothing out.
- **It writes no `usalign_status` / `usalign_path`.** Status is per-comparison
  (`usalign_<prefix>_status`), which suppresses the framework's usual leaf status stamp. So do not
  look for `usalign_status`, and note that collect never skips already-done rows and re-reads every
  TSV. Idempotent and safe.
- **`--force` on the run** re-submits every design, printing `Re-running USalign for all designs
  (including those with status OK).` Don't use it to fix a failed comparison until you have read the
  reason out of the status column.

## What it collects

Columns are leaf-prefixed *and* prefix-keyed: `usalign_<prefix>_<field>` (with `-l`, the leaf becomes
`usalign_<label>`).

| Column | Meaning | How to read it |
| --- | --- | --- |
| `_status` | `OK`, `ERROR: <reason>`, or `missing` | the only success signal for this tool |
| `_ID1`, `_ID2`, `_IDali` | sequence identity by A, by B, over the alignment | **trust columns — read first** |
| `_L1`, `_L2` | lengths of A and B | check against the parent table |
| `_Lali` | aligned length over all chains | **on a complex, `_Lali − L_target` is the pose answer** |
| `_TM1`, `_TM2` | TM-score normalised by A and by B | fold quality **only** under `--mm 0`; target-dominated under `--mm 1` |
| `_RMSD` | RMSD over the aligned residues, Å | fold quality under `--mm 0`; **inverted and misleading under `--mm 1`** |
| `_sup_path` | the superposition PDB USalign wrote | open it when a number surprises you |

On disk: `<run_dir>/<table>/usalign/<prefix>/<name>.tsv` plus `<name>.pdb` (the superposition).
*(The collector's docstring names an older layout — trust this one.)*

## Filters you can write against it

```python
# fold gate (--mm 0, binder extracted): the prediction folds back onto its backbone
def apply_filter(df):
    return df[(df["usalign_fold_status"] == "OK") & (df["usalign_fold_RMSD"] < 2.0)]

# pose gate (--mm 1, whole complex, L_target = 193 here): the binder actually aligned
def apply_filter(df):
    ok = df["usalign_pose_status"] == "OK"
    return df[ok & ((df["usalign_pose_Lali"] - 193) > 40)]
```

The pose filter needs `L_target` as a literal, because the fixed target length is a property of the
campaign and not of the table. **Record the number you used in the campaign log** — a `Lali`
threshold without its `L_target` is unreadable six weeks later.

## What it is blind to

- **Everything chemical.** Geometry only: no charges, no solvation, no packing quality.
- **Whether the pose is physically possible in the real assembly.** A binder can align beautifully in
  the target frame and still clash with a symmetry-related protomer you trimmed away — that is
  `ringfit`'s question.
- **Which target residues the binder touches.** `Lali − L_target` says it stayed, not *where*.
  Hotspot recall on the predicted pose is a separate gate and neither substitutes for the other.
- **Whether the target itself landed.** Under a forced template that is the first gate of all, and a
  usalign comparison against the reference cannot distinguish a good binder from a collapsed target.
  Check target geometry before reading anything here.

## `_RMSD` IS OVER ALIGNED RESIDUES ONLY — never threshold on it without `_Lali`

**This is the single easiest way to build a filter that selects the worst designs in a pool,
and it looks completely reasonable while doing it.** `_RMSD` is computed over the `_Lali`
aligned pairs, *not* over the whole structure. So the fewer residues align, the smaller the
RMSD — a structure that diverged so badly it barely aligns reports a **beautiful** number.

*Measured, 272 binder complexes, whole-complex `--mm 1` prediction vs designed backbone
(193-residue templated target + a 63–90 residue binder):*

| pose RMSD bin | n | mean hotspot recall | median distance to nearest hotspot |
| --- | --- | --- | --- |
| **< 1.0 Å** | 26 | **0.144** | **21.17 Å** |
| 2.5–3.0 Å | 55 | 0.650 | 2.73 Å |

**Spearman(RMSD, hotspot recall) = +0.407 — positive.** The 26 "best" rows by RMSD were
binders sitting **21 Å off the epitope**: they failed to align at all, so the alignment was
target-only and the RMSD was measured over 193 perfectly-templated target residues.
**23 rows had `Lali - L_target == 0` — literally zero binder residues aligned — at a mean
RMSD of 0.51 Å.**

What to do instead:

- **Threshold on `_TM` , not `_RMSD`.** TM is length-normalised over the full structure, so a
  non-aligning region costs you. Same data: Spearman(TM, recall) = **+0.529**, and every row
  above TM 0.8 had non-zero recall.
- **Always report `_Lali` beside any RMSD**, and on a complex report `_Lali` minus the fixed
  chain's length — "how much of the thing I care about actually aligned". It correlated
  **+0.675** with recall, better than TM. Caveat: `_Lali` is aligned pairs summed over *all*
  chains, so that subtraction is a lower bound unless you confirm the fixed chain aligned fully.
- **When both structures have identical sequence and residue count** (e.g. the same binder in
  two conformations), USalign's alignment-based RMSD is the wrong tool entirely — do a forced
  1:1 Kabsch over every residue instead. Same trap, opposite direction: a *larger* real
  conformational change can report a *smaller* `_RMSD`.
