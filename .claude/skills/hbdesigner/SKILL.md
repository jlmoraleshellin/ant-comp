---
name: hbdesigner
description: How to run the custom hbdesigner tool on Modal — RosettaCommons/HBDesigner, which designs buried hydrogen-bond networks onto an existing backbone (GNN + PyRosetta). Covers monomer and one-sided/two-sided interface design, every flag, the poly-glycine output convention, the zero-output success, the variable fan-out, the graft-back step for omitted chains, and the fixed_positions hand-off to proteinmpnn. Load before composing an hbdesigner run or reading its columns.
---

# hbdesigner

**Custom tool** (lives in `tools/hbdesigner/`, not bundled with prosapia). Wraps
[RosettaCommons/HBDesigner](https://github.com/RosettaCommons/HBDesigner) (MIT;
PyRosetta is free for academic use), pinned at commit
`ed65fa053786394a33efcedd5624c80ffbfef12b` (HEAD of `main`, 2026-08-31) and baked
into the Modal image together with its ~160 MB of weights — **no weights Volume and
nothing to install by hand**.

**`action: create`** — each input backbone yields up to `--top-k` ranked networks,
which are new entities, so it mints a **child table** with rows `<parent>_hb1 ..
<parent>_hb<k>` (the `_hb<i>` convention mirrors proteinmpnn's `_f<i>`). GPU by
default (L4 + 16 CPUs); packing is the CPU-parallel half.

## What it does, and what it is not

HBDesigner picks `--n-res` (2–6) positions on an **existing** backbone and assigns
polar residues there so they form one connected, **buried** hydrogen-bond network,
then packs and scores it with PyRosetta and keeps the best `--top-k`.

**It designs those positions and nothing else. Every other position in the output
structure is written as GLYCINE** (`to_pdb(unk_to_gly=True)`). That is a design
convention, not corruption — but it means the output is a *network stub on a
backbone*, not a foldable protein. Anything that computes on the sequence (a
predictor, a composition filter, a hydrophobicity score) will see poly-Gly and
produce a confident wrong answer. The intended chain is:

```
rfdiffusion3  ->  hbdesigner  ->  proteinmpnn (network fixed)  ->  boltz  ->  usalign/pyrosetta
```

The `hbdesigner_fixed_positions` / `hbdesigner_fix<i>` columns exist to make that
third step exact.

Its numbers are **meaningless** if: you read the output as a designed protein; the
backbone has no buried volume (nothing to bury a network in — you get
`error: no networks`, not a bad score); or you designed one side of an interface
and skipped the graft-back (the partner chain is then poly-glycine, so every
downstream interface metric is computed against a chain that is not there).

It does **not** design the rest of the sequence, relax or validate the fold,
predict a structure, or score an interface. Its only judgement is the filter set
(`--max-buns`, `--min-sat`, …), and those are yours to set.

## Invocation

```bash
# monomer: a 3-residue buried network on each rfd3 backbone, 5 ranked designs each
sapia run hbdesigner <run_dir> -t table0 --n-res 3 --n-samples 200 --top-k 5
# poll .exit as usual, then collect into the CHILD table the run reserved:
sapia collect hbdesigner <run_dir> -t table1

# one-sided interface design: the network must include target residue B5, chain B
# is not designable, and chain B's real sequence is grafted back automatically
sapia run hbdesigner <run_dir> -t table0 -i rfdiffusion3_path \
    --n-res 3 --n-samples 200 --anchor-res B5 --omit-chains B

# two-sided interface design: nothing special — a multi-chain input is enough
sapia run hbdesigner <run_dir> -t table0 --n-res 3 --n-samples 200
```

### Flags

| Flag | Default | Notes |
| --- | --- | --- |
| `-i/--input-column` | `rfdiffusion3_path` | The backbone to design onto. The default is the column rfd3 actually writes (**not** proteinmpnn's broken `rfdiffusion_path`), and a column the table does not have is a **loud submit-time error**, never a silent zero-task run. |
| `--design-model` | `design_020` | `design_020` (moderate noise, general purpose) or `design_002` (low noise, conservative). Both ship in the image. |
| `--n-res` | `2` | Network size in residues, 2–6. The only thing designed. |
| `--n-samples` | `100` | Networks sampled before packing/scoring. Only a fraction survive (upstream: ~10 of 200). Upstream's scaling advice: `n_res 2 -> 100`, `3 -> 200`, `4/5 -> 500`, `6 -> 1000`. |
| `--top-k` | `5` | Ranked designs **kept per backbone**. You can get fewer — see the traps. |
| `--t-range LO HI` | `0.1 1.0` | Sampling temperature range (`--T_range`). |
| `--seed` | none | Upstream's own caveat: runs are not bit-reproducible even with a seed (parallel packing + Rosetta). |
| `--min-burial` | `0.0` | Minimum sidechain-neighbour burial for a designable position. Raise to force the network into the core. |
| `--min-core-res` | `0` | Minimum core residues per network. |
| `--guide-res` | none | Residues whose Cb centroid a virtual guide atom sits at, so the network forms near them: `A3,A26`. **PDB numbering.** `{expr}` islands resolved per design. |
| `--guide-radius` | `1e6` (off) | Hard Cb-distance cap from the guide atom. |
| `--anchor-res` | none | Residues **every** returned network must contain: `B5`. PDB numbering, `{expr}` resolved. This is what makes one-sided interface design work. Must be fewer than `--n-res`; upstream rejects hydrophobic anchors. |
| `--guide-seq` | all `X` | Per-position AA conditioning, one entry per `--n-res`: `S,N,T` / `X,T,X` / `S,N\|Q,T`. **No spaces** (the token is word-split). SER/THR networks pack far more often — `S,X,X` is a cheap success-rate win. |
| `--omit-aa` | none | AAs the network may not use: `R,K` (`--omit_AA`). |
| `--omit-chains` | none | Chains HBDesigner **sees but may not design into** (except via `--anchor-res`). Chain mini-language (`B`, `A:C`, `A,C`). The one-sided interface flag. The omitted chain comes back poly-glycine → grafted back by default. |
| `--sel-chains` | none | Chains to run on; the rest are **removed before** design and concatenated back after, keeping their sequence. Different from `--omit-chains`: the model never sees them, so it cannot avoid clashing with them. |
| `--symm-chains` | none | Symmetrize the network across chains: `A,B` or `A,B;C,D`. Experimental upstream. |
| `--symm-file` | none | Rosetta `.symm` file for strict symmetry. Requires `--symm-chains`; the path must be readable inside the task container (put it in the run_dir). |
| `--max-buns` | `0` | Max buried unsatisfied **heavy** atoms (`--max_BUNs`). Higher = more permissive. |
| `--max-buphs` | `5` | Max buried unsatisfied **polar hydrogens** (`--max_BUPHs`). Higher = more permissive. |
| `--min-sat` | `0.5` | Minimum saturation, range 0–2. **Higher = stricter** — the opposite direction from the two `--max-*`. |
| `--max-hb-energy` | `0.0` | Rosetta `score:hb_max_energy` used while scoring. |
| `--max-hb-score` | `0.0` | Max HB_Score a returned network may have. More negative = stricter. |
| `--n-workers` | CPU count | Packing workers. Defaults to `-c/--cpus-per-task`, else 16. **Must be ≥ 1** — see the traps. |
| `--cpu` | off | CPU inference (`*_cpu.yaml` configs) and `gpus_per_task = 0`. Much slower; debugging only. |
| `--graft-chains` | `auto` | `auto` = graft exactly the `--omit-chains` chains; `none` = never; or an explicit list. See below. |
| `--set TOKEN` | none | Forward a raw `run_hbdesigner` flag verbatim, repeatable (e.g. the two custom-checkpoint flags). No whitespace inside a token. |

Resources: one task = one backbone. Allow real time — `--n-samples 200` with
PyRosetta packing is minutes, not seconds (`RESOURCES` timeout is 4 h).

## Columns (leaf-prefixed `hbdesigner_`)

| Column | Meaning |
| --- | --- |
| `path` | The designed structure — the **grafted** file when chains were grafted, else the raw rank PDB. Poly-glycine outside the network. |
| `rank` | 1 = best by HBDesigner's own order: fewest `buried_heavy_unsats`, then fewest `buried_unsat_hpol`, then highest `saturation`, then lowest `hb_score_full`. |
| `hb_score_full` | ref2015 energy of the network per network residue, relative to the same backbone as poly-Gly. **More negative is better.** |
| `hb_score_hb` | The same difference under an hbond-only score (fa_rep 0.55 + hbond_sc + hbond_bb_sc). More negative is better; this is the one that is about hydrogen bonds. |
| `avg_burial` | Mean sidechain-neighbour count over the network residues. Higher = more buried; a buried network is the whole point. |
| `saturation` | Fraction of the network's polar groups that are satisfied (0–2). Higher is better. |
| `buried_heavy_unsats` | Buried unsatisfied heavy atoms left in the network. Lower is better; `--max-buns` gates it. |
| `buried_unsat_hpol` | Buried unsatisfied polar hydrogens. Lower is better; `--max-buphs` gates it. |
| `network` | HBDesigner's own description: `<chain><resnum><aa>` joined by `:`, e.g. `A12S:A16T:B5N`. **PDB numbering.** Includes `--anchor-res` residues. |
| `network_seq` | The network's amino acids in `fixed_positions` order, e.g. `STN`. |
| `n_network_res` | How many residues the network has. Usually `--n-res`, but anchors and symmetrization change it — read it, do not assume. |
| `network_chains` | Chains the network touches, comma-joined, in PDB order. **This is what `--chains-to-design` must list** for `fixed_positions` to line up. |
| `fixed_positions` | The network in proteinmpnn's `--fixed-positions` syntax: `12,26/16` (`,` between positions, `/` between chains, one group per chain in `network_chains` order). **1-based within each chain**, mapped through the output PDB's residue order. |
| `fix1 … fixK` | The same positions, one integer per column, for the per-row `{expr}` hand-off. |
| `resnum_offset` | Per network chain, `mpnn position − PDB resnum` when constant over that chain, else `<chain>:var`. All zeroes ⇒ the resnums in `network` *are* the mpnn positions. |
| `resnum_shift_max` | Max \|mpnn position − PDB resnum\| over the network residues. **0 = no renumbering happened.** The filterable form of the trust check. |
| `grafted` | `yes` if the omitted chains were restored from the input PDB, else `no`. |
| `graft_identity` | Fraction of the grafted chains now matching the input. 1.0 = fully restored; NA when nothing was grafted. **Below 1.0 means the structure is wrong, not merely imperfect.** |
| `n_res_total` | Residues in the output structure. Checked equal to the input's; a mismatch is an error status. |
| `status` | `OK`, or `error: …`. Gate on it — always. |

### Filters you can now write

**`-f` takes a path to a module, never an inline expression** — see trap 12. Each
filter below is a small `.py` file written into the run_dir and passed by path.

```python
# <run_dir>/filters/good_networks.py
# networks with no buried unsats and a well-buried, well-satisfied core
def apply_filter(df):
    return df[
        (df["hbdesigner_status"] == "OK")
        & (df["hbdesigner_buried_heavy_unsats"] == 0)
        & (df["hbdesigner_saturation"] >= 0.8)
        & (df["hbdesigner_avg_burial"] >= 5.0)
    ]
```

```python
# <run_dir>/filters/trusted_mapping.py
# only designs whose residue mapping is provably unshifted, and whose graft worked
def apply_filter(df):
    return df[
        (df["hbdesigner_resnum_shift_max"] == 0)
        & (df["hbdesigner_graft_identity"] == 1.0)
    ]
```

```python
# <run_dir>/filters/best_network.py
# the best network per backbone only
def apply_filter(df):
    return df[df["hbdesigner_rank"] == 1]
```

```bash
sapia run proteinmpnn <run_dir> -t table1 -i hbdesigner_path \
    -f <run_dir>/filters/good_networks.py
```

## The hand-off to proteinmpnn

The point of the tool. Keep the network fixed and design everything else:

```bash
sapia run proteinmpnn <run_dir> -t table1 -i hbdesigner_path \
    --chains-to-design A \
    --fixed-positions '{hbdesigner_fix1},{hbdesigner_fix2},{hbdesigner_fix3}/'
```

Two things make this work, and one limits it:

- **`--fixed-positions` is 1-based within each parsed chain**, not PDB numbering.
  The collector never copies HBDesigner's resnums; it maps them through the output
  PDB's own residue order. `resnum_offset` / `resnum_shift_max` report whether a
  shift existed. (rfd3 output is numbered 1..L per chain, so the shift is usually
  0; a trimmed crystal target is where it bites.)
- **`--chains-to-design` must equal `hbdesigner_network_chains`**, in that order,
  or the groups land on the wrong chains. For a chain you design but that has no
  network residues, add an empty group: `'{hbdesigner_fix1}/'` designs chains `A,B`
  with only A's position fixed.
- **`fixed_positions` is a string, so it cannot be templated per row.** prosapia's
  `{expr}` islands resolve **integers only**. Use the `fix<i>` columns for the
  per-row route; use the `fixed_positions` string when you are reading the table or
  running one design. A `{hbdesigner_fix3}` that is NA on some row raises a
  per-design error at submit time (by design — it is loud): filter to a uniform
  network size first, with a module whose `apply_filter` returns
  `df[df["hbdesigner_n_network_res"] == 3]`, passed as `-f <path to that module>`.

Upstream does this with LigandMPNN and "fix everything that is not glycine"; the
`fix<i>` columns are the same idea, computed from the reported network rather than
from the poly-Gly sequence, so anchors and grafted chains do not confuse it.

## Traps

**1. A run that finds nothing still exits 0.** `rank_and_save` returns early —
writing **no PDB and no CSV** — when nothing passes scoring or symmetrization
(`"No valid networks passed symmetrization. Ending run."`), and the process exits
0. This tool records that as a child row **`<parent>_hb0`** with an `error: no
networks …` status, an empty path and NA metrics. **`_hb0` is a failure record,
not a design.** Never read it as one; always gate on `hbdesigner_status == 'OK'`.
If you see many of them: raise `--n-samples`, lower `--n-res`, relax `--min-sat` /
`--max-buphs`, or try `--guide-seq S,X,X`.

**2. The fan-out is not `--top-k` per parent.** Upstream keeps
`min(top_k, n_surviving)`. A parent may contribute 5 rows, 1 row, or only the
`_hb0` record. Never compute an expected row count; read the table.

**3. Everything outside the network is GLYCINE.** By design
(`to_pdb(unk_to_gly=True)`). Do not fold `hbdesigner_path` with boltz, do not score
it with pyrosetta, do not read a sequence off it. Run proteinmpnn first.

**4. One-sided interface design needs the graft-back, and it is automatic.** With
`--omit-chains B`, chain B is seen but not designed — and it still comes back
poly-glycine, because the scaffold's sequence is cleared before design. The task
script runs upstream's `graft_seq.py` against **this tool's own input PDB** for
exactly the chains `--omit-chains` named, writing into `<design>/grafted/`, and
`hbdesigner_path` then points at the grafted file. It only overwrites positions
that came back as glycine, so designed and anchor residues survive. Check
`hbdesigner_grafted == 'yes'` and `hbdesigner_graft_identity == 1.0`; without the
graft the downstream structure is simply wrong.
**`--graft-chains` is refused together with `--sel-chains`**: `--sel-chains`
appends the unused chains at the **end** of the output, so output and input no
longer share a residue order, and `graft_seq.py` copies **by position** — it would
write the wrong residues. `--sel-chains` needs no graft anyway (it returns those
chains with their own sequence), but note it can **reorder the chains** relative to
your input.

**5. `--n-workers 0` is not "no parallelism" — it crashes.** Upstream builds its
packing DataLoader with `persistent_workers=True`, which raises on
`num_workers=0`, and `Pool(0)` raises too. The builder refuses anything below 1.
Use `--n-workers 1`.

**6. `wandb` is a hard import dependency.** `inference_hbdesigner` imports the
training module, which imports `wandb` at module level. The image and the task
script both set `WANDB_MODE=disabled` / `WANDB_SILENT=true` so a container never
tries to phone home. Do not remove either.

**7. HBDesigner speaks PDB numbering; proteinmpnn does not.** `--guide-res`,
`--anchor-res` and the reported `network` are all PDB chain+resnum.
`--fixed-positions` is 1..L per chain. See the hand-off section — this is the most
expensive mistake available here, because a wrong mapping produces a plausible
number instead of an error.

**8. Insertion codes are folded into the residue number** when HBDesigner writes a
PDB (`res.id[1] + insertion_code_offset`), so output numbering can differ from a
crystal input's. The mapping is computed from the output file, so
`fixed_positions` stays right; `resnum_offset` will say `var` for that chain.

**9. The input is copied in as `<design>.pdb`.** Upstream names every output after
the input file's stem, so the task script stages the backbone under the design name
— that is what makes `<design>_HBDes_rank_<i>.pdb` predictable and collision-free.
CIF inputs are converted up front (`ensure_pdb`); HBDesigner reads PDB only.

**10. A design with no result TSV is skipped, not failed.** If a task was killed or
is still queued, collect prints `no result TSV … (skipping)` and mints nothing —
no evidence is not the same as failure. Re-collect after the task finishes.

**11. Upstream's `pyproject.toml` does not list everything the code imports, and
`--no-deps` believes it.** At the pinned commit `[project].dependencies` names only
pyrosetta, numpy, omegaconf, biopython, wandb, networkx, pandas and pebble. An AST
scan of all 29 files under `hbdesigner/` plus a transitive import trace from
`inference_hbdesigner` found four more, now installed explicitly in the image:

- **`scipy`** — module-level in `data/hbnet.py`, `data/protein.py`,
  `inference/protein_ops.py`. Signature:
  `ModuleNotFoundError: No module named 'scipy'` from `hbnet.py` during the
  `run_hbdesigner --help` build step. This killed the first ever build.
- **`git` (GitPython)** — module-level in `train/trainer.py`. `inference_hbdesigner`
  transitively imports the *training* module (the same quirk as trap 6), so an
  inference-only run needs **both `wandb` and GitPython**. This would have been the
  next build failure.
- **`biotite` and `hydride`** — **lazy** imports, inside `biotite_hbond_detect()`
  and `run_hydride()` in `data/hbnet.py`, and in upstream's `[dependency-groups].dev`
  rather than `[project].dependencies`. `biotite_hbond_detect()` is the function
  that returns `HBScore_Full`, `HBScore_HBond`, `Avg_Sc_Neighbors`, `saturation`,
  `buried_heavy_unsats`, `buried_unsat_Hpol` — i.e. it *is* the source of this
  tool's `hb_score_full`, `hb_score_hb`, `avg_burial`, `saturation`,
  `buried_heavy_unsats`, `buried_unsat_hpol` columns. Because the imports are lazy,
  **a missing biotite/hydride passes the whole image build and then crashes a GPU
  task mid-scoring** — the expensive failure mode. The fourth smoke test exists
  solely to convert that into a build-time error.

(`tqdm` is also undeclared, but only in `scripts/merge_networks.py`, off the
inference path; it is installed anyway.) The lesson generalises: for this repo,
`--no-deps` plus upstream metadata is not a dependency list — a source scan is.

**12. `-f/--filter` is a module path, not an expression.** prosapia declares it as
`type=Path` and loads it with `importlib.util.spec_from_file_location`; the module
must define `apply_filter(df) -> df`, or it raises
`AttributeError: … does not define an 'apply_filter' function`. Passing a
pandas-style string dies at **submit** time with
`ImportError: Could not load module from <your expression>` — measured, and it
costs a submit. Write the filter into the run_dir (on Modal, through the
workstation with a `--cmd` heredoc — a local path does not exist in the container)
and pass that path. Because the module is applied to the whole DataFrame before the
manifest is built, it is also the supported way to run on **one named design** or a
handful — `return df[df["name"] == "design_3"]` — which is how you buy a cheap probe
before committing a batch.

## On-disk layout

```
<run_dir>/<child_table>/hbdesigner/
├── <design>.tsv                               one row per kept network (+ the _hb0 record)
└── <design>/
    ├── <design>.pdb                           staged input = the graft reference
    ├── <design>_HBDes_rank_<i>.pdb            raw outputs
    ├── <design>_HBDes_stats.csv               upstream's metrics table
    └── grafted/<design>_HBDes_rank_<i>.pdb    only when grafting ran
```

## Image notes

Python **3.10 strictly** (upstream pins `>=3.10,<3.11`), so it cannot share the
built-in pyrosetta tool's 3.12 image. Torch 2.6.0+cu124 with the matching prebuilt
`torch-cluster` / `torch-scatter` wheels (upstream's `gpu-cu124` extra) — those are
binary wheels **only** for that exact torch, so nothing is allowed to re-resolve
them. PyRosetta comes from upstream's pinned public wheel URL
(`PyRosetta4.Release.python310…2024.39`), not `pyrosetta-installer`. The repo is
installed **editable** on purpose: `inference_hbdesigner` resolves its weights as
`Path(__file__).parents[2]/model_weights/*.pt`, which only exists when the package
runs from the checkout, and setuptools' package discovery would otherwise drop
`hbdesigner/scripts/` (no `__init__.py`) and with it `graft_seq.py`. First build is
long (torch + PyRosetta); later runs are seconds.

**Upstream under-declares its dependencies**, so the image installs the package
with `--no-deps` and curates the dependency set by hand (see trap 11). On top of
upstream's declared list it adds `scipy>=1.10,<1.16` (the `<1.16` is load-bearing:
scipy 1.16 dropped Python 3.10), `GitPython`, `biotite`, `hydride` and `tqdm`, and
re-asserts `numpy==1.26.4` in that same `pip_install` call so the resolver cannot
drift numpy while solving scipy/biotite. A **fourth smoke test**,
`python -c 'import scipy, git, biotite.structure, hydride, tqdm'`, exists because
`run_hbdesigner --help` only exercises module-level imports — it cannot catch a
lazily imported runtime dependency. **When the pinned commit is bumped, re-run the
import scan**; upstream's metadata will not tell you what is missing.
