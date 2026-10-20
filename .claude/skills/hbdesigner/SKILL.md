---
name: hbdesigner
description: >-
  HBDesigner on Modal AND vib: designs buried hydrogen-bond networks onto an existing backbone, as monomer or one-/two-sided interface design, and writes `hbdesigner_fold_sequence` — the input binder's sequence with the network written in, which is what a predictor folds. Load before composing an hbdesigner run or reading its columns. Traps: a network must touch EVERY chain of the scaffold or the run returns nothing, the structure output is poly-glycine, zero output counts as success, fan-out is variable, the two backends run DIFFERENT checkouts, and `--guide-res` needs 2 distinct residues and is not reproducible.
---

# hbdesigner

**Custom tool** (lives in `tools/hbdesigner/`, not bundled with prosapia). Wraps
HBDesigner (MIT; PyRosetta is free for academic use). Runs on **Modal and vib**,
and the two backends run **different checkouts** — see the next section before
pooling their rows.

**`action: create`** — each input backbone yields up to `--top-k` ranked networks,
which are new entities, so it mints a **child table** with rows `<parent>_hb1 ..
<parent>_hb<k>` (the `_hb<i>` convention mirrors proteinmpnn's `_f<i>`). GPU by
default (L4 + 16 CPUs); packing is the CPU-parallel half.

## Two backends, two HBDesigners

| | Modal | vib |
| --- | --- | --- |
| checkout | RosettaCommons/HBDesigner, commit `ed65fa053786394a33efcedd5624c80ffbfef12b` | **Kuhlman-Lab** checkout, commit `b21ae46c` |
| how it is provided | baked into `modal_image.py` with its ~160 MB of weights — no Volume, nothing to install | the lab conda env `hbdesigner_gpu`, activation `$SAPIA_VIB_SHARED_ACTIVATION/hbdesigner.sh` (`SAPIA_ACTIVATE_HBDESIGNER`), Python 3.10 |
| gemmi / biotite / hydride | all present | **none of the three** |
| the worker's PDB reader | gemmi | **Biopython** (`hbdesigner_pdb_parser` says which) |

**The vib scoring path**, traced through the installed source (`cli()`
inference_hbdesigner.py:1085 → `HBDesRunner.run` 266 → `validate_inputs` sets
`opts.packer = "hbpacker"` **unconditionally** at 263-264, with no CLI flag for it
→ `pack_with_hbpacker` 644-720 → `Pool.map(minimize_and_score_network)` 704-714 →
`minimize_and_score_network` inference/packing.py:35-155):

- `score_network` (packing.py:122, defined 404-448) writes `HB_Score_full`,
  `HB_Score_hb`, `Avg_Burial`;
- `get_satisfaction` (packing.py:124, defined data/hbnet.py:67-198) writes
  `saturation`, `buried_heavy_unsats`, `buried_unsat_Hpol`;
- `filter_packs` (887-892) applies `--max_BUNs`, `--max_BUPHs`, `--min_sat`,
  `--max_hb_score`; `rank_and_save` (382-474) writes the stats CSV.

**biotite and hydride are not needed on vib**: nothing under `inference/` imports
either, `biotite_hbond_detect` (hbnet.py:281) has no caller anywhere, and
`run_hydride` is reachable only through `pack_with_rosetta(mode="hydride")`, whose
only callers are the training scripts. `hbnet.score_protein` and
`hbnet.rosetta_hbond_detect` are likewise **training-only**
(scripts/train_hbpacker.py:466, scripts/train_hbdesigner.py:804) and are **not**
the vib scoring path — an earlier version of this skill said they were; they are
not.

**The poly-glycine premise is established on the vib checkout only.** Everything
`fold_sequence` does rests on HBDesigner clearing the whole scaffold's sequence
except the anchors and writing `unk_to_gly=True`. If the Modal pin instead
preserves the input sequence, every Modal row would come back with **NA fold
columns and the glycine note** (the independent check would see non-glycine
residues everywhere and refuse) — a loud failure, not a wrong sequence, but a
total one. **The first Modal run of this revision must check that**: look at
`hbdesigner_fold_seq_note` on a handful of rows before trusting or debugging
anything else.

Every flag this tool emits exists in the vib version's `--help`, and the stats CSV
it writes carries every column the collector reads. What is **not** established is
whether the Modal pin scores through these same functions. **Treat rows from the
two backends as separate populations** — `hbdesigner_pdb_parser` is the in-table
marker — until a run has compared them on the same backbones.

The env is not ours to modify; nothing may be installed into it. The tool's worker
therefore imports only the standard library at module level (no gemmi, no
Biopython, no prosapia) and imports its PDB reader lazily.

**Submit time is a different environment and hard-depends on gemmi.** The
manifest builder now reads every input structure (to check the chain span and
that each `--anchor-res` / `--guide-res` residue exists), and it runs in the
prosapia venv, which has gemmi but **no Biopython** — so the Biopython fallback
is not available there and a gemmi failure is a failed submit. A `ResidueMap`
error on any one input — a duplicate chain+resnum, insertion codes folded onto
one number — **aborts the whole submit**, not just that row. That is deliberate
(the design would fail in the task anyway), but it means one malformed structure
in a table blocks the batch until you filter it out with `-f`.

## The chain-span rule — read this before any interface run

Upstream's `filter_packs` keeps a packed network only when
**`r["n_chains"] == total_chains`**: a network must touch **every chain of the
scaffold as HBDesigner sees it** (after `--sel-chains` extraction). Residues can
only land on chains that are designable (not in `--omit-chains`) or named by
`--anchor-res`. So a three-chain target with `--omit-chains B,C --anchor-res B30`
can never return anything — chain C is unreachable, every sampled network is
discarded, and the run ends with `Successes: 0` after the GPU has loaded the
models. **Nine runs were lost to this.** The builder now refuses it at submit,
naming the scaffold, the unreachable chains and the three ways out. Monomer and
two-sided runs (no `--omit-chains`) are never affected.

### The one-sided recipe that works

```bash
# 1. cut the input to the anchor's protomer plus the binder, ANCHOR CHAIN FIRST
sapia run chainsel <run_dir> -t table0 -i <structure column> --chains B,D
# 2. design: the network spans both chains, the anchor pins the target residue
sapia run hbdesigner <run_dir> -t table0 -i chainsel_path \
    --omit-chains B --anchor-res B138 --n-res 3 --guide-seq H,H,X --n-samples 200
```

Why each part:

- **`--chains B,D` with the anchor chain first.** `--sel-chains` is not an
  alternative: HBDesigner resolves the anchor to an absolute residue index on the
  **full** scaffold and applies it *after* `--sel_chains` has shrunk the scaffold,
  so unless the anchor is in the file's first chain the task dies with
  `IndexError: index 261 is out of bounds for axis 0 with size 231`. The builder
  refuses that combination. Cut the file instead.
- **`--anchor-res B138` on a polar donor/acceptor.** Hydrophobic anchors are
  rejected upstream; the alphabet is `R,N,D,Q,E,H,K,S,T,W,Y`.
- **`--n-res 3` = the anchor plus 2 designed residues.** The anchor counts
  towards `--n-res`.
- **`--guide-seq H,H,X` to force one histidine.** `guide_seq` is an **unordered
  multiset** over `--n-res` entries, scaled by the number of DESIGNED positions.
  With one anchor only two positions are designed, so `H,X,X` is a soft
  preference, while `H,H,X` reliably yields one His. Observed: `H,X,X` gave His in
  1 of 3 networks, `H,H,X` in 4 of 4.

### What it produced

On six parents, **single-anchor** runs gave His-on-binder networks for `Ser`,
`Tyr`, `Thr` and `Asn` anchors, with the ring nitrogen 2.5-3.0 A from the anchor's
donor in some (e.g. `B110Y:D17H:D19W`, ND1 to B110 TYR OH at 2.94 A;
`B138D:D19H:D59K`, ND1 to B138 ASP OD2 at 3.05 A). **Two-anchor** runs — one
anchor per protomer, the other way round the chain-span rule — gave either
nothing or no binder histidine: forcing a residue onto each protomer spends the
network on the target and leaves too little of it on the binder. Prefer the cut
over the second anchor.

### The limitation this recipe carries

Cutting to one protomer **hides the other protomer the binder also contacts**. A
step that REDESIGNS the binder in that context cannot see those contacts, and the
pilot evidence is that it loses them:

- **Observed**: 0 of 248 refolded redesigns reached ipSAE 0.7 — but **only 57 of
  the 248 rows carried an ipSAE at all**, so the denominator for the rate is 57,
  not 248, and the other 191 are "not measured", not "failed".
- **Observed, paired control** (`report_runH_matched_control.txt`, in
  `campaigns/.analysis/20261013_142959_pool_v3_stage/reports/`): Run A children
  of the four shared parents — `kecheng_bd_t016_ri005_c6_4`, `luca_f1_j0_n_535`,
  `luca_f2_j1_n_518`, `luca_f2_j23_n_544` — are **32 rows, 27 with an ipSAE
  value**, reaching ipSAE >= 0.5 on **12**, >= 0.6 on **10**, >= 0.7 on **8**.
  The network arm over the same parents is **176 rows, 44 valued**, at
  **4 / 1 / 0** for the same three thresholds. Quote the denominators with the
  rates: both arms are mostly unvalued.
- **INFERENCE, from four parents only**: the loss tracks how many binder residues
  touch *only* the omitted protomer — 7-9 such residues for the three parents
  that lost, 2 for the one that did not. Four parents is not a result; treat it
  as a hypothesis that explains the pairing, not as an established mechanism.

**This is why `fold_sequence` exists** — it keeps every other binder residue
exactly as the input had it and changes only the network positions, so the
contacts to the hidden protomer survive by construction rather than by luck.

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
| `--guide-res` | none | Residues whose Cb centroid a virtual guide atom sits at, so the network forms near them: `A3,A26`. **PDB numbering.** `{expr}` islands resolved per design, and the checks below run **after** that resolution. **At least 2 DISTINCT residues** (a centroid of one is undefined, and `A12,A12` has the same centroid as `A12`) — fewer is a submit-time error. Tokens are `<chain letter><number>` only: **no ranges** (`A12-15`), no insertion codes (`A12A`), no negative numbers (`A-5`); all three are submit-time errors naming the flag and the token. For a single residue use `--anchor-res`. |
| `--guide-radius` | `1e6` (off) | Hard Cb-distance cap from the guide atom. **Refused on any design whose `--guide-res` resolved to nothing** (there is no guide atom to measure from, so upstream would ignore it), and refused at `<= 0`. Checked per design, after `{expr}`. |
| `--anchor-res` | none | Residues **every** returned network must contain: `B5`. PDB numbering, `{expr}` resolved, same token format as `--guide-res` (no ranges). **One residue is fine here** — this is the single-residue flag. This is what makes one-sided interface design work. Upstream rejects hydrophobic anchors (`restype_hb_sc` is `R,N,D,Q,E,H,K,S,T,W,Y`). |
| `--guide-seq` | all `X` | Per-position AA conditioning, **exactly one entry per `--n-res`** (a submit-time error otherwise; upstream's default is `["X"] * n_res`, not the literal `X,X,X` its own `--help` claims): `S,N,T` / `X,T,X` / `S,N\|Q,T`. Every alternative must be `X` or one of `R,N,D,Q,E,H,K,S,T,W,Y` — checked at submit. **No spaces** (the token is word-split). SER/THR networks pack far more often — `S,X,X` is a cheap success-rate win. |
| `--omit-aa` | none | AAs the network may not use: `R,K` (`--omit_AA`). Must be a subset of `R,N,D,Q,E,H,K,S,T,W,Y` (the only residues HBDesigner ever places), and must not forbid a letter `--guide-seq` asks for — both checked at submit. |
| `--omit-chains` | none | Chains HBDesigner **sees but may not design into** (except via `--anchor-res`). Chain mini-language (`B`, `A:C`, `A,C`). The one-sided interface flag. The omitted chain comes back poly-glycine → grafted back by default. |
| `--sel-chains` | none | Chains to run on; the rest are **removed before** design and put back once design has finished, **in the input's own chain order** (measured) but **as poly-glycine on vib, not with their sequence** (trap 18). Different from `--omit-chains`: the model never sees them, so it cannot avoid clashing with them. **Refused together with `--anchor-res` unless the anchor's chain is the file's FIRST chain AND is itself selected** — only then does the anchor's absolute index survive the extraction. The refusal is deliberately conservative: `--sel-chains A,B --anchor-res B5` would in fact work (every chain before B is selected), but it is refused too, because the general rule is easy to get subtly wrong and the failure mode when it is wrong is a silently mis-pinned residue rather than an error. On vib, prefer cutting the input with `chainsel` and using `--omit-chains`. |
| `--symm-chains` | none | Symmetrize the network across chains: `A,B` or `A,B;C,D`. Experimental upstream. |
| `--symm-file` | none | Rosetta `.symm` file for strict symmetry. Requires `--symm-chains`; the path must be readable inside the task container (put it in the run_dir). |
| `--max-buns` | `0` | Max buried unsatisfied **heavy** atoms (`--max_BUNs`). Higher = more permissive. |
| `--max-buphs` | `5` | Max buried unsatisfied **polar hydrogens** (`--max_BUPHs`). Higher = more permissive. |
| `--min-sat` | `0.5` | Minimum saturation, range 0–2. **Higher = stricter** — the opposite direction from the two `--max-*`. |
| `--max-hb-energy` | `0.0` | Rosetta `score:hb_max_energy` used while scoring. |
| `--max-hb-score` | `0.0` | Max HB_Score a returned network may have. More negative = stricter. |
| `--n-workers` | CPU count | Packing workers. Defaults to `-c/--cpus-per-task`, else 16. **Must be ≥ 1** on both backends — see the traps. The vib version's `--seed` help text suggests `n_workers 0`; its own code asserts `>=1`. Ignore the help text. |
| `--cpu` | off | CPU inference (`*_cpu.yaml` configs) and `gpus_per_task = 0`. Much slower; debugging only. |
| `--graft-chains` | `auto` | `auto` = graft exactly the `--omit-chains` chains; `none` = never; or an explicit list. See below. |
| `--set TOKEN` | none | Forward a raw `run_hbdesigner` flag verbatim, repeatable (e.g. the two custom-checkpoint flags). No whitespace inside a token. |

Resources: one task = one backbone. Allow real time — `--n-samples 200` with
PyRosetta packing is minutes, not seconds (`RESOURCES` timeout is 4 h).

## Columns (leaf-prefixed `hbdesigner_`)

**Trust columns first — read these before any science column.**

| Trust column | Meaning |
| --- | --- |
| `resnum_shift_max` | Max \|mpnn position − PDB resnum\| over the network residues. **0 = no renumbering happened**, so the resnums in `network` *are* the mpnn positions. Any non-zero value: read `resnum_offset` before using `fix<i>`. |
| `resnum_offset` | Per network chain, `mpnn position − PDB resnum` when constant over that chain, else `<chain>:var`. All zeroes ⇒ no shift. |
| `pdb_parser` | Which reader built the residue mapping: `gemmi` (Modal, or anywhere gemmi is importable) or `biopython` (the fallback, which is what vib uses). The in-table marker of **which backend produced the row**. The two agree on the standard 20 amino acids and on every well-formed layout checked, but **not on every residue name**: `UNK` and the D-amino acids are counted by gemmi and skipped by Biopython, and `PCA` reads as `E` to gemmi and `Q` to Biopython. A file containing one of those gets a reader-dependent `n_res_total`. Gate any cross-backend comparison on this column. |
| `guide_atoms_removed` | HETATM lines removed from this rank's PDB because they matched HBDesigner's `--guide_res` marker atom — **resname `ORI` AND atom name `V1`** (element `V`), both required so an unrelated `ORI` ligand is not deleted. `0` = checked, none present (the normal case without `--guide-res`). `1` = the marker was taken out and the file as HBDesigner wrote it kept under `<design>/guide_raw/`. Upstream adds exactly one, so **a value above 1 means the output does not match the code this tool was built against** — it is reported, not capped. **`NA` = the strip pass did not run or failed, so the file may still carry the marker — NA is not 0 here.** Collect prints how many OK rows are NA. |
| `grafted` / `graft_identity` | See below; `graft_identity` below 1.0 means the structure is wrong, not merely imperfect. |

**The foldable output.**

| Column | Meaning |
| --- | --- |
| `fold_sequence` | **The sequence a predictor may fold — not `path`.** The STAGED INPUT's binder-chain sequence with the network's binder positions overwritten by the network residues: every residue the input already had is kept, only the network changes. One-letter, in the chain's own order. **NA** (never an empty string) when the binder chain cannot be identified — read `fold_seq_note`. Only as good as the staged input: if that was poly-glycine, so is this. Feed it on with `mkcomplex -i hbdesigner_fold_sequence`, then a predictor. |
| `n_fold_changes` | How many of the network's binder positions differ from the staged input's residue there. `0` is legitimate (the network re-chose what was already there). Not a quality score, and neither direction is "better" — it is how far the binder moved. |
| `fold_chain` | The chain `fold_sequence` belongs to: the single network chain **not** in `--omit-chains`. NA with the other two. On a monomer run the one network chain is trivially the binder, so a monomer DOES get a sequence. |
| `centers` | The designed **histidines on the binder**, in hisbond's centres format: `'<index0>:HIS'` joined by `;`, e.g. `'11:HIS;19:HIS'`. **`index0` is 0-based into `fold_sequence`** (hence into the binder chain in file order), **not** a residue number — `fold_sequence[index0] == 'H'` is checked before the column is written. The token is always a bare `HIS`: hbdesigner designs the residue, not its protonation state, so `HIS-S`/`HIS-P` would assert something it never measured. **NA, never `''`**, when the network put no histidine on the binder — hisbond *raises* on an empty centres cell, so a `""` here would look hisbond-ready and is not. |
| `n_centers` | How many such histidines. **`0` is a measurement** (this network designed none) and pairs with `centers` NA. Both are NA together only when `fold_sequence` is NA, i.e. when there is no binder chain to index into at all. Neither direction is "better": one designed His is the usual pH-switch target, two is not twice as good. |
| `fold_seq_note` | NA when the sequence was built, else why not. The commonest reason is a network on two chains with neither omitted, where which one is the binder is genuinely unknown. A row can be `OK` and still have no fold sequence — the block is non-fatal and never costs a row its science columns. The two checks that can refuse a sequence on a *well-formed* run are worth knowing: the output and the staged input must have the **same ordered residue numbers** for the binder chain, and **every non-network binder position must be glycine in the rank PDB** (HBDesigner writes everything it did not design as `G`). The second is independent of the string that was built, so it is the one that catches a stale rank PDB or a mismatched staged input — and when it fires **treat the whole row as suspect, not only the fold columns**: if the rank PDB is from an earlier run, then `network`, the scores, `fixed_positions` and `path` describe that run too. **The check only runs when a binder chain was identified.** On a two-sided design, or whenever the omit filter leaves zero or two or more chains, the fold block exits early and no glycine check happens — there the per-task clear (trap 20) is the only defence against a stale rank PDB. |

| Column | Meaning |
| --- | --- |
| `path` | The designed structure — the **grafted** file when chains were grafted, else the raw rank PDB. Poly-glycine outside the network, so **do not fold it**; fold `fold_sequence`. |
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
# only designs whose residue mapping is provably unshifted, whose graft worked,
# and whose file is known to be free of the --guide_res marker atom.
# graft_identity is NA when nothing was grafted (a monomer run, or no
# --omit-chains), so testing `== 1.0` alone selects ZERO rows there -- allow NA.
# guide_atoms_removed is NA when the guide-atom strip did not run, so this
# also excludes those rows.
def apply_filter(df):
    gi = df["hbdesigner_graft_identity"]
    return df[
        (df["hbdesigner_status"] == "OK")
        & (df["hbdesigner_resnum_shift_max"] == 0)
        & (gi.isna() | (gi == 1.0))
        & (df["hbdesigner_guide_atoms_removed"].notna())
    ]
```

```python
# <run_dir>/filters/best_network.py
# the best network per backbone only
def apply_filter(df):
    return df[df["hbdesigner_rank"] == 1]
```

```python
# <run_dir>/filters/foldable.py
# rows that actually carry a sequence a predictor can fold.
# .notna(), NOT != "" : the column is NA when no binder could be named, and an
# empty string is never written, so a "" test silently keeps nothing useful.
def apply_filter(df):
    return df[
        (df["hbdesigner_status"] == "OK")
        & (df["hbdesigner_fold_sequence"].notna())
    ]
```

**Under `-l <label>` every column above becomes `hbdesigner_<label>_*`** — so a
second variant in the same run_dir is filtered as
`df["hbdesigner_hisrun_fold_sequence"]`, not `df["hbdesigner_fold_sequence"]`.

```bash
sapia run proteinmpnn <run_dir> -t table1 -i hbdesigner_path \
    -f <run_dir>/filters/good_networks.py
```

## Two hand-offs: pick the right one

- **Designing a binder from scratch** (an rfd3 backbone with no sequence yet):
  hand off to **proteinmpnn** with the network fixed, below. Everything outside
  the network still has to be designed.
- **Adding a network to a binder that ALREADY has a sequence**: hand off
  `hbdesigner_fold_sequence` to **mkcomplex**, then a predictor —
  `sapia run mkcomplex <run_dir> -t <child table> -i hbdesigner_fold_sequence …`.
  Redesigning the whole binder here would throw away the sequence that already
  works (and, after a one-protomer cut, the contacts to the hidden protomer — see
  the recipe section). The graft column changes only the network positions.

## The hand-off to hisbond and propka (a designed-His network)

**THIS HAND-OFF IS CURRENTLY BROKEN AND NEEDS A HUMAN DECISION.** `hisbond` no
longer takes a centres COLUMN: its `--centers` wants comma-separated 1-based
binder positions, supplied per design through `{expr}` islands, and one island
resolves to one integer. `hbdesigner_centers` is a packed `'11:HIS;19:HIS'`
string (0-BASED), which cannot ride through an island, and hbdesigner publishes
no numeric `centre<k>_idx1` columns of its own. Until it does, either add them to
hbdesigner (the way `atomiumph` has them) or read the positions off
`hbdesigner_centers` by hand and pass literals. Remember the +1: hbdesigner's
indices are 0-based, `--centers` is 1-based.

```bash
# once hbdesigner publishes numeric positions, the call reads:
sapia run hisbond <run_dir> -t <child table> -i <predicted structure column> \
    --centers '{<a 1-based position column>},{<another>}' \
    --binder-chain D \
    --target-chains A,B,C \
    --expect-binder-len-col hbdesigner_fold_sequence
```

- **`--target-chains` is required and has no default** — for the triD complexes
  it is the three target copies, `A,B,C`. So is `--binder-chain`.
- **`-i` must be the predicted structure of THIS ROW's `fold_sequence`.** The
  centre indices are positions in that sequence; pointed at the parent's
  structure, or at a sibling's, they index the wrong molecule. **hisbond cannot
  always catch that**: when `n_centers` <= `n_fold_changes` every centre is a
  position the network actually changed, so a wrong `-i` usually lands on a
  non-His and errors — but when **`n_centers` > `n_fold_changes`** at least one
  centre was *already* a histidine in the input, and there hisbond would happily
  measure the parent's own His as though it were the designed centre and return
  a plausible number. Check that `-i` is the prediction made from this row.
- **`--binder-chain` is the chain in the PREDICTED structure, not
  `hbdesigner_fold_chain`.** The letter is set by **mkcomplex's segment order**:
  in the triD complexes the three target copies are prepended, so the binder
  comes out as `D`. `hbdesigner_fold_chain` names the chain in hbdesigner's own
  output and need not be the same letter. Note the length check only
  discriminates when the lengths differ — with a target protomer the same length
  as the binder it would pass on the wrong chain.
- **`centre<k>` columns are in ascending index order**, the same order as the
  `;`-joined `centers` string, so `centre1` is the lowest-index histidine.
- Then `propka --resnum-cols` off hisbond's `centre<k>_idx1` columns, as in the
  `ph-switch` skill.

The order for a grafted design is: `mkcomplex -i hbdesigner_fold_sequence` → a
predictor → relax → `hisbond` / `propka`.

## The hand-off to proteinmpnn

For a binder with no sequence yet. Keep the network fixed and design everything
else:

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
- **`--chains-to-design` must list the chains `hbdesigner_network_chains` names,
  in that order**, or the groups land on the wrong chains. For a chain you design
  but that has no network residues, add an empty group: `'{hbdesigner_fix1}/'`
  designs chains `A,B` with only A's position fixed.
- **But do not copy `network_chains` into `--chains-to-design` blindly.** With
  `--anchor-res B5 --omit-chains B` the anchor is a TARGET residue, so
  `network_chains` is `A,B` — and `--chains-to-design A,B` would hand the target
  chain to proteinmpnn to **redesign**, fixing only the one anchor position. That
  is a redesigned target, which is almost never what you want. Either design the
  binder chain alone (`--chains-to-design A`, and drop B's group from
  `--fixed-positions`) and accept that the anchor is simply not redesigned, or
  pass **every** position of chain B as fixed. Read `hbdesigner_network_chains`
  to find out whether the network crossed the interface; do not feed it straight
  through.
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
removes the unused chains before design and puts them back once design has
finished, and `graft_seq.py` copies **by position** between two files it asserts have the
same residue count — so a graft is only correct if the two still line up residue
for residue, which this tool does not establish. Pass `--graft-chains none`.
Two older claims here were wrong and are withdrawn: the unused chains do **not**
keep their own sequence (they come back poly-glycine, trap 18), and the chain
ORDER is **not** changed — the observed `--sel_chains D` output was chains
A,B,C,D, the input's own order. Because they come back poly-glycine and the graft
is refused, `--sel-chains` leaves you with no way to restore them through this
tool at all.
**A FAILED graft is now an error row, not a dead task.** `graft_seq.py` asserts
`p_target.n_res == p_ref.n_res`, so a stray atom aborts it; the task script runs
the loop under `set +e`, records the exit code, and the worker writes
`error: graft failed …`. The task still exits non-zero, so it is red in the
ledger AND present in the table. Before this it died under `set -e`, wrote no
TSV, and collect skipped the design silently (trap 10).

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
**If the fold makes two residues share a chain+resnum, that is an error row**, not
a best guess: chain+resnum is how HBDesigner names a network residue, so an
ambiguous one cannot be mapped to a ProteinMPNN position. Both readers refuse it
(gemmi path by an explicit duplicate-key check, Biopython by `PERMISSIVE=False`,
which reports it as `Blank altlocs in duplicate residue`). Same for an altloc
point mutation or a file numbered twice.

**9. The input is copied in as `<design>.pdb`.** Upstream names every output after
the input file's stem, so the task script stages the backbone under the design name
— that is what makes `<design>_HBDes_rank_<i>.pdb` predictable and collision-free.
CIF inputs are converted up front (`ensure_pdb`); HBDesigner reads PDB only.

**10. A design with no result TSV is skipped, not failed.** If a task was killed or
is still queued, collect prints `no result TSV … (skipping)` and mints nothing —
no evidence is not the same as failure. Re-collect after the task finishes.

**11. Upstream's `pyproject.toml` does not list everything the code imports, and
`--no-deps` believes it.** **This trap is about the MODAL pin only, and does not
describe vib.** The vib env has neither biotite nor hydride and does not need
them: nothing under `inference/` imports either, and `biotite_hbond_detect` has
no caller anywhere in that checkout — so the lazy-import crash described here
cannot happen on vib. Read the rest of this trap as a note on the Modal image.
At the pinned commit `[project].dependencies` names only
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
  `buried_heavy_unsats`, `buried_unsat_Hpol` — i.e. **on Modal** it *is* the source
  of this tool's `hb_score_full`, `hb_score_hb`, `avg_burial`, `saturation`,
  `buried_heavy_unsats`, `buried_unsat_hpol` columns. (On vib those same six
  columns come from `score_network` and `get_satisfaction` on the hbpacker path
  — see the backend section. Whether the Modal pin reaches those same functions
  is not established, which is why the two backends' numbers are not assumed
  comparable.) Because the imports are lazy,
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
handful — `return df[df.index == "design_3"]` (the frame is indexed by name; `df["name"]` raises KeyError) — which is how you buy a cheap probe
before committing a batch.

**13. `--guide-res` with ONE residue is refused at submit time, and that is the
point.** The guide atom sits at the **Cb centroid** of the guide residues, so one
residue has no centroid. Before this check a single residue reached a GPU task and
died there: signature `ValueError: Only 1 provided. You must provide at least 2
guide res for centroid calculation.`, job dead in ~13 s with the GPU already
allocated. The builder now refuses, all after `{expr}` resolution and per design,
and a refusal **aborts the whole submit** rather than skipping that row (a
silently dropped design reads as a filter problem three steps later): fewer than
**two distinct** residues — note upstream counts *tokens*, so `A12,A12` would
pass there; requiring two **distinct** residues is this tool's policy, on the
grounds that the centroid of one residue listed twice is that residue and the
repeat is a typo;
`--guide-radius` on a design whose `--guide-res` resolved to nothing; and a range
/ insertion code / negative number in either `--guide-res` or `--anchor-res` —
HBDesigner parses each token with `re.split('(\d+)')` and unpacks exactly two
pieces, so `A12-15`, `A12A` and `A-5` all fail inside the task; and **a residue
that is not in the input structure** — upstream's `Parsing error! 0 matches found
for res B999` otherwise fires after the GPU is allocated, and the numbering it
wants is the input PDB's own (an rfd3 backbone restarts at 1 per chain, so a
crystal number is usually wrong). The same submit
pass also checks `--guide-seq` has exactly `--n-res` entries drawn from
`R,N,D,Q,E,H,K,S,T,W,Y` or `X`, and that `--omit-aa` is a subset of that alphabet
and does not forbid a letter `--guide-seq` asks for. **If you want to steer a
network at ONE residue, that is `--anchor-res`, not `--guide-res`.**

**14. `--guide-res` makes HBDesigner write a marker atom into every output PDB.**
`add_guide_atom` appends a HETATM — atom `V1`, residue `ORI`, element `V` — in a
**new chain after the last one**, because `rank_and_save` writes with
`to_pdb(unk_to_gly=True, no_hetatm=False)`. It is a record of where the network
was steered, not part of the structure, and the input has none, so a residue count
taken over HETATM differs between output and input. The task therefore **strips it
before anything else touches the file** (notably before upstream's `graft_seq.py`,
which asserts `p_target.n_res == p_ref.n_res`) and keeps the file as HBDesigner
wrote it under `<design>/guide_raw/`. Read `hbdesigner_guide_atoms_removed`: `0`
normal, `1` stripped, **`NA` means the strip did not run and the file may still
carry it** — in which case the next tool to read `hbdesigner_path` sees a one-atom
chain of element V in a poly-glycine backbone.
Two details worth knowing: `guide_raw/` is **refreshed from the current rank PDB
on every strip pass**, so it always holds this run's output and a `--force` rerun
in the same leaf cannot resurrect the previous run's coordinates; and the strip
round-trips the text, which normalises line endings to `\n`, so for a CRLF input
the kept copy is content-identical but not byte-identical (a file with nothing to
remove is not rewritten at all).

**This pass runs on Modal too.** The edit that added it is **additive, not
invisible**: Modal tasks now also write `guide_strip.tsv` and (with
`--guide-res`) `guide_raw/`, and every Modal row gains `pdb_parser` and
`guide_atoms_removed`. The six metric columns and the structure files are
unchanged when there is no marker to remove.

**15. The same tool, two different programs.** Modal runs the RosettaCommons pin;
vib runs the Kuhlman-Lab checkout `b21ae46c` out of the lab conda env, where the
six metric columns come from `score_network` + `get_satisfaction` on the
hard-coded `hbpacker` path (see the backend section for the file:line trace).
Flags and stats-CSV columns line up between the two; **the numbers have not been
shown to**, and whether Modal reaches the same scoring functions is unverified.
`hbdesigner_pdb_parser` (`gemmi` vs `biopython`) is the in-table marker of which
backend wrote a row. Give a backend switch its own `-l` label and do not pool the
two until a run has compared them on the same backbones.

**16. `--guide-res` is not reproducible, and the guide atom is not in the
scoring.** The guide atom used during SAMPLING is the Cb centroid of the guide
residues plus Gaussian noise of **sigma 4.0 A, resampled on every run**
(`guide_atom_sigma = 4.0`), so the designable set moves run to run: **35** (probe2,
`ledger_brief_runH_hbd_probe2.md`), **29 and 6** (probe3,
`ledger_brief_runH_hbd_probe3.md`) and **5** (case E5,
`ledger_brief_runH_hbd_probe4.md`) designable positions over one backbone — the
briefs are under `campaigns/.analysis/20261013_142959_pool_v3_stage/reports/`.
Do not treat two
`--guide-res` runs as replicates, and do not conclude a backbone is hopeless from
one. The marker written into the output PDB is a *different* draw with
`sigma=0.01` (the near-exact centroid) and is added **only when the output PDBs
are written** — it is never present during packing or scoring, so it cannot
influence a number.

**17. The omitted chain is scored as a poly-glycine backbone.** Burial,
saturation and the unsatisfied counts over an `--omit-chains` target see backbone
only — the target's real side chains are not there, except at the `--anchor-res`
position, which keeps its identity. So `avg_burial` and `saturation` on a
one-sided interface run are **optimistic about the target side**: a network that
looks buried against a glycine wall may clash with the real side chains. Check
the packed interface with a tool that sees the real target (`pyrosetta`,
`ifacegeom`, `epitope`) after the graft.

**18. `--sel-chains` does NOT return the unused chains with their sequence** on
the vib checkout — they come back **poly-glycine**. Observed: a `--sel_chains D`
run produced chains A,B,C,D with non-glycine residues only at the three network
positions. Note what that same observation shows about ORDER: the output chains
are A,B,C,D, **the input's own order** — an earlier version of this skill said
`--sel-chains` appends the unused chains at the end and may reorder them, and
that is withdrawn. The sequence loss is the real problem, and `--graft-chains` is
refused alongside `--sel-chains`, so there is no way to restore them here.
Combined with the anchor `IndexError` above, there is little reason to use
`--sel-chains` on vib at all: cut the input with `chainsel` and use
`--omit-chains`.

**19. The two PDB readers do not agree on every residue name.** gemmi's
amino-acid table has 129 entries, Biopython's 1032. They agree on the standard 20
— the whole of this tool's domain, since HBDesigner writes `unk_to_gly=True` — and
produced identical residue maps on every well-formed file checked. Outside that:
`UNK`, `SEC`, `PYL`, `ASX`/`GLX` and all 18 D-amino acids are counted by gemmi and
**skipped** by Biopython; 942 names go the other way; and `PCA` is `E` to gemmi and
`Q` to Biopython. **On vib the consequence is an error row, not a quietly
reader-dependent `n_res_total`**: the counts do not cancel, because HBDesigner
rewrites the sequence — an input containing `UNK` comes back as glycine, which
Biopython *does* count, so the output has more residues than the input and the
row becomes `error: output has N residues, input has N-k`. The message names the
parser and tells you to `grep -c UNK` the staged input. The same file passes
under gemmi, which counts `UNK` on both sides. `PCA` is the one that changes a
letter without changing a count; it surfaces through the residue-identity check.
Gate cross-backend comparisons on `hbdesigner_pdb_parser`.

**20. A rerun clears the design's old output, but a re-COLLECT does not retract
old children.** The task now deletes `<design>.tsv` and `<design>/` before
running, because HBDesigner writes no CSV and no PDB when nothing survives — so
without the clear, a `--force` rerun that found nothing would have had the worker
read the PREVIOUS run's stats CSV and rank PDBs and emit full `OK` rows with the
old parameter set's metrics. That is fixed. What is NOT fixed, and cannot be:
**`<parent>_hb<k>` children minted by an earlier collect are not retracted by a
later one.** Collectors mint rows; they do not delete them. A rerun that yields
fewer networks leaves the surplus children from the previous collect in the
table, pointing at files that no longer exist. So: **a changed parameter set gets
a new `--table-label`**, not a `--force` into the same child table. Use `--force`
only to re-run rows under the same parameters.

## On-disk layout

```
<run_dir>/<child_table>/hbdesigner/
├── <design>.tsv                               one row per kept network (+ the _hb0 record)
└── <design>/
    ├── <design>.pdb                           staged input = the graft reference
    ├── <design>_HBDes_rank_<i>.pdb            outputs, guide atom stripped
    ├── <design>_HBDes_stats.csv               upstream's metrics table
    ├── guide_strip.tsv                        per rank PDB, guide atoms removed
    ├── guide_raw/<design>_HBDes_rank_<i>.pdb  only when a guide atom was removed
    └── grafted/<design>_HBDes_rank_<i>.pdb    only when grafting ran
```

`<design>.tsv` and the whole `<design>/` directory are **deleted at the start of
every task** (trap 20), so what you see always belongs to the last run of that
design. The child ROWS from an earlier collect are not deleted — see trap 20.

## Image notes (Modal only — vib uses the lab conda env, see the table above)

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
