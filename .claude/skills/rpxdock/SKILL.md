---
name: rpxdock
description: How to run the custom rpxdock tool on Modal — RPXdock rigid-body docking of a scaffold into a ONE-COMPONENT symmetric architecture, scored with precomputed residue-pair motif tables. SCOPE: only CYCLIC (C2–C17, CxSTACK) docking is supported and verified; the dihedrals Dx_y and the one-component cages are deliberately OUT OF SCOPE and unverified — read the scope section before using them. Covers the premise (origin-centred single-chain monomer input), every flag, the standing campaign decision to use the afilmv_ehl hscore alias and when to revisit it, the hscore Volume and its three aliases, the trust metrics that catch a silently meaningless score, the numpy pin and what it costs, the submit-time input guard (chain count, backbone, models, numbering) that turns a measured >10-minute willutil hang into an instant error, and the failure signatures found while building it. Load before composing an rpxdock run or reading its columns.
---

# rpxdock — symmetric rigid-body docking with RPX motif scores

**What it answers.** *Given this single-chain scaffold, is there a rigid-body
placement that makes a good **C2** (or C3, … C17) interface, and what does that
assembly look like?* RPXdock enumerates the placements compatible with a
point-group symmetry
and scores the interface each one creates by looking its residue pairs up in
precomputed motif tables ("hscore"). Paper: Sheffler et al., *PLoS Comput Biol*
2022, [10.1371/journal.pcbi.1010680](https://doi.org/10.1371/journal.pcbi.1010680).

**`action: create`.** Each kept pose is a new entity: one parent scaffold row fans
out to up to `--nout-top` child rows in a child table, named `<parent>_d1`,
`<parent>_d2`, … (rank 1 = best scoring).

**CPU-only.** The manifest builder forces `gpus_per_task = 0`, so you do **not**
need `-g 0`.

---

## ⚠ SCOPE: only CYCLIC (`Cx`) docking is supported and verified

**Deliberate decision, 2026-10-02, by the campaign owner.** Everything
non-cyclic — the dihedrals `Dx_y` and the one-component cages
`T2 T3 O2 O3 O4 I2 I3 I5` — is **out of scope, unverified, and must not be used
without doing the work listed below first.** The code paths still exist and will
still run; nothing stops you. That is exactly the risk: they will return
plausible-looking numbers.

**Use only `C2`…`C17` (and `CxSTACK`), with a SINGLE-CHAIN monomer input.**

### What is actually verified (cyclic only)

| Claim | Evidence |
| --- | --- |
| `C2` docks real scaffolds end to end | 100-res rfd3 monomers → 6 rows, 58–65 s; 119-res fixture → 167 poses, 74 s |
| row count matches disk | 6 rows ↔ 6 `.pdb`; 5 ↔ 5 on the fixture |
| `n_res` matches the parent's length | 100 = 100 on every row |
| `score = 1.0·rpx + 0.01·ncontact` | reproduces to 4 dp on all 7 rows across two independent runs |
| bad input errors, never fabricates | `error: no CA atoms found in …`, 0 child rows minted |
| the submit-time guard | 140 local assertions, plus a real-data check on cyclic |

### What must be done BEFORE enabling dihedral or cage docking

Do not simply pass `--architecture D3_2` and trust the output. Each of these is a
real, open gap:

1. **The chain-count guard is laxer there, and that is untested.** Cyclic refuses
   `>1` chain; dihedral/cage refuse only `>2`, because a `Dx_2` scaffold is
   legitimately a dimer on a 2-fold. **But a 2-chain input is exactly what caused
   the measured >10-minute willutil hang** (see the trap below), and causation was
   never proven. **So the dihedral/cage path may still hang on a 2-chain input and
   the guard will not stop it.** The test that would settle this — the same
   154-residue dimer into `D3_2`, watching for the `pdbfile.py:295` stall — was
   deliberately NOT run. Run it first.
2. **Cages and dihedrals have never been verified end to end**, in any run. See
   *Cages and dihedrals are not verified end to end* below for the numpy-2
   sampler bug, the `bad strides` empty-beam signature, and why upstream's own
   `test_onecomp.py` does not cover this configuration either.
3. **The input contract is different and nothing checks it.** Cages and dihedrals
   need the **asymmetric unit of a cyclic oligomer, pre-aligned with its symmetry
   axis on Z** — not the whole oligomer, which RPXdock regenerates. No code here
   can verify the axis is where you think it is. Get it wrong and you get a
   well-formed dock of a meaningless building block.
4. **`disp` only applies there.** It is NA for cyclic, so no reader has yet
   exercised that column on real data.

When you do enable them, update this section rather than deleting it — record
what you verified and how.

---

## Premise — the conditions under which the numbers are meaningless

Read this before trusting any score. RPXdock will happily dock the wrong thing and
report plausible numbers.

1. **One component only.** This tool drives RPXdock's `--inputs1` and nothing else.
   Supported: `C2`…`C17` (and `CxSTACK`), the one-component cages
   `T2 T3 O2 O3 O4 I2 I3 I5`, and the dihedrals `Dx_y` (x ∈ 2 3 4 5 6 8, y = 2 or
   x). Everything else — `T32`, `O43`, `I32`, `ASYM`, `AXLE_*`, `PLUG_*`, and the
   layers — needs a second/third input list that is not wired, and is **refused by
   name at submit time**, before any task is queued.
2. **The input must be centred on the origin.** RPXdock does not centre inputs and
   its samplers assume an origin-centred body. The collected `input_com_dist` says
   how far off it was; `--recenter-input` fixes it (see the flag's caveat).
3. **For `Dx_y` and the cages, the input must be the ASYMMETRIC UNIT of a cyclic
   oligomer, pre-aligned with its own symmetry axis on Z** — one protomer, not the
   whole oligomer, which RPXdock regenerates itself. Nothing in the tool can check
   this. Get it wrong and you get a well-formed dock of a meaningless building
   block.
4. **Every chain of the input is collapsed into one and renumbered from 1.**
   `Body.set_pose_info` sets `chain` to a constant and `resno` to `arange(n)`. So a
   multi-chain input is docked as one long concatenated chain, and `--allowed-residues`
   speaks in **positions in the file (1..N)**, never in author residue numbers.
   `n_chains_in` is collected so you can filter on it.
   **Since the multi-chain hang this is enforced at submit time:** a cyclic dock
   refuses more than **1** chain, a cage or dihedral more than **2** — see the trap
   below. For `Cx` the input you want is one protomer; use `chainsel` to get it.
5. **Scores are motif-table lookups, not energies.** Comparable within one batch at
   one `--hscore-files` setting; not across settings, and not against `pyrosetta`'s
   `if_dG` or anything physical. The `hscore` column records which table produced
   each number precisely so two batches cannot be compared by accident. **Which
   alias to pass is a campaign decision — see the next section; the default is
   helix-only and wrong for beta-containing designs.**

**What it does NOT do:** no sequence design (the dock is your scaffold's backbone,
repeated by symmetry — send it to `proteinmpnn`/`atomium` next), no relaxation, no
interface-area or packing metrics (that is `cms` / `pyrosetta` on the dumped
structure), no multi-component assemblies, and no check that your input's symmetry
axis is where you think it is.

---

## Choosing the hscore alias — a standing campaign decision, READ FIRST

> ### For this campaign the alias is **`afilmv_ehl`**, always.
>
> **Every run must pass `--hscore-files afilmv_ehl` explicitly.** The tool's
> built-in default is `ilv_h`, whose tables are built from **helix pairs** and
> applied SS-**in**dependently — the wrong table set for these designs, whose
> targets may contain **beta-sheets**. `afilmv_ehl` is the **SS-dependent** set
> covering all secondary-structure types.
>
> *Measured, and the reason this is a preference rather than a correctness bug:*
> `ilv_h` does **not** return ~0 on a β-rich scaffold — it scored `rpx` 78.90 on a
> 0 %-helix / 65 %-sheet design, against `afilmv_ehl`'s 107.74 on the same one.
> So a forgotten flag gives you plausible, wrongly-scaled numbers rather than an
> obvious failure. See *`frac_helix` is NOT a gate* further down.
>
> **A run that forgot the flag silently used `ilv_h`, and its scores should be
> DISCARDED — not re-interpreted, not rescaled.** It will not error, it will not
> warn, and it will return plausible numbers; the only thing that says so
> afterwards is the collected `hscore` column.

### The revisit condition — this is a choice, not a law

The alias follows the **secondary structure of the designs**, not habit, and not
this page. Revisit it when the designs change:

| If the designs are | Use | Why |
| --- | --- | --- |
| beta-containing, or mixed, or unknown | **`afilmv_ehl`** | SS-dependent, all SS types. The current campaign. |
| **all-helical** | `ilv_h` | Cheaper and much faster: **365 MB vs 5.70 GB**, and a measured **33–49 s** table load instead of a far larger one. Nothing is gained by paying for beta tables that no design uses. |
| anything | **not `ailv_h`** | see below |

A future all-helical campaign should drop back to `ilv_h` deliberately and say so
here. Do not inherit `afilmv_ehl` just because this section exists.

### `ailv_h` — do not use, not tested

* **~1.4 GB**, holding BOTH the `.txz` tarballs (~475 MB) and the `.txz.pickle`
  sidecars (~923 MB).
* **The documented rename-to-`.bak` workaround has NOT been applied** on the
  Volume, so `get_hscore_file_names` will load the **pickles**, not the tarballs
  (pickles win — see the Volume section).
* **The pickles are untested** and may fail to unpickle under this image's
  Python 3.12 / numpy 1.26: they were generated elsewhere, under an older stack.
* **The user has explicitly deprioritised testing it.** Anyone reaching for
  `ailv_h` must test the sidecars first — a single-design smoke run — before
  trusting or fanning out.

### ⚠ Operational consequences of `afilmv_ehl`

**Memory: `--mem 64G` is verified sufficient**, for one scaffold, both with and
without the pickles. The image default is `cpu 4 / memory 16G` and the memory must
hold the **decompressed** tables, so 16G is not enough — always pass `--mem`.
Peak RSS was never measured, so whether something lower than 64G would do is
unknown; 96G was used for the one-off pickle generation and also sufficed.

**Table-load time dominates — but the pickles cut it ~3×. All measured, same
scaffold (`denovo_diff_0`, 100 res, C2), same flags, `--mem 64G`:**

| | `ilv_h` (tarball) | `afilmv_ehl` tarball | `afilmv_ehl` **pickle** |
| --- | --- | --- | --- |
| `hscore_seconds` (load) | 30.5 s | **617.6 s** | **199.2 s** |
| total `seconds` | 64.0 s | 655.1 s | **236.9 s** |
| search (total − load) | 33.5 s | 37.4 s | 37.7 s |
| load as % of runtime | 48 % | **94 %** | 84 % |

* **The pickles are 3.10× faster to load and 2.76× faster end to end.** The search
  itself is untouched (37.4 → 37.7 s) — every saving is in the load, as expected.
* Loading scaled **20.2×** between `ilv_h` and `afilmv_ehl` tarballs, slightly
  *worse* than their 15.6× size ratio.
* The cost is still paid by **every container** (one scaffold per task; nothing is
  amortised). For a 20-design batch that is **1.1 h of table loading with pickles,
  against 3.4 h without.**
* **Correctness was verified, not assumed:** the pickle run reproduced the tarball
  run **bit-identically** — `score` 109.62184143066406, `rpx` 107.7418441772461,
  `ncontact` 188, `n_docks` 138. rpxdock is deterministic, so any drift there would
  have meant the pickles were not equivalent.

**Re-baseline every timing rule.** Any "a dock should finish in under N minutes"
guidance — including the 10-minute figure used to recognise the multi-chain hang
below — was set against **`ilv_h`**. Under `afilmv_ehl` a **healthy** task spends
**~200 s (pickles) or ~620 s (tarballs) loading before docking starts**, which looks
exactly like a hang. Distinguish them with the evidence table in that trap: an
`.err` showing the `using hscore …` line and nothing further, with an empty design
dir, is **normal** for the first few minutes. A fresh alias or a changed `--mem`
deserves a single-design probe before any fan-out.

### Scores from two aliases must never be compared

Cross-reference of the rule in Premise 5: **a motif score is a table lookup, so it
is comparable only within one `--hscore-files` setting.** The collected `hscore`
column carries the alias on **every row** precisely so two batches cannot be mixed
by accident. **A campaign that switches alias midway has two incomparable score
scales in its lineage** — rank within an alias, never across, and if you must
switch, re-dock the earlier scaffolds under the new alias rather than rescaling.

---

## Invocation

```bash
# C2 from a table of rfd3 backbones, 5 docks kept per scaffold.
# --hscore-files afilmv_ehl is MANDATORY for this campaign (see the section above);
# the built-in default ilv_h is helix-only and wrong for beta-containing designs.
sapia run rpxdock <run_dir> -t table0 -i rfdiffusion3_path \
    --architecture C2 --nout-top 5 \
    --hscore-files afilmv_ehl --mem 64G
sapia collect rpxdock <run_dir> -t table1      # the table the run reserved
```

**`--mem 64G` above is a placeholder, not a measured value** — size it from the
single-design timing probe described in the alias section.

`-i/--input-column` is **effectively required**: `default_input_column` is the
`"not applicable"` sentinel (as in `chainsel`/`cms`/`usalign`), because the scaffold
may live in `rfdiffusion3_path`, `boltz_path` or `chainsel_path` and silently
picking one would make the others a no-op. Omitting it is a named error at submit
time.

Use `-l/--dir-label` to run two architectures off one parent table without the
leaves colliding, and `--table-label` to keep the two child tables apart.

### The column prefix is the DIR-LABEL, not the tool name

**Measured, and it silently breaks every filter.** Columns are leaf-prefixed, and
the leaf is `<tool>[_<dir_label>]`, not `<tool>`. So:

| Run | Output dir | Columns |
| --- | --- | --- |
| `sapia run rpxdock …` | `<table>/rpxdock/` | `rpxdock_score`, `rpxdock_status`, … |
| `sapia run rpxdock … -l verification` | `<table>/rpxdock_verification/` | `rpxdock_verification_score`, `rpxdock_verification_status`, … |
| `sapia run rpxdock … -l c3` | `<table>/rpxdock_c3/` | `rpxdock_c3_score`, … |

A filter written against `rpxdock_score` on a labelled run **matches nothing and
raises nothing** — `-f` just returns an empty frame and the next step submits zero
tasks. Read the child table's header before writing a filter, and **collect with
the same `-l` you ran with** or the collector writes a second, unlabelled leaf.

Everything below is written with the bare `rpxdock_` prefix. On a labelled run,
substitute `rpxdock_<label>_` throughout.

---

## Flags

| Flag | Default | What it does |
| --- | --- | --- |
| `--architecture` | **required** | The symmetry, in RPXdock's spelling. One-component only (see Premise 1). |
| `--nout-top` | `10` | Docks kept per scaffold = child rows per parent. Fewer may come back: redundancy filtering happens first. There is deliberately **no `--nout-each`** — one scaffold per task means RPXdock's "top overall" and "top per dock" are the same set. |
| `--hscore-files` | `ilv_h` | Motif-table alias (see the Volume below). **The default is helix-oriented; this campaign requires `afilmv_ehl` — pass it on every run.** Changing it changes the score scale. |
| `--hscore-data-dir` | `/rpxdock_files` | Where the alias is looked up. This is RPXdock's own default **and** the Volume mount point, so it needs no flag on Modal. |
| `--docking-method` | `hier` | `hier` = hierarchical coarse-to-fine (what the paper uses, what `--beam-size` tunes); `grid` = flat enumeration. |
| `--max-trim` | `0` | Residues RPXdock may trim from a terminus to relieve a clash. `0` disables trimming and is markedly faster. When >0 the dumped structure may be **shorter** than the scaffold — check `reslb`/`resub`. |
| `--beam-size` | `100000` | Samples carried into each hierarchical stage after the first. The main runtime knob. |
| `--recenter-input` | off | Move the backbone CoM to the origin before docking. Fixes Premise 2, but the transforms RPXdock reports are then relative to the **recentred** pose — pre-centre the input instead if you plan to use them. |
| `--score-only-ss` | `EHL` | Score only residues of these SS types. Mostly redundant with the helix-only `ilv_h` tables; matters with `afilmv_ehl`. |
| `--weight-rpx` | `1.0` | Weight of the motif term in `score`. `rpx` is also collected unweighted. |
| `--weight-ncontact` | `0.01` | Weight of the contact-count term in `score`. `ncontact` is also collected unweighted. |
| `--max-bb-redundancy` | `3.0` | Minimum backbone separation (≈ non-aligned RMSD, Å) between two kept docks. **This is what decides how different your `--nout-top` rows are from one another.** |
| `--allowed-residues` | `""` (all) | Residues allowed in the new symmetric interface, in the prosapia positions mini-language (`'20:60'`, `'{h1_start}:{h1_end},{h2_start}:{h2_end}'`). One chain group only. Positions are 1..N **in the file**. |
| `--use-orig-coords` | off | Write the input's real atoms into the dumped docks (see the atom-content trap below). |
| `--set FLAG` | — | Raw flag forwarded verbatim to `python -m rpxdock`, repeatable: `--set '--max_longaxis_dot_z 0.7'`, `--set '--max_cluster 1000'`. Refused if it duplicates a flag above. |

Resource flags (`-c`, `--mem`, `-T`) override the image's `cpu 4 / 16G / 02:00:00`.
**`--mem` is not optional under `afilmv_ehl`** — see the alias section.

---

## The hscore Volume

Modal Volume **`rpxdock-hscore`** (`SAPIA_MODAL_VOLUME_RPXDOCK_HSCORE` in `.env`),
mounted at **`/rpxdock_files`**. Populated by hand; **never delete it**, nothing
repopulates it. Three lowercase alias directories:

| Alias | Size | Contents | Note |
| --- | --- | --- | --- |
| `ilv_h` | 348 MiB, `.txz` only | ILV residues, **helix pairs only**, SS-**in**dependent | the tool default. Correct only for **all-helical** designs. Not pickled |
| `ailv_h` | **1.30 GiB** (453 MiB `.txz` + 879 MiB `.txz.pickle`) | adds Ala, SS-independent | **do not use, untested** — its sidecars came from an unknown stack and may carry the base-pickle defect below |
| `afilmv_ehl` | **14.55 GiB** (5.31 GiB `.txz` + 9.24 GiB `.pickle`) | all SS types, SS-**dependent** | **this campaign's alias. PICKLED 2026-10-02** — loads in **199 s**, was 618 s. `--mem 64G` verified |

**Which one to pass is a campaign decision, not a default** — see
[Choosing the hscore alias](#choosing-the-hscore-alias--a-standing-campaign-decision-read-first) above.

### Pickling an alias — and the upstream bug that makes the naive route fail

Done once for `afilmv_ehl` (2026-10-02): **618 s → 199 s per task, 3.10×**, with
output verified **bit-identical** to the tarballs. The script is
`scripts/rpxdock_pickle/gen_pickles.py` (`generate`, `regen_base`, `check`,
`fresh`, `basecheck`, `promote`). Read this before pickling another alias.

**⚠ Upstream's `--generate_hscore_pickle_files` produces an UNLOADABLE base
table.** The five `hier*` pickles are fine; the `base.rpx` one raises in any fresh
process:

```
AttributeError: Can't get attribute '_PickleWorkaround.flush_only_netcdf_file'
on <module 'xarray.backends.scipy_'>
```

`ResPairScore.rotspace` is an xarray Dataset still backed by the scipy netCDF
store at dump time, so the pickle records a by-name reference to a class xarray
creates **dynamically** inside `_open_scipy_netcdf(..., flush_only=True)`. Nothing
registers that class in a process that only loads pickles. The data is intact —
the failure is serialisation, and it fails **loudly on load**, not silently.

**`.load()` alone does NOT fix it.** Materialising the arrays leaves `ds._close`
bound to `ScipyDataStore.close`, and *that hook* drags the class in. The fix is to
**rebuild the Dataset** from `np.array` copies of its data vars and coords
(keeping attrs), after which `_close` is `None`.

**Verification that actually proves it — in this order:**

1. **Scan the pickle bytes before writing** for `scipy_`, `_PickleWorkaround`,
   `flush_only`, `xarray.backends`. Any hit means the reference leaked; abort
   without touching the staged file.
2. **Load it in a genuinely FRESH process, with no `.txz` loaded first.** This is
   the step that matters: the broken file loads perfectly in a *primed* process,
   so testing the convenient way reproduces the illusion.
3. **Compare content against the original `.txz`** — array shapes, dtypes and md5s
   (79 entries for `afilmv_ehl`'s base; expect 0 differences).
4. **After promoting, re-dock a known scaffold and require bit-identical scores.**
   rpxdock is deterministic, so this is the real end-to-end equivalence test.

**Procedure, and why each step is shaped this way:**

* Generate into a staging dir **outside** the alias dir (`/rpxdock_files/.pickle_staging_<alias>/`).
  A partial set inside the alias dir silently drops every `.txz` (see the trap
  above), so nothing enters until all files are verified.
* **`cd` into the staging dir first** — the flag writes `os.path.basename(f) + '.pickle'`,
  a *relative* path, so it lands in the process cwd. In a prosapia task that is
  `/runs`, i.e. ~10 GiB on the wrong Volume.
* The flag needs `--architecture C2 --inputs1 <any shipped pdb>` to satisfy
  argparse, even though it `sys.exit()`s **before** docking.
* prosapia never commits `rpxdock-hscore` (only the runs Volume), so this needs a
  **one-off Modal function with an explicit `Volume.commit()`**, not a `sapia run`.
* Promote by moving all files **as a set**, then commit. Reverting means moving all
  six back out, again as a set.
* Cost: ~509 s to generate at `--mem 96G`, plus ~9.24 GiB of storage.

**`ailv_h` is plausibly the same defect** — its sidecars came from an unknown
Python/numpy and have never been loaded here. If it is ever needed, run step 2
against its base pickle before trusting it; regenerating locally is likely easier
than diagnosing someone else's dump.

Two traps, both from `get_hscore_file_names` (`rpxdock/score/rpxhier.py:270`):

> #### ⚠ A PARTIAL pickle set silently produces WRONG SCORES
>
> Read before creating any `.pickle` on this Volume. From the source
> (`rpxhier.py:270`): `picklefiles = picklefiles1 or picklefiles2`, then
> `fnames = txzfiles; if len(picklefiles): fnames = picklefiles`. **If even ONE
> `.pickle` exists in an alias dir, RPXdock uses the pickles and ignores every
> `.txz` entirely** — no error, no warning. An alias needs all 6 files; leave 3
> pickles there after a failed or interrupted write and RPXdock scores against a
> third of the motif tables and returns plausible numbers.
>
> The guard that looks like it catches this does not:
> `assert len(txzfiles) in (0, len(txzfiles))` is a tautology. So:
> **generate into a scratch directory, verify the count is 6, and only then move
> all six in as one step.** Never write pickles directly into a live alias dir.
>
> #### ⚠ `--generate_hscore_pickle_files` writes to the CURRENT DIRECTORY
>
> It does `newf = os.path.basename(f) + '.pickle'` — a **relative** path — so the
> files land in the process's cwd, *not* beside the `.txz` files. Upstream's own
> message says "move these into same directory as original .txz files". In a
> prosapia task the cwd is `/runs`, so a naive invocation **dumps ~10 GiB onto the
> runs Volume**. Always `cd` to the intended destination first.
>
> It is a flag on `python -m rpxdock`, not a utility: it loads the alias's tables,
> dumps each one, prints `Faster but non-portable .pickle cache files generated`,
> and calls `sys.exit()` **before any docking**. So it needs no scaffold — but it
> does pay the full table load first (~10 min and the same memory for
> `afilmv_ehl`), plus the write.
>
> Note also that prosapia **only commits the runs Volume** after a task
> (`sapia_modal_task.py` commits `runs_volume` in its `finally`); it never commits
> `rpxdock-hscore`. A write there relies on Modal's implicit end-of-container
> persistence, so do it deliberately with an explicit `commit()`, not as a side
> effect of a `sapia run`.

* **Pickles win over tarballs.** The function globs `*.pickle` *and* `*.txz` and
  then does `fnames = txzfiles; if len(picklefiles): fnames = picklefiles`. So
  `ailv_h` loads its `.txz.pickle` sidecars, never the tarballs. Upstream calls
  those "faster but non-portable"; they were generated elsewhere under an older
  Python/numpy and may fail to unpickle in our 3.12 image. That failure would be
  **alias-specific**: `ailv_h` breaks while `ilv_h` and `afilmv_ehl` work.
  *Workaround:* rename the sidecars on the Volume (`*.txz.pickle` →
  `*.txz.pickle.bak`) so the glob misses them and it falls back to the `.txz`.

  **Measured on the Volume (corrects an earlier figure):** `ailv_h` is **~1.4 GB**,
  not ~475 MB. It holds BOTH the `.txz` tarballs (~475 MB — that is where the old
  number came from) AND the `.txz.pickle` sidecars (~923 MB). **No `.bak` files
  exist: the workaround above has not been applied.** The pickles are therefore
  still what `ailv_h` would load, and they are **untested** — nobody has run a dock
  against this alias in this image. Budget the disk and treat the first `ailv_h`
  run as the test.
* **`WARNING: using slower, portable tarball format…` on every task is NORMAL.**
  It is the `else` branch of that same function, printed for any alias without
  pickle sidecars — i.e. for `ilv_h` and `afilmv_ehl`, always. It is not an error
  and not a misconfiguration.

**`small_ilv_h` is a TEST fixture**, 860 KB, shipped inside the package at
`rpxdock/data/hscore/`. Fine for a smoke test, **never** for a reported score. To
point at it:

```bash
sapia run rpxdock <run_dir> -t table0 -i rfdiffusion3_path --architecture C2 \
  --hscore-files small_ilv_h \
  --hscore-data-dir /usr/local/lib/python3.12/site-packages/rpxdock/data/hscore
```

---

## Columns collected (prefix `rpxdock_`, or `rpxdock_<dir_label>_`)

**The dock**

| Column | Read it as |
| --- | --- |
| `path` | the dock structure (PDB) — the symmetric assembly, all protomers |
| `status` | `OK` or `error: …` |
| `score` | the weighted RPXdock score this row is ranked by, **higher is better** |
| `rpx` | the unweighted motif component |
| `ncontact` | the unweighted contact-count component |
| `model` | index into the run's `_Result.txz` — the pose's stable identity |
| `rank` | 1 = best of this scaffold's kept docks |
| `reslb` / `resub` | 0-indexed first/last scaffold residue present. With `--max-trim 0`: `0` and `n_res-1` |
| `disp` | displacement along the symmetry axis — **cages/dihedrals only, NA for cyclic** |

**Provenance**

| Column | Read it as |
| --- | --- |
| `architecture`, `nfold` | the symmetry docked into, and the order of the axis |
| `hscore` | which motif tables produced `score`/`rpx`. **Never compare rows whose `hscore` differs.** |
| `scaffold_path` | the structure actually docked (post CIF→PDB) |
| `result_path` | the `_Result.txz`, for rescoring upstream-side |
| `n_docks` | poses RPXdock returned after redundancy filtering, before the `--nout-top` cut. Same on every row of a scaffold. `n_docks == nout_top` means **the cut bound** and there may be more worth keeping. |

**Trust metrics — these are the point**

| Column | Why it exists |
| --- | --- |
| `frac_helix` / `frac_sheet` / `frac_loop` | SS of the Body RPXdock scored, from willutil's pure-Python DSSP (no PyRosetta). **Descriptive, NOT a gate — see the measurement below.** Useful for reading a batch (an all-loop body really is a bad scaffold, and `body.py:151` refuses one outright), but do not filter on it. |
| `n_res` | residues in that Body |
| `n_chains_in` | chains in the input **file**. Anything >1 means a concatenation was docked (Premise 4). With the submit-time guard in place this reads `1` on every cyclic dock, and `1` or `2` on a `Dx_y`/cage dock. A `2` on a `Cx` row means the dock predates the guard, or the worker was driven by hand. |
| `input_com_dist` | Å from the origin of the input file's CA centre of mass, measured **before** any recentring. Large + `recentered == False` ⇒ docked in a frame nobody intended. |
| `recentered` | whether `--recenter-input` was applied (and therefore whether the reported transforms are in the recentred frame) |
| `hscore_seconds`, `seconds` | table-load and total wall time, for sizing. Measured: `hscore_seconds` is **33–49 s of a 58–75 s** single-chain dock — most of a short dock is loading motif tables, not searching. |

**NA is not zero.** `disp` is blank for cyclic because it does not apply. A `score`
of 0 is a real, bad value.

### Filters you can now write

**`-f` is a path to a module, not an expression** — see the trap below. Each of
these is a `.py` file in the run_dir, passed as `-f <run_dir>/filters/<name>.py`.

**These examples assume an UNLABELLED run.** With `-l verification` every name
below becomes `rpxdock_verification_…`; see the prefix rule above. Writing them
with a prefix the table does not have is a silent no-op.

```python
# <run_dir>/filters/good_docks.py
# docks worth looking at. NOTE: deliberately does NOT filter on frac_helix --
# see "frac_helix is not a gate" below. Thresholds are alias-specific; these
# suit ilv_h, and afilmv_ehl runs higher.
def apply_filter(df):
    return df[
        (df["rpxdock_status"] == "OK")
        & (df["rpxdock_score"] > 50)
        & (df["rpxdock_ncontact"] > 60)
        & (df["rpxdock_n_chains_in"] == 1)
    ]
```

> ### `frac_helix` is NOT a gate — measured, and it falsifies the old advice
>
> This skill used to claim that because `ilv_h` scores "helix pairs only", a
> non-helical scaffold "scores ~0 everywhere with no error", and the filter above
> used to carry `frac_helix > 0.4`. **Both were wrong.** Measured on
> `denovo_diff_0`, a scaffold with **`frac_helix` 0.00 / `frac_sheet` 0.65**:
>
> | alias | `rpx` on that scaffold | `n_docks` |
> | --- | --- | --- |
> | `ilv_h` | **78.90** — the highest in its batch | 117 |
> | `afilmv_ehl` | **107.74** | 138 |
>
> A 0 %-helix scaffold scored *above* an 81 %-helix one (the `C3_1na0` fixture,
> `rpx` 77.64) under the supposedly helix-only tables. The `_h` in `ilv_h` refers
> to how the tables were GENERATED; as the alias table says, they are
> SS-**in**dependent at scoring time. **A `frac_helix > 0.4` filter would have
> discarded the best-scoring design in the batch.**
>
> What the numbers do support: `afilmv_ehl` scores this β-rich scaffold ~37 %
> higher than `ilv_h`, which is consistent with SS-dependent tables having more to
> say about sheets — one more reason for the campaign's alias choice, but on
> different grounds than "`ilv_h` returns zero".
>
> Still open: nobody has checked willutil's DSSP against an independent
> assignment, so an alternative reading is that the SS call itself is wrong on
> rfd3 backbones. The fixture's clean 0.81 on a known helical bundle argues
> against that, but it is untested.

```python
# <run_dir>/filters/best_dock.py
# only the best dock of each scaffold
def apply_filter(df):
    return df[df["rpxdock_rank"] == 1]
```

```python
# <run_dir>/filters/best_dock_labelled.py
# the same filter on a run submitted with `-l verification`
def apply_filter(df):
    return df[
        (df["rpxdock_verification_status"] == "OK")
        & (df["rpxdock_verification_rank"] == 1)
        & (df["rpxdock_verification_n_chains_in"] == 1)
    ]
```

```python
# <run_dir>/filters/cut_bound.py
# scaffolds where the --nout-top cut bound (more poses may be worth keeping)
def apply_filter(df):
    return df[df["rpxdock_n_docks"] > df["rpxdock_rank"].max()]
```

```python
# <run_dir>/filters/suspicious_frame.py
# suspicious frame
def apply_filter(df):
    return df[
        (df["rpxdock_input_com_dist"] > 5)
        & (df["rpxdock_recentered"] == False)  # noqa: E712
    ]
```

```bash
sapia run <tool> <run_dir> -t table1 -f <run_dir>/filters/good_docks.py
```

---

## Traps

### A multi-chain input HANGS for >10 minutes before docking even starts
**The most expensive failure found in this tool, and now guarded.** A 2-chain,
154-residue input made every array task sit in willutil's PDB parser —
`willutil/pdb/pdbfile.py:295`, `df.ri.iloc[i]` in a per-atom pandas loop — for over
ten minutes, **before any rigid-body sampling began**, and the whole array's wall
clock went with it.

**Signature, if you ever meet it again (e.g. on a worker driven by hand):**

| Evidence | What you see |
| --- | --- |
| `<script>_<id>.out` | the opening `[HH:MM:SS] task N: rpxdock on <name> (<path>)` line and nothing else |
| `<script>_<id>.err` | only the hscore `WARNING: using slower, portable tarball format…` and the `using hscore …` line |
| `.exit` file | **absent** — the task is neither done nor failed |
| `<out_dir>/<name>/` | exists and is **empty** — no dumped PDB, no `_Result.txz` |
| a stack dump / `py-spy` | stuck at `willutil/pdb/pdbfile.py:295` in `df.ri.iloc[i]` |

It was also **scientifically wrong**: RPXdock collapses every chain into one chain
renumbered from 1, so a 2-chain dimer docked into `C2` becomes a single
concatenated pseudo-chain of double the length — not a `C2` dock of anything.

**The guard — submit time, in the workstation, seconds. The limit is per
architecture family:**

| Architecture | Protocol | Chains allowed | Over the limit |
| --- | --- | --- | --- |
| `C2`…`C17`, `CxSTACK` | `cyclic`, `cyclic_stack` | **1 only** — a single-chain monomer | **hard error** |
| `Dx_y` (`D3_2`, `D3_3`, …) | `onecomp` | **1 or 2** — 2 is legitimate, a `Dx_2` scaffold is a dimer on the 2-fold | **hard error at 3+** |
| `T2 T3 O2 O3 O4 I2 I3 I5` | `onecomp` | **1 or 2** (same rule) | **hard error at 3+** |

**There is no 2-chain warning.** Two chains into a `Cx` is refused outright; two
chains into a `Dx_y` or cage passes silently. Both errors name the design, the
observed chain count and the architecture, and the run aborts — nothing is
submitted.

To fix an offending input, run `chainsel` first and dock its `chainsel_path`.

**Reference timings, measured, single chain — and these are `ilv_h` timings:**
100 residues and 119 residues both dock end to end in **58–75 s**. Under
`afilmv_ehl` (this campaign's alias) they do not apply at all; re-baseline from a
single-design probe before calling anything a hang — see the alias section. Of that, **loading the `ilv_h` motif tables alone
is 33–49 s** (`rpxdock_hscore_seconds`) — so most of a short dock is table loading,
not searching. Raising `--beam-size` buys search time against a fixed ~40 s
overhead; batching many scaffolds per container would amortise it, which this tool
deliberately does not do (one scaffold per task).

### Other inputs refused at submit time
The same guard opens the file RPXdock will actually read (gzip handled) and
refuses, by name, with the design and the path:

| Refused | Why |
| --- | --- |
| a suffix that is not `.pdb`, `.pdb.gz`, `.cif`, `.cif.gz`, `.mmcif` | the input column is not pointing at a structure |
| no parseable `ATOM` records in the first model | empty, HETATM-only, or not a PDB |
| `ATOM` records that do not parse as fixed-width PDB | truncated or hand-edited file |
| more than one model | RPXdock docks model 1 and silently drops the rest |
| any residue missing `N`, `CA` or `C` | a CA-only or gapped backbone either crashes in `body.py` or scores ~0 with no error |
| residue numbering going backwards, a residue split across two blocks, or a jump > 1000 within a chain | the renumbered 1..N positions `--allowed-residues` speaks in would not be the ones you counted |

**Warning only** (printed, still submitted): a missing backbone `O` — willutil
guesses it with `add_bb_o_guess`, but the secondary structure, and therefore the
motif score, is then partly computed on guessed atoms.

A **missing input file** is still a skip with a printed message, not an error —
that is a gap in the table, not a wrong instruction.

### `.pdb.gz` is handed to RPXdock as-is
Routing converts **only** mmCIF. `prosapia.utils.ensure_pdb` short-circuits on
`suffix == ".pdb"`, and a `.pdb.gz`'s suffix is `".gz"`, so calling it
unconditionally round-tripped every gzipped PDB through gemmi for nothing. RPXdock
reads gzipped PDB natively — its own shipped fixture is `C3_1na0-1_1.pdb.gz`.
`rpxdock_scaffold_path` tells you which file was actually docked.

### A guard refusal is INVISIBLE to `$?`, and leaves a phantom table

Two things measured while verifying the guard (2026-10-02). Both matter to anyone
automating a submit, and neither is rpxdock's fault alone.

**1. `sapia modal-shell --cmd` returned exit code 0 even though the inner
`sapia run` raised `ChainCountError` and queued nothing.** The refusal was visible
only in the output TEXT.

**Confirmed by direct probe, not inference:** `--cmd 'exit 3'` → `rc=0`, and
`--cmd 'ls /nonexistent'` → `rc=0` with the real `No such file or directory` on
stderr. prosapia is not at fault (`modal_shell_from_args` is literally
`raise SystemExit(subprocess.run(...).returncode)`); the code is lost inside
`modal shell --cmd` or prosapia's base64 `bash <(...)` wrapper. CLAUDE.md's old
claim that `--cmd` "exits with its code" was wrong and has been corrected.

**So never conclude a submit succeeded because `$?` was 0.** Check the output for
`Submitting N designs`, or — authoritative, and the only check available once the
stdout is gone — for `<out_dir>/<script>_logs/<script>_modal.json` holding
`{"app_id", "n_tasks"}`. No file means nothing was queued. An orchestrator that
trusts the exit code will march straight on to `collect` and report "0 rows" as
though the designs were simply bad.

**2. A refused submit still registers the child table.** The table name is
reserved before the manifest is built, so a run that then refuses leaves:

* a row in `_registry.tsv` (e.g. `table1_guard_dimer_c2`, parent `table0`, tool
  `rpxdock`),
* an output dir holding only `.meta.json` and empty `configs/` and
  `<leaf>_logs/`,
* and **no `<table>.tsv` at all.**

So `_registry.tsv` can list a table that does not exist. When auditing lineage,
**a registry row with no matching `.tsv` is a refused or abandoned submit, not a
table with zero rows.** It costs no compute; it is just debris. Delete the leaf
dir and leave the registry row, or clean both — but know which you are looking at.

### `-f/--filter` is a module path, not an expression
prosapia declares `-f` as `type=Path` and loads it with
`importlib.util.spec_from_file_location`; the module must define
`apply_filter(df) -> df` or it raises `AttributeError: … does not define an
'apply_filter' function`. A pandas-style string dies at **submit** time with
`ImportError: Could not load module from <your expression>` — measured. Write the
filter into the run_dir (on Modal, through the workstation with a `--cmd` heredoc;
a local path does not exist in the container) and pass that path. Since the module
sees the whole DataFrame before the manifest is built, it is also the supported way
to dock **one named scaffold** — `return df[df["name"] == "design_3"]` — before
committing a batch.

### A failed scaffold mints **no row at all**
This is a `create` tool: a scaffold whose dock failed has no entity to key a row
to, so it simply does not appear in the child table (same contract as
`bindcraft2`). The failure *is* recorded — the task exits non-zero, and
`<out_dir>/<scaffold>.tsv` holds `status: error: …` — but the table alone will not
show it. **To audit a run, compare the child table's distinct `parent_name` values
against the parent table's rows.** A parent with no child either errored or
produced no non-clashing dock.

### `RuntimeError: bad strides, strides not supported` means *no viable dock*
Raised from `bvh_collect_pairs_range_vec`. Verified cause: the hierarchical beam
emptied — at some stage nothing scored above zero, so an empty `(0,4,4)` array
reached the C++ binding. It is a **"nothing docked" outcome wearing a crash's
clothes**, not a broken install. Reproduced deterministically when docking a cage
against the toy `small_ilv_h` table. If you see it with real tables, the scaffold
has no interface the motif tables recognise (check `frac_helix`), or the
architecture does not suit it.

### `assert … 'body is all loops and not sub-body!!'`
`body.py:151`. The input's SS came back entirely loop. That is a coordinate or
backbone problem in the input (missing O atoms, a Cα-only file), **not** a docking
problem. Note willutil fills a missing O with `wu.chem.add_bb_o_guess`, so this
fires only when the backbone is genuinely unusable.

### A dumped dock is **backbone-only** by default
`make_pdb_from_bodies` writes `N CA C O CB` plus a pseudo-atom named `CEN` —
measured: 6 distinct atom names, 1428 atoms for a 119-residue C2. That is fine for
a backbone look and for redesign, but it is **not an all-atom structure**: anything
that counts atoms, packs side chains or computes buried area will be misled. Pass
`--use-orig-coords` for the real atoms. Residue numbering in the dump is
**1..N per chain**, not the input's.

> #### How to tell whether `--use-orig-coords` actually took
>
> **Use atoms per residue and the `CEN` count. Do NOT use a distinct-atom-name
> count** — an earlier revision of this page claimed "83 atom names, 3752 atoms"
> for an all-atom dump, and that is wrong. 3752 atoms over 238 residues is ~15.8
> per residue, which is only reachable if **hydrogens** were counted; these dumps
> carry none. ~35 is the correct number of distinct heavy-atom names across the 20
> amino acids, so the name count barely moves and is **not diagnostic**.
>
> | | atoms/residue | `CEN` count |
> | --- | --- | --- |
> | backbone-only (default) | **~5** (N CA C O CB) | one per residue |
> | `--use-orig-coords` | **~8.4** (real side chains) | **0** |
>
> *Measured 2026-10-02* on `outputs/20261002_143419_dimer_phase2`: a 111-residue
> C2 dump gave 1832 ATOM records over 222 residues = **8.25 atoms/residue, 0 `CEN`**.
> A distinct-name count on the same file gave 35, not 83.

### Cages and dihedrals are **not verified end to end**
Cyclic (`C2`, `C3`) is verified in the real image, start to finish. The
one-component cage/dihedral path is not, and has a history:
* under **numpy 2** it dies for every input at `sampling/xform_hier.py:59`
  (`int(np.ceil(ang / angresl))` on a 1-element array);
* under the pinned **numpy 1.26** it builds the sampler and runs, but the only
  tables available for testing (`small_ilv_h`) emptied the beam and produced the
  `bad strides` signature above;
* RPXdock's own `tests/search/test_onecomp.py` opens with
  `pytest.importorskip('pyrosetta')`, so it does not cover this configuration
  either.

**Run one scaffold before fanning out a cage campaign.** Same for
`--docking-method grid` (works, but returns thousands of poses — the C2 smoke test
gave 3222 — so `--max-bb-redundancy` and `--nout-top` matter much more there).

### `--score_only_sspair` is not exposed, on purpose
The image pins `numpy<2` (otherwise cages/dihedrals are dead on arrival).
`Body.filter_pairs` uses `np.ones(..., dtype=np.bool)`, and `np.bool` exists in
numpy 2 but **not** in 1.26 — so that option raises
`AttributeError: module 'numpy' has no attribute 'bool'` under the pin. It is
reachable through `--set '--score_only_sspair HH'` if you accept that. (Measured
separately: with the cyclic protocol it appeared to be a no-op even where it ran,
producing identical scores.)

### `pytest --pyargs rpxdock` is weak evidence
Several search tests (`test_cyclic.py`, `test_onecomp.py`, `test_multicomp.py`,
`test_plug.py`, `test_axle.py`, `test_asym.py`, `test_result.py`) begin with
`pytest.importorskip('pyrosetta')` and **skip** in this image. A green run proves
very little. The real acceptance check is a dock that produces a `_Result.txz` with
non-degenerate scores — which is why the image build runs one.

### No PyRosetta, by design
`--use_rosetta` is a `store_true` on RPXdock's CLI, so the CLI already defaults to
the PyRosetta-free path (verified in-container). Consequences:
* SS comes from willutil's `wu.dssp` on the backbone N/CA/C/O, not Rosetta's;
* the helix-termini options (`--term_access*`, `--termini_dir*`) print
  `no pyrosetta, helix termini stuff diabled` and are skipped;
* **`rp.search.result_from_tarball` cannot reopen a result's bodies** —
  `AttributeError: module 'rpxdock.rosetta.triggers_init' has no attribute
  'pose_from_file'` (verified). This is why the worker runs the dock in-process and
  dumps structures from the live `Result` rather than reopening the tarball. If you
  want to re-dump poses from a `result_path` later, do it somewhere with PyRosetta.
* **Footgun for anyone editing the worker:** `Body.__init__`'s own signature
  default is `use_rosetta=True` — the opposite of the CLI's. Code that constructs a
  `Body` directly must pass `use_rosetta=False` explicitly.

### Image build is slow the first time
RPXdock's C++ extensions are JIT-compiled by `cppimport` on first import. The image
pre-builds them and then sets `CPPIMPORT_RELEASE_MODE=1` so nothing compiles at run
time.

**Measured, whole image build: ~258 s (4.3 min)** — not the ~7 min quoted earlier
in this skill, which was the bare extension-compile time timed by hand in a
container. Its self-verification passed:

```
extensions OK: 12
rpxdock image OK: numpy 1.26.4 ndock 18 best 10.831188201904297
```

Still allow 10+ minutes for the first `sapia run rpxdock -e modal` (build plus
queue), and seconds thereafter.

If a future base image moves to g++ 14 or newer, the build will fail in
`rpxdock/extern/Eigen` with `'Eigen::Transpose<…>' has no member named 'derived'`
(measured on g++ 15). Debian bookworm's g++ 12 compiles it cleanly. Pin the
compiler rather than patching Eigen.

---

## What to do with the docks

A dock is a symmetric assembly built from the *parent's* backbone, so the obvious
next steps are:

| Question | Next |
| --- | --- |
| What sequence makes this interface real? | `proteinmpnn` / `atomium` with `-i rpxdock_path` (and `--set '--ca_only'` is **not** needed — the dump has N/CA/C/O/CB; pass `--use-orig-coords` if a tool needs side chains) |
| How big is the new interface? | `cms` with `-i rpxdock_path`, naming two protomer chains |
| Is it well packed? | `pyrosetta` with `-i rpxdock_path` |
| Did the designed sequence fold back to it? | `boltz` on the designed sequence, then `usalign --col-a boltz_path --col-b rpxdock_path` |

---

## Sanity-checking a change to this tool

There is no test suite in this workspace. What was run while building it, and what
to re-run after an edit:

```bash
# discovery + flags (from the repo root, so PROSAPIA_TOOLS_DIR resolves to tools/)
./.venv/bin/sapia run rpxdock --help
./.venv/bin/sapia collect rpxdock --help

# the pure logic, no image needed — import the modules directly and check:
#   classify_architecture  over every accepted form (C2, c3, C17, C3STACK, T3,
#     O4, I5, D3_2, D3_3, D8_8) and every refused one (T32, I32, O43, ASYM,
#     AXLE_3, PLUG_C3, P6_632, P4M_4, F_32, C1, C99, D3, D3_4, XY, "")
#   parse_set_flags        passthrough, multiple tokens, the OWNED_RPX_OPTIONS
#                          collision, a token that is not a flag
#   build_rpxdock_manifest 3 fields per row, gpus_per_task forced to 0, the
#                          config JSON's protocol/nout_top/argv, the missing
#                          -i error, --nout-top 0, and {expr} resolving to a
#                          different allowed-residues file per design
#   collect_rpxdock        a real worker TSV -> N Collected with parent set,
#                          floats/ints/bools cast, `disp` NA on cyclic, and a
#                          failure TSV yielding zero rows
```

`tools/rpxdock/rpxdock_structure.py` is **stdlib-only and imports nothing from
prosapia, numpy or the image** — that is deliberate, so the guard can be exercised
against handwritten PDB text in a scratch dir. It is also the single parser behind
`rpxdock_input_com_dist` and `rpxdock_n_chains_in`, shared with the worker, so a
change to it moves two collected columns. Re-check, as pure functions over a file
path:

```
#   input_kind(path)            -> 'pdb' for .pdb and .pdb.gz (case-insensitive),
#                                  'cif' for .cif/.cif.gz/.mmcif,
#                                  InputStructureError for anything else
#   resolve_input_pdb(...)      a .pdb and a .pdb.gz come back UNCHANGED and no
#                                  run_dir/.cif_to_pdb/ is created; a .cif comes
#                                  back as a converted file under .cif_to_pdb/
#   scan_structure(path)        n_chains / chain_order / n_residues / ca_coords on
#                                  a plain and a gzipped file; first model only;
#                                  HETATM ignored; MODEL records counted past the
#                                  first ENDMDL
#   validate_structure(scan)    RAISES InputStructureError on: an empty file, a
#                                  HETATM-only file, a truncated ATOM line, two
#                                  MODELs, a CA-only file (missing N/C), residue
#                                  numbering that goes backwards, a residue split
#                                  across two non-adjacent blocks, a within-chain
#                                  jump > MAX_RESNUM_JUMP (1000).
#                               RETURNS one warning string for a missing backbone
#                                  O, and [] for a clean file. A single
#                                  MODEL/ENDMDL pair is clean.
#   check_chain_count(...)      cyclic + cyclic_stack: 1 chain passes, 2 and 3
#                                  RAISE ChainCountError;
#                               onecomp (Dx_y and the cages): 1 and 2 pass, 3
#                                  RAISES. Returns None — there is no warning case.
#   check_input_structure(...)  the two above in one call over a path; warnings
#                                  returned, refusals raised
#   build_rpxdock_manifest      a 2-chain input + --architecture C2 aborts the
#                                  whole run (ChainCountError, zero rows written);
#                                  the same input + --architecture D3_2 submits;
#                                  a MISSING file is still SKIPPED with a printed
#                                  message, not raised
#   rpxdock_worker.input_frame_metrics   unchanged numbers for a known file —
#                                  it now reads through scan_structure
```

The end-to-end dock needs the image. The image build itself runs a C2 dock on
RPXdock's shipped `C3_1na0-1_1.pdb.gz` + `small_ilv_h` and **fails the build** if
it produces no scoring pose, if any of the 12 compiled extensions is missing, if
PyRosetta turns up, or if `tar`/`bzip2` are absent — so a green build is already a
meaningful check. To exercise the worker itself without Modal, build the image
locally from `modal_image.py`'s steps and run:

```bash
python tools/rpxdock/rpxdock_worker.py --name d0 \
  --input <site-packages>/rpxdock/data/pdb/C3_1na0-1_1.pdb.gz \
  --config <a config json like the manifest builder writes> --out-dir out/
```
