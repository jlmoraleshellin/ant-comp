---
name: ssprofile
description: >-
  How compact is this chain, and is it helix, strand or loop? Takes any structure column and a
  REQUIRED --design-chains selection, and writes per design the CA radius of gyration
  (MASS-UNWEIGHTED) plus 3-state DSSP helix/beta/loop fractions, measured on the selected chain(s)
  ALONE. Load before composing an ssprofile run, before checking whether a helicity weight or SS
  preset actually changed the designs, and before reading any ssprofile_ fraction. The traps:
  the selection is assigned in isolation, so a binder strand paired only with a TARGET strand
  reads as loop; and ssprofile_sequence must be checked first, because a wrong chain letter
  (an AUTHOR chain ID) gives tidy numbers for the wrong protein. Covers the OK / warn / error contract (warn = Rg measured, fractions NA),
  why the fractions match rpxdock's frac_* basis, and why the skill offers no thresholds.
---

# ssprofile

The CA radius of gyration and 3-state secondary-structure composition of a chosen chain
selection. **`action: update`**: it annotates the table it runs on and mints no child.
`-t` is required. CPU-only: the builder forces `gpus_per_task = 0`, so no `-g 0` is needed.

## Premise, and when its numbers are meaningless

It assumes one structure file per design, and a `--design-chains` selection that is
**protein**: each selected chain is **one continuous polypeptide with a complete N/CA/C/O
backbone**. The numbers describe **that selection's own backbone, and nothing else**.

When that premise does not hold, here is what you get:

| Situation | What it does |
| --- | --- |
| A selected chain letter is not in the file | `error:`, every measured column NA (`input_column` and `source_path` are still recorded). It never falls back to another chain, and never returns a smaller selection |
| **Any** selected chain has no amino-acid residue with a CA (empty chain, ligand/ion-only chain, unreadable file) | `error:`, every column NA. With `A,L` where `L` is ligand-only, the whole row is an error, not a quietly smaller selection |
| An amino-acid residue is missing N, CA, C or O | `warn:` naming the residue (`A16(CA)`, `A12(O)`). Rg, `n_res` and `sequence` are measured, **fractions NA**. A missing O is never guessed (willutil would guess it; this tool refuses). A CA-less residue is left out of `n_res` and Rg, **including one at a chain end**, which the break test alone could not see |
| A C(i)–N(i+1) gap > 2.5 Å inside a selected chain | `warn:`. Rg, `n_res` and `n_chain_breaks` are measured, **fractions NA** |
| The DSSP call itself raises | `warn:`, fractions NA |
| **The WRONG chain is selected** (the target instead of the binder) | **`OK`, with plausible numbers for the wrong protein.** The tool cannot detect this. Read `sequence` (and `n_res`) |
| **The structure is a prediction whose fold you have not checked** | **`OK`.** The fractions then describe a misfold faithfully |
| A one- or two-residue selection | `OK`, with Rg 0 (or tiny) and all loop. `n_res` tells you why |

In short: **`OK` is the only status where all four science columns are measured.** It does
not mean the selection was the right one.

## Trust columns: read these before the science

| Column | Expected value | What a surprise means |
| --- | --- | --- |
| `ssprofile_chains_used` | The chain(s) you meant, expanded (`A:C` comes back as `A,B,C`) | A different list means a different flag was passed than you think |
| `ssprofile_sequence` | **The selected chain's sequence**, one-letter, chains joined by `/`. Modified residues are lower-case (`MSE` becomes `m`); anything without a letter is `X` | **The wrong-chain guard.** `n_res` cannot catch a wrong chain when binder and target have similar lengths; the sequence can. A target's known sequence here means you profiled the target. On a backbone generator's raw output this is whatever names it wrote (often all `G`), which still tells binder from target. Compare it against the designed sequence column when there is one |
| `ssprofile_n_res` | **The length you expect the selected chain to have**: the binder length the generator was asked for, or the chain length you know | A target-sized `n_res` means you profiled the target. A count short of the expected one, on an `OK` row, means the chain itself is shorter; CA-less residues give a `warn:`. Record the expected length in the campaign log |
| `ssprofile_n_chain_breaks` | `0` on every `OK` row | `> 0` only on a `warn:` row. NA on a missing-backbone-atom `warn:` (the gap test needs N and C) |
| `ssprofile_input_column`, `ssprofile_source_path` | The `-i` column you passed, and the file it resolved to on that row | **Provenance, written on every row including `error:` rows.** Two runs that shared a label, or a column filled partly by one run and partly by another, show up here as rows that disagree. `source_path` is the file actually read: a `.cif` input appears as its staged PDB under `<run_dir>/.cif_to_pdb/` |

**Chain letters are AUTHOR chain IDs (`auth_asym_id`), as a viewer shows them, for both PDB
and mmCIF input.** gemmi reads author names, and a `.cif` is staged to PDB with those names.
mmCIF also carries `label_asym_id`, which can differ: 1TUP's author chains `E,F,A,B` are label
chains `A,B,C,D`. A tool or spec that speaks in label IDs (boltzgen's epitope frame, for
example) does not name the same chain with the same letter. Pass the author letter here.

**Chain letters are a convention of the generator, not of the biology.** Find the binder's
letter per source before running:

- bindcraft2: read `bindcraft2_binder_chain`, and check `bindcraft2_binder_chain_src`. The
  `last_chain` fallback is an under-count for a multi-chain binder (see the bindcraft2 skill).
- rfd3, boltz, boltzgen: whatever the run's own spec assigned. Confirm with `sequence` and
  `n_res`.

## The question it answers, and the threshold

**Is this chain compact, and how much of it is helix, strand or loop?** Typical uses:

- checking that a helicity weight or SS preset (for example bindcraft2's
  `weights_binder_helicity` or `mixed_topology`) actually moved the distribution;
- comparing a designed backbone with its prediction (run twice with different `-l`);
- spotting extended or unfolded chains, which show a high Rg for their `n_res`.

**This skill gives no threshold, deliberately.** House rule: a tool records scores and never
filters, and a skill ships no cutoff. Rg grows with chain length, so a bar that holds at one
`n_res` is wrong at another. SS fractions depend on what the campaign asked for. **Read the
distribution in your own table and choose a bar yourself**, then record it and the reason in
the campaign log.

### Exact definitions

- **`radius_of_gyration`**: CA-based and **mass-unweighted**. Take the selection's CA
  coordinates `r_i` (one per amino-acid residue), let `centroid = mean(r_i)`, and compute
  `Rg = sqrt(mean(|r_i − centroid|²))`, in Å. Every CA weighs the same, and residue identity
  plays no part. This matches the informal Rg comparison used earlier in this campaign; keep
  it that way so later comparisons stay on one basis. With several chains selected it is the
  Rg of their **union**, not an average of per-chain values. Reported to 4 decimals.
- **`helix_fraction` / `beta_fraction` / `loop_fraction`**: counts of mdtraj simplified DSSP
  codes, divided by `n_res`. **H includes 3-10 (`G`) and π (`I`) helix. E includes isolated
  β-bridges (`B`). Loop is turn `T`, bend `S` and coil.** They sum to 1.0 on an `OK` row.

### Why these fractions sit on the same basis as rpxdock's `frac_*`

rpxdock's `frac_helix`, `frac_sheet` and `frac_loop` come from willutil's `wu.dssp`. Read at
willutil commit `ded19f6` (`willutil/chem/dssp.py`), that function is **not** a DSSP
implementation, despite rpxdock's docs calling it "pure-Python". It writes N/CA/C/O to a
temporary PDB as **poly-GLY** and calls `mdtraj.compute_dssp(simplified=True)`, then renames
`C` to `L`. ssprofile makes the same call on the same poly-GLY N/CA/C/O input, built in memory
instead of through a file. The tests check it against a re-implementation of willutil's file
round-trip, character by character.

Two consequences of that convention:

- **Proline is treated as an H-bond donor**, which strict DSSP does not do. Measured against
  the same call with real residue names on 1UBQ, 2GB1, 1CRN and 4HHB (319 residues,
  15 prolines): **one residue differed in total** (in 1CRN).
- **The mdtraj version behind rpxdock is not recorded.** rpxdock's image installs willutil
  from git with mdtraj unpinned; this image pins `mdtraj==1.11.0`. So the parity is "same code
  path", not yet "same version".

## Invocation and flags

```bash
sapia run ssprofile <run_dir> -t table0 -i rfdiffusion3_path --design-chains B -l bb
sapia collect ssprofile <run_dir> -t table0 -l bb
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `-i/--input-column` | **effectively required** | The default is the sentinel `"not applicable"` (the cms/chainsel/epitope convention), and the builder **raises** unless you name a column. Any structure column works: `rfdiffusion3_path`, `boltz_path`, `bindcraft2_path`, `chainsel_path`, and so on. The column is **resolved row first, then up the lineage**, so `-i rfdiffusion3_path -t table1` reads the parent's backbone. `.cif`/`.cif.gz` are staged to PDB up front (cached in `<run_dir>/.cif_to_pdb/`) |
| `--design-chains` | **required, no default** | The chains to profile, in the chain mini-language: `,` separates and `:` is an inclusive letter range (`D`, `A,B`, `A:D`). An empty or repeated chain raises at submit time. With several chains, Rg is of the union and the fractions are over all their residues together, with DSSP seeing every selected chain (strands paired *between* selected chains count as E) |
| `--designs-per-task` | `20` | Designs per task. Each design takes well under a second |
| `-l/--dir-label` | none | **Required whenever ssprofile measures a second thing on one table** (backbone vs prediction, or two chain selections). See *Reruns* below for what happens without it |
| `-f/--filter`, `--force` | | The usual. Already-`OK` rows are skipped on a rerun unless `--force` |

**Reruns without a new `-l`** write to the same leaf dir and the same columns. What happens
depends on the rows, not on intent:

- **No `--force`, every row already `OK`:** nothing is submitted (`No designs to submit.`). The
  old values stand, even if you changed `-i` or `--design-chains`.
- **No `--force`, some rows not `OK`:** only those rows rerun, **with the new settings**, while
  the `OK` rows keep the old run's values. The result is **a single column holding two
  different measurements**. `ssprofile_input_column`, `ssprofile_source_path` and
  `ssprofile_chains_used` disagree between the two groups of rows; check them after any rerun.
- **`--force`:** every row reruns, and every value is replaced by the new measurement.
- **Stale files cannot leak into a rerun.** The task script deletes each design's `<name>.tsv`
  and `<name>_residues.tsv` before the worker starts. A task killed mid-run (OOM, timeout,
  pre-emption) therefore leaves `error: worker crashed` or `missing`, never an earlier run's
  `OK`.

To measure something *different* on the same table, use a new `-l`.

Resources: 2 CPU, 4 GB, 30 min per task. No GPU.

**Executors.**

- **Modal:** `modal_image.py` is `debian_slim` + `gemmi`, `numpy` and `mdtraj==1.11.0`.
  mdtraj pulls in scipy. The build checks that the compiled DSSP kernel imports.
- **vib: NOT runnable yet.** There is no `SAPIA_ACTIVATE_SSPROFILE` in `.env`, and the
  workspace venv that `activation/vib/prosapia_python.sh` activates **does not have mdtraj**.
  Porting it means adding mdtraj to that venv (via `scripts/vib_bootstrap.sh`) and then adding
  the `.env` line. Until then, the tool appears in `sapia run --help` on vib, and every task
  dies at activation with `set SAPIA_ACTIVATE_SSPROFILE`.

## Gotchas, each with its signature

- **The selection is assigned in isolation.** A binder β-strand that pairs only with a target
  strand gets **`beta_fraction` 0 and a higher `loop_fraction`**, with status `OK`.
  Signature: a binder that looks like a strand in PyMOL (on the full complex) reads as loop
  here. That is by design: the question is the chain's *own* structure. If you want the strand
  counted, select both chains (`--design-chains A,B`), and accept that the numbers are then
  about both.
- **Wrong chain letter.** Signature: `OK` everywhere, `sequence` equal to the target's
  sequence on every row (often with `n_res` equal to the target's length), and fractions that
  barely vary across designs (the target is fixed). Check `sequence` before anything else:
  `n_res` alone misses this when binder and target have similar lengths. Also check that you
  passed an author chain ID, not an mmCIF `label_asym_id`.
- **An unexpanded chain range reaching the worker.** The builder expands `A:B`. If you call
  `ssprofile_worker.py` standalone with `A:B`, you get
  `error: chain(s) A:B absent ... (have: A,B)`. Pass `A,B` there.
- **A refused submit leaves a phantom table** (a CLAUDE.md known gap). A missing `-i`, an
  empty `--design-chains`, or a column that resolves on no row raises while the manifest is
  built. That happens *after* the registry row and leaf dir exist. Signature: no
  `Submitting N designs` line, and no `<leaf>_logs/<script>_modal.json`.
- **Unresolvable rows are submitted, not dropped.** A design whose input column is empty on
  its row and on every ancestor still goes into the manifest, with an empty path. Signature:
  `<name>: could not resolve <col> (row or lineage) -- submitting anyway` at submit, and
  `error: could not resolve a structure path from column '<col>'` in the table. If **every**
  row fails to resolve, the submit raises instead.
- **A `warn:` row is not a failure of the worker.** Signatures:
  `warn: N residue(s) missing backbone atoms, SS not assigned: A12(O),...` and
  `warn: N backbone break(s) (C-N > 2.5 A after A57,...), SS not assigned`. Rg is real on
  those rows; the fractions are absent on purpose.
- **Multi-model files** (NMR ensembles, some prediction dumps): only the **first model** is
  read, silently. Alternate locations keep the first conformer; hydrogens are removed.
- **Lowercase or blank chain IDs cannot be selected.** The chain mini-language upper-cases
  letters, so a file whose binder is chain `b` gives `error: chain(s) B absent ... (have: A,b)`.
  Rename the chain upstream (`chainsel --rename-to`).
- **A worker crash** (for example mdtraj not importable) writes
  `error: worker crashed` for every design of that task, and the task exits non-zero.
- **The collector re-checks the contract.** An `OK` TSV whose fractions do not sum to 1, or a
  `warn:` TSV carrying fractions, is collected as `error: ...` with NA. Seeing one means the
  worker or the file was tampered with or broken. Do not hand-edit result TSVs.

## What it collects

Column names are `ssprofile_<field>` (or `ssprofile_<label>_<field>` with `-l`).

| Column | Meaning | How to read it |
| --- | --- | --- |
| `ssprofile_status` | `OK` / `warn: ...` / `error: ...` / `missing` | **`OK` is the only status where all four science columns are measured.** `warn:` means Rg is measured and the fractions are NA. `error:` and `missing` mean every measured column is NA (provenance is kept on `error:`) |
| `ssprofile_chains_used` | **Trust.** The chains profiled, comma-joined | Should be exactly what you meant |
| `ssprofile_sequence` | **Trust.** One-letter sequence of the selection, chains joined by `/` | The wrong-chain guard (see *Trust columns*). Measured on `OK` and `warn:` rows, NA on `error:` |
| `ssprofile_input_column` | **Provenance.** The `-i` column this row was measured from | Recorded on every row the worker wrote, `error:` included. NA only for `error: worker crashed` and `missing` |
| `ssprofile_source_path` | **Provenance.** The file actually read (after CIF-to-PDB staging) | Same. Empty when the column resolved to nothing on the row or its ancestors |
| `ssprofile_n_res` | **Trust.** Amino-acid residues with a CA in the selection | The CA count behind Rg and the denominator of the fractions. Compare with the expected chain length |
| `ssprofile_n_chain_breaks` | **Trust.** In-chain C–N gaps > 2.5 Å | `0` on `OK`. `> 0` explains a break-`warn:`. NA on a missing-atom `warn:` |
| `ssprofile_radius_of_gyration` | CA Rg in Å, mass-unweighted (definition above) | Scales with chain length: **compare only at similar `n_res`**. A high Rg for its `n_res` means an extended or partly unfolded chain |
| `ssprofile_helix_fraction` | DSSP H+G+I / `n_res` | 0.0 on an `OK` row is a measured zero. Empty means not measured |
| `ssprofile_beta_fraction` | DSSP E+B / `n_res` | Same. **Blind to strands paired only with unselected chains** |
| `ssprofile_loop_fraction` | Everything else / `n_res` | Same |
| `ssprofile_path` | `<out_dir>/<name>_residues.tsv`: one row per residue (`chain, resnum, resname, ss`) in the order Rg and DSSP saw them | Open it to see *where* the helices are. `ss` is empty on a `warn:` row. Empty on `error:` |

On disk: `<run_dir>/<table>/ssprofile[_<label>]/` holds `<name>.tsv`,
`<name>_residues.tsv`, `ssprofile_tasks/task_<k>.tsv` (the sub-manifests) and
`ssprofile_logs/`.

## A filter you can write

It filters on status only. The house convention is that no metric bar ships in a skill; pick
yours from the observed distribution.

```python
# <run_dir>/filters/ssprofile_ok.py
def apply_filter(df):
    return df[df["ssprofile_status"] == "OK"]
```

```bash
sapia run <next_tool> <run_dir> -t table0 ... -f 'ssprofile_bb_status == "OK"'
```

To *see* the distribution before choosing a bar, ask `table-analyst` to describe
`ssprofile_<label>_radius_of_gyration` against `ssprofile_<label>_n_res`, and the three
fractions, over the `OK` rows. Also ask it to count the `warn:` and `error:` rows.

## What it is blind to

- **Everything outside the selection.** No interface, no clash, no pose, no contact. For
  where a binder landed and which of its residues touch the target, use **`epitope`**
  (`n_iface_design_res`, `iface_design_resnums`). For interface area or shape, use **`cms`**.
  For interface energy, use **`pyrosetta`** (`if_dG`, `if_dSASA`).
- **Energetics, packing and designability.** A compact all-helix chain can still be a bad
  design. For energy and packing, use **`pyrosetta`** (`total_score`, `packstat`). For whether
  the sequence folds to the backbone, run a predictor and then **`usalign`**.
- **Whether a predicted structure is right.** It profiles the coordinates it is given. On a
  prediction, read it after the self-consistency gate (`usalign`), not before.
- **Whether DSSP is right on these backbones.** Neither this mdtraj DSSP nor rpxdock's has
  been cross-checked against an independent SS assignment on this repo's own data (the
  rpxdock skill flags the same open question). On four crystal structures it gave the
  textbook picture, but generated backbones are not crystal structures. Treat the fractions
  as a consistent descriptor, not ground truth.
- **SS within a chain.** Fractions say how much, not where. The per-residue file says where.

## Status

Built 2026-10-09. **Not yet verified on real data.** It needs a `tool-reviewer` pass and a
worker run under a `verification` label before a campaign relies on its columns (see
`verifying-a-tool`). Assumptions that only real outputs can settle:

- the binder sits on the chain letter callers pass;
- generator outputs carry O on every residue (a missing one is a `warn:`);
- generator outputs have no in-chain gaps over 2.5 Å.

Unit tests: `uv run --with mdtraj==1.11.0 python tools/ssprofile/test_ssprofile.py`.
