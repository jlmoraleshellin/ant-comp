---
name: mkcomplex
description: >-
  How do I make a predictor fold the whole binder-target complex instead of the binder alone? Rebuilds a design's multi-chain complex sequence by putting the fixed target chains back around a binder sequence (`--prepend`/`--append`, `--repeat` for homo-oligomeric targets, `--design-copies` for a homo-oligomeric DESIGN chain, `-C` at table width). Load before composing a mkcomplex run or any Boltz/AF3 run on a binder table. Trap: `proteinmpnn_sequence` is the binder monomer alone, and off a proteinmpnn table `-i` is mandatory or nothing is submitted.
---

# mkcomplex

**Custom tool** (`tools/mkcomplex/`, not bundled with prosapia). **`action: update`**: adds a new sequence column to the table it reads, in place; `-t` is required. Pure string work: no structure is read, no GPU, no dependencies beyond the standard library, so the whole per-design step lives in `mkcomplex.sh` (1 CPU, 1 GB; the builder sets `gpus_per_task = 0`, no `-g 0`). On vib `SAPIA_ACTIVATE_MKCOMPLEX` is optional (the script sources it only when set).

## Why it exists: read this before any binder prediction

**ProteinMPNN redesigns only the binder chain, so `proteinmpnn_sequence` holds the binder monomer alone.** Feeding that column to Boltz or AlphaFold3 folds the binder *in isolation*: confident-looking pLDDT for a monomer that was never the question and no interface at all. Nothing errors; the numbers answer a different question. mkcomplex writes a **new** sequence column that puts the fixed target chain(s) back around the design's own sequence in the `/`-separated chainbreak form predictors parse (chains sharing a sequence collapse into one homo-oligomer entity, distinct sequences become separate entities; see `prosapia.utils.chains`).

**Rule: on a binder table, mkcomplex goes between the sequence designer and the predictor, and the predictor's input column is `mkcomplex_sequence`, never `proteinmpnn_sequence`.**

## Invocation

```bash
sapia run mkcomplex <run_dir> --table table1 --prepend-seqs <target sequence> --repeat 2 -C 400
sapia collect mkcomplex <run_dir> --table table1
```

**Off a proteinmpnn table, `-i/--input-column` is mandatory.** The default is `proteinmpnn_sequence`; on any other table that column does not exist, so **zero rows are ready**, the run prints `No designs to submit.` and queues nothing (on Modal `modal-shell --cmd` exits 0 regardless, so it reads as success). Pass the real column (`-i atomium_sequence`, `-i trifuse_sequence`, …) and check the submit line.

| Flag | Default | Notes |
| --- | --- | --- |
| `--input-column` | `proteinmpnn_sequence` | The design's own sequence; nothing else assumes this column |
| `--prepend-seqs` | `""` | Literal fixed chain sequences placed **before** the design's |
| `--prepend-fasta` | `""` | Same, from a FASTA on the runs volume (resolved into the task container's `/runs`) |
| `--append-seqs` / `--append-fasta` | `""` | Same, placed **after** |
| `--repeat` | `1` | N consecutive copies of **each** fixed chain, the way to build a homo-oligomeric target (`--prepend-seqs <protomer> --repeat 2` gives chains A,B). Never touches the design's own chain; with two distinct prepended sequences `--repeat 2` duplicates each consecutively, it does not tile the pair |
| `--design-copies` | `1` | N consecutive copies of **the design's own** sequence (a homo-oligomeric *design* chain, different on every row); see below. `0` or negative is refused at submit |
| `--separator` | `/` | Chainbreak character; `/` is what Boltz and AF3 parse |
| `-C/--max-concurrent` | `40` (prosapia default) | **Raise it to the row count** (`-C 400`): one task is a second of string work, and the default throttle only chops a table into waves of scheduling overhead. At or above the number of ready rows it means a single wave |

**One source per side, never both:** `--prepend-seqs` *or* `--prepend-fasta`, not both; same for append.

**Concurrency.** Cheapest step in a binder chain: do not throttle it. On Modal `-C` becomes `max_containers`, and an over-large value **degrades to queueing, it does not fail** (the plan cap is shared with the whole group: `NO_COLOR=1 modal container list` shows what else is live); on vib it is the array throttle.

## `--design-copies`: repeating the design's OWN chain

`--prepend-seqs`/`--append-seqs` are **one literal for the whole run** and `--repeat` multiplies those fixed chains only, so none can emit the design's own sequence more than once. `--design-copies N` emits N consecutive copies of *this row's* sequence (read per row from `--input-column`):

```bash
sapia run mkcomplex <run_dir> -t <table> -i <fused-sequence column> --prepend-seqs T,T,T --design-copies 3   # T/T/T/F/F/F
```

**The layout is always: prepend chains → the N design copies → append chains**, each `--repeat`/`--design-copies` block consecutive and the flags independent (`--repeat 2 --design-copies 2` with one prepended and one appended chain gives `P/P/F/F/A/A`). At the default (`1`) every existing column keeps its exact name, value and meaning.

`design_chain_index` is the 0-based index of the **FIRST** design copy (`--design-copies 3` with three prepended chains gives `3`; the design occupies chains 3, 4, 5, letters D, E, F). `n_chains`, `chain_lens` and `total_len` all count the copies: `T/T/T/F/F/F` with `len(T)=152`, `len(F)=170` gives `n_chains 6`, `chain_lens 152,152,152,170,170,170`, `total_len 966`, `design_chain_index 3`.

**What N>1 does to the one-chain consumers.** Several downstream tools are built around ONE binder chain; with N identical copies they still answer about **a single protomer**, a different question from the oligomer's:

- `boltzconf` and `af3conf` `--binder-chains` require exactly one chain; their binder pLDDT, `chain(s)_ptm`, interface PAE and ipSAE describe that copy only.
- `propka --binder-chains` cuts the named chains out together as the "free" state: name every protomer of a homo-oligomeric design (`A,B,C`) to compare the oligomer with itself alone; naming one gives a trimer-vs-monomer `dpka`, not the unbinding you probably mean.
- `hisbond` uses one binder chain, so a histidine pointing at a **sister copy** is invisible unless you pass the sister chains (E, F in `T/T/T/F/F/F`) as target chains, which makes the sister interface the thing measured.
- `cms` and `ifacegeom` take chain lists and can cover all N copies, but only if you list them: a list naming one copy silently reports one copy.

Work out the chain letters from `design_chain_index` and `design_copies` before writing any of those flags.

## Traps

- **Chain order decides chain letters downstream.** Prepending puts the target first (chains A, B, …) and the binder last; anything that names chains afterwards (`ringfit --target-chains`, a pose RMSD, a per-chain BSA) must agree. Decide prepend vs append once.
- **Identical chains are collapsed into ONE entity by the predictor.** The N design copies are identical strings by construction, so Boltz and AF3 read `F/F/F` as a single entity with **3 chain ids**, as they read `T/T/T`: you ask for a 3-copy homo-oligomer and the predictor imposes that symmetry. Three *different* chains are three different sequences, which `--design-copies` cannot express. *Signature of the misunderstanding:* `n_chains 6` in the table but two entities in the predictor's input spec or confidence JSON (`chain_pair_iptm` 2x2, not 6x6).
- **`--design-copies` multiplies downstream cost superlinearly**: `T/T/T/F/F/F` at 152+170 is 966 residues, closer to nine times the single-copy job than three. Check `total_len` before launching Boltz/AF3. **Hard ceiling of 26 chains**: chains map positionally onto A–Z downstream (`select_chains` zips onto `string.ascii_uppercase` and truncates silently), so mkcomplex refuses a 27th with an error row.
- **The target sequence must match the target structure you score against.** If the target was trimmed (a membrane belt deleted, say), hand mkcomplex the *trimmed* sequence: a mismatch is the classic off-by-N, the predictor models residues the reference does not have, every later pairing shifts, and the numbers look plausible and are wrong.
- **A collect over a filtered run** reports `Collected N` for every row with the input column and stamps unselected rows `missing`; collect also re-serialises other columns (`3` to `3.0`).
- **A submit-time raise** (e.g. `--design-copies 0`) queues nothing and leaves an empty leaf dir with a `.meta.json`; as an `update` tool it registers no phantom table (only `create` tools leave a `_registry.tsv` row with no `.tsv`).

## Columns collected (`mkcomplex_` prefix)

- `sequence`: the full `/`-joined complex. **The predictor's input column.**
- `n_chains`, `total_len`, `chain_lens`: shape of what was built; all three **count the design copies**.
- `design_chain_index`: 0-based index of the design's own chain; with `--design-copies N` of the **FIRST** copy (the copies are `design_chain_index` … `design_chain_index + N − 1`).
- `design_copies`: the N that run used; read with `design_chain_index` it lets a later reader reconstruct the layout from the table alone.
- `status`: `OK` or `error: …`; every value column is NA on an error row, including `design_copies`.

**`design_copies` is `NA`, never `1`, on rows collected from a run made before the flag existed** (the collector reads it from the per-design TSV and only the current `mkcomplex.sh` writes the column; an invented `1` would be indistinguishable from a recorded one). So `NA` means "old run *or* failed row": disambiguate with `mkcomplex_status`.

Check `n_chains`, `chain_lens` and `design_copies` on the first row before an expensive prediction. Filters:

```bash
-f 'mkcomplex_design_copies == 3 and mkcomplex_status == "OK"'                 # rows built as a 3-copy design oligomer
-f 'mkcomplex_n_chains == 6 and mkcomplex_design_chain_index == 3'              # the exact T/T/T/F/F/F shape
-f 'mkcomplex_status == "OK" and mkcomplex_total_len < 1200'                    # inside the predictor budget
```

## Error rows, not plausible strings

Each is an `error: …` status with **every value column empty (`NA`, never `0`)**; the task still exits 0, so the rest of the array collects.

| Input | Result |
| --- | --- |
| `--design-copies 0` or negative | Refused **at submit**; nothing queued |
| Input sequence containing the separator | `error: input sequence already contains the separator` (it is already a multi-chain string) |
| Non-uppercase / non-standard residue in a **fixed** chain | `error: invalid residue in fixed chain …` |
| Non-uppercase / non-standard residue in the **design's own** sequence | `error: invalid residue in design sequence …` |
| Empty / `NA` / `None` input sequence | `error: empty input sequence` |
| An **empty element** of a non-empty fixed-chain list (interior `T,,T`, or a FASTA record with a header and no sequence) | `error: empty chain in the fixed-chain list`; a *trailing* comma (`T,`) is harmless, and giving a side no chains at all is fine |
| More than 26 chains in total | `error: … exceeds the 26-chain ceiling` |

The design sequence is held to the same 20-uppercase-code alphabet (`ACDEFGHIKLMNPQRSTVWY`) as the fixed chains; older versions validated only the fixed chains and would copy a lowercase or `X` residue through (replicated N times by `--design-copies`). A previously-`OK` row containing `X` or lowercase now errors. **Known `X` source:** `trifuse_worker.one_letter()` can emit `X` for a residue gemmi does not tabulate, so a `trifuse_sequence` with `X` becomes an error row: count non-`ACDEFGHIKLMNPQRSTVWY` characters in the input column before a run rather than discovering a wave of error statuses.
