---
name: mkcomplex
description: >-
  Rebuilds a design's full multi-chain complex sequence by putting the fixed target chains back around a ProteinMPNN binder sequence, so a predictor folds the complex instead of the binder alone. Load before composing a mkcomplex run, or before any Boltz/AF3 run on a binder table. Covers `--prepend`/`--append`, `--repeat` for homo-oligomeric targets, and launching at full table width (`-C 400`).
---

# mkcomplex

**Custom tool** (lives in `tools/mkcomplex/`, not bundled with prosapia).
**`action: update`** — adds a new sequence column to the table it reads, in place.
`-t` is required. Pure string work: no structure is read, no GPU, no dependencies
beyond the standard library, so the whole per-design step lives in `mkcomplex.sh`.

## Why it exists — read this before any binder prediction

**ProteinMPNN redesigns only the binder chain, so `proteinmpnn_sequence` holds the
binder monomer alone.** Feeding that column straight to Boltz or AlphaFold3 folds the
binder *in isolation* — you get confident-looking pLDDT for a monomer that was never
the question, and no interface at all. This is a silent failure: nothing errors, the
numbers just answer a different question.

mkcomplex writes a **new** sequence column that puts the fixed target chain(s) back
around the design's own sequence, in the `/`-separated chainbreak form the predictors
already parse. Chains sharing a sequence collapse into one homo-oligomer entity;
distinct sequences become separate entities (see `prosapia.utils.chains`).

**Rule: on a binder table, mkcomplex goes between proteinmpnn and boltz, and the
predictor's input column is `mkcomplex_sequence`, never `proteinmpnn_sequence`.**

## Invocation

```bash
sapia run mkcomplex <run_dir> \
    --table table1 \
    --prepend-seqs MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ \
    --repeat 2 \
    -C 400
sapia collect mkcomplex <run_dir> --table table1
```

| Flag | Default | Notes |
| --- | --- | --- |
| `--input-column` | `proteinmpnn_sequence` | The design's own sequence. Nothing assumes this column. |
| `--prepend-seqs` | `""` | Literal fixed chain sequences placed **before** the design's. |
| `--prepend-fasta` | `""` | Same, from a FASTA **on the runs volume** (resolved into the task container's `/runs`). |
| `--append-seqs` / `--append-fasta` | `""` | Same, placed **after**. |
| `--repeat` | `1` | Emit N consecutive copies of **each** fixed chain — this is how you build a homo-oligomeric target. `--prepend-seqs <protomer> --repeat 2` gives chains A,B. |
| `--separator` | `/` | Chainbreak character. `/` is what Boltz and AF3 parse. |

**One source per side, never both:** `--prepend-seqs` *or* `--prepend-fasta`, not both.
Same for append.

## Concurrency: launch the whole table at once

**One task per design, and a task is a second of string work** — 1 CPU, 1 GB, no GPU (the manifest builder sets `gpus_per_task = 0` itself, so you never pass `-g 0` here), no weights, no image build worth the name. The only thing that makes mkcomplex slow is the default throttle: `-C/--max-concurrent` is `40`, which chops a 400-row table into ten waves of scheduling overhead for work that would otherwise finish in one.

**Raise it to the row count.** `-C 400` is fine on **both** backends, and on Modal there is room above that. Size `-C` to the table rather than guessing: at or above the number of ready rows it means a single wave.

**What actually caps it on Modal.** `-C` becomes `max_containers` on the Function (`concurrency_limit` server-side) with no client-side validation, so three ceilings apply, in this order:

| Ceiling | Value | What happens at it |
| --- | --- | --- |
| Per-Function hard limit | **4,000** concurrent containers | Modal's own limit for a single Function; `-C` above it is pointless. |
| Workspace plan cap | **100** containers (Starter) / **5,000** (Team) / custom (Enterprise) | Shared by the whole workspace. |
| `spawn_map` enqueue rate | inputs sent 512 per call | Not a cap: on `RESOURCE_EXHAUSTED` the client retries with a warning about "rate limits or function backlog limits". |

**An over-large `-C` degrades to queueing, it does not fail.** Tasks past the cap simply wait for a slot, so the cost of guessing high is zero — which is why `-C 1000` is a reasonable default on a Team-plan workspace and merely a no-op above 100 on Starter.

**The Modal workspace is shared with the lab.** `vubmodal` has several members, so the plan cap is a *group* entitlement, exactly like `--max-gpu-fraction` on vib. mkcomplex containers are 1 CPU and live seconds, so 400–1000 of them is not the problem; the judgement call is a long-running GPU tool, not this one. Check what else is live with `NO_COLOR=1 modal container list` before claiming a four-figure slice.

```bash
sapia run mkcomplex <run_dir> -t table1 -C 400 --prepend-seqs <target> --repeat 2
```

This is the cheapest step in any binder chain. Do not throttle it to be polite — the expensive neighbours (boltz, bindcraft2, rfd3) are where concurrency limits earn their keep, and holding mkcomplex back just delays them.

## Traps

**Chain order decides chain letters downstream.** Prepending puts the target first, so
the target becomes chains A, B… and the binder is the last chain. Anything you run
afterwards that names chains — `ringfit --target-chains`, a pose RMSD, a per-chain BSA —
must agree with that order. Decide prepend vs. append once and keep it.

**`--repeat` copies each fixed chain, not the whole set.** With two distinct prepended
sequences and `--repeat 2` you get each duplicated consecutively, not the pair tiled.

**The target sequence must match the target structure you scored against.** If you
trimmed the target (deleting a membrane belt, say) the sequence handed to mkcomplex is
the *trimmed* one. A mismatch here is the classic off-by-N: the predictor models
residues the reference does not have, every later residue pairing shifts, and the
numbers look plausible and are wrong.

## Columns collected (`mkcomplex_` prefix)

- `sequence` — the full `/`-joined complex. **This is the predictor's input column.**
- `n_chains`, `total_len`, `chain_lens` — shape of what was built.
- `design_chain_index` — 0-based index of the design's own chain among the chains.
  Use it to work out which chain letter the binder is.
- `mkcomplex_status`.

Check `n_chains` and `chain_lens` on the first row before launching an expensive
prediction — they are the cheapest possible confirmation that the complex is what you
think it is.
