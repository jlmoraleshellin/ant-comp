---
name: proteinmpnn
description: >-
  ProteinMPNN sequence design on Modal. Load before composing a proteinmpnn run. Traps: from an rfdiffusion3 table pass `-i rfdiffusion3_path` or nothing is submitted; a bare flag in `--set` needs `--set=`; `--fixed-positions` requires a chain list; fixing a whole target chain can substitute for mkcomplex.
---

# proteinmpnn

Designs sequences for input backbones. **`action: create`** — mints a child table, one row per designed sequence.

Full reference: `docs/tools/proteinmpnn.md` in the prosapia repo, and
`sapia run proteinmpnn --help`.

IMPORTANT: prosapia's bundled proteinmpnn is a wrapper of the original. Check the github repo for all information: https://github.com/dauparas/ProteinMPNN/tree/main

## Verified invocation

```bash
sapia run proteinmpnn <run_dir> -t table0 -i rfdiffusion3_path --num-seq-per-target 2
```

5 backbones → 10 sequences in `table1`, one task on an L4, ~1.5 min including a first-time
image build.

## The input-column trap

`default_input_column` is **`rfdiffusion_path`** — RFdiffusion, *not* rfdiffusion3. Coming
from an rfd3 table you **must** pass `-i rfdiffusion3_path`, or the run finds no ready rows
and submits nothing. Check the parent table's columns before composing the command.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--num-seq-per-target` | — | Sequences sampled per input backbone. |
| `--sampling-temp` | — | Lower = more conservative. |
| `--designs-per-task` | 10 | Designs bin-packed per task. Designs with identical params share one batched call. Lower for more parallelism. |
| `--chains-to-design` | all | Chain mini-language: `A:C,E` → `A B C E`. |
| `--fixed-positions` | none | Keep positions fixed. `/` separates chains, `,` fragments, `start:end` inclusive 1-indexed, `{expr}` resolves up the lineage. E.g. `9:23/10,11,18:20`. **Requires `--chains-to-design` or `--symmetry`** — see below. |
| `--tied-positions` | none | Tie positions across chains. Same language. Mutually exclusive with `--symmetry`. |
| `--symmetry` | none | Homo-oligomer convenience: `auto` or an integer chain count. Ties all chains, and lets `--fixed-positions` describe a single asymmetric unit. **A plain number, never a point group** (`12`, not `C12`). |
| `--bias-aa` | none | `'D:1.39 E:1.39 K:1.39'`. |
| `--set` | none | Raw flag passed through, repeatable. **A bare flag needs `=`:** `--set=--use_soluble_model`. See the quoting trap below. |
| `--seed`, `--batch-size` | — | |

Default Modal resources: **L4**, 8 CPU, 8 GiB, 1 h timeout. Weights ship inside the image,
so there is no cache Volume and no download step.

## The `--set` quoting trap

`--set '--use_soluble_model'` **fails** with `--set: expected one argument`. The submitter's argparse treats a quoted value that starts with `-` as a new option, so the flag is consumed and `--set` is left without an argument. Nothing is queued.

Rule: **single-token flag → `--set=`; flag with a space-separated value → either form.**

`--set` is `action="append"`, so it is **repeatable — stack as many flags as you need**, each
in whichever form that flag requires:

```bash
sapia run proteinmpnn <run_dir> -t table0 -i rfdiffusion3_path \
    --set=--use_soluble_model \
    --set '--omit_AAs C' \
    --set=--bias_by_res_jsonl=/runs/inputs/bias.jsonl
```

All `--set` values are joined with spaces onto the `protein_mpnn_run.py` argv tail. Two flags
in **one** `--set` also work (`--set '--ca_only --use_soluble_model'`) — the space is what
keeps argparse from mistaking it for an option. Corollary of that same joining: the tail is
expanded **unquoted** in `proteinmpnn.sh`, which is how `'--omit_AAs C'` becomes two argv
entries — so a path containing a space will break. Keep Volume paths space-free.

Per-position bias goes through `--set` (`--bias_by_res_jsonl`); the typed `--bias-aa` is
global-only.

## `--fixed-positions` needs a chain list

The submitter refuses `--fixed-positions` / `--tied-positions` unless `--chains-to-design`
or `--symmetry` is also given, because the position groups map **one-to-one onto that chain
list, in order**.

This is a guardrail against a silent no-op, not an arbitrary limit. ProteinMPNN's
`make_fixed_positions_dict.py` looks each chain up by its index in `--chain_list`; with no
chain list nothing matches, the fixed dict comes back empty for *every* chain, and the whole
structure is redesigned at exit 0. Group count matters too: one group too few is an
`IndexError` (loud), but a group of the wrong *length* silently frees the residues it misses.

### Designing all chains with the target fixed (instead of mkcomplex)

The FASTA holds **only the designed chains**, `/`-joined — that is why
`proteinmpnn_sequence` is the binder monomer alone and why `mkcomplex` exists. Making the
target a *designed* chain with every position fixed should bring it back into the sequence:

```bash
--chains-to-design A,B --fixed-positions 1:609/      # A = target (609 aa), B = binder free
```

The trailing empty group leaves the binder designable. Inside the model this is a
near-identity with the usual `--chains-to-design B`: `sample()` computes
`chain_mask = chain_mask * chain_M_pos * mask`, so an all-fixed designed chain *is* a visible
chain — same conditioning on the true target sequence, same fixed-first decoding order, same
cost — and `mask_for_loss` still excludes fixed positions, so `score` and `seq_recovery` stay
binder-only and stay comparable to tables built the other way.

**Unverified in this workspace.** Read from upstream source, never run here. Test on one
backbone and one sequence, checking for a target+binder-length record, before trusting it.

**Prefer `mkcomplex` anyway**, unless you actually need part of the target redesigned:

- It keeps the folded target sequence **independent of the design coordinates** — swapping a
  crystallographic construct variant for wild-type, restoring a trimmed region, or
  `--repeat`-ing an oligomer is a flag change, not a re-run of sequence design.
- It preserves the categorical check: `fixed_chains=['A'], designed_chains=['B']` in every
  FASTA header and no target-length record anywhere. Under fixed positions every record is
  target+binder length and "was the target redesigned?" becomes a string diff you must run
  yourself.
- It saves nothing but one submit+collect of pure string work — no GPU.

## What it collects

Child rows named `<parent>_f1`, `<parent>_f2`, … each linked to its parent backbone.

Columns (leaf-prefixed `proteinmpnn_`): `sequence`, `score` (lower is better),
`global_score`, `seq_recovery`, `T` (sampling temp), `sample`, `iteration`, plus `_status`
and `_path` (the FASTA).

`proteinmpnn_sequence` is what the structure predictors consume — it is Boltz's default
input column, so a Boltz run straight after needs no `-i`.
