# Status, and the gap: LPC is not yet a prosapia tool

**Date:** 2026-10-07

## What exists right now

### In `ant-comp` (this repo)

```
campaigns/tnf/
├── 01_recon_human_mouse.md       conservation analysis, epitope rationale
├── 02_lpc_inputs.md              what LPC needs, field by field
├── 03_status_and_gap.md          this file
├── targets/
│   ├── 1tnf_human.pdb            cleaned trimer, chains A/B/C, 6-157
│   ├── 2tnf_mouse.pdb            cleaned trimer, chains A/B/C, 9-157
│   └── targets_entry.yaml        the config block, with rationale
└── scripts/                      7 rerunnable analysis scripts
```

`tools/` is **unchanged** — still atomium, bindcraft2, chainsel, cms, mkcomplex,
ringfit. There is no `tools/laproteina/`.

### Outside the repo

`/home/lpignatti/Proteina-Complexa` — LPC cloned (branch `dev`), `--minimal`
environment built (7.3 GB venv), `.env` and `env.sh` configured, both target
entries registered and listed by `complexa target list`.

Validation state:

| check | result |
| --- | --- |
| `.env`, `DATA_PATH`, config parse | pass |
| target entries (both) | pass |
| `complexa.ckpt` / `complexa_ae.ckpt` | **absent** (6.5 GB, not downloaded) |
| Foldseek, `sc` | absent, optional, evaluation-only |

## The gap

The campaign plan is:

```
LPC backbones  ->  atomium (tied positions, human+mouse)  ->  Boltz
```

**atomium is a prosapia tool.** It reads a structure-path column out of a
prosapia table and writes a child table of sequences. LPC, as configured, is a
standalone program that writes PDB files into its own `./inference` directory.
It knows nothing about `run_dir`s, tables or lineage.

So the two halves do not currently connect. Something has to put LPC's backbones
into a prosapia table before atomium can design on them.

## Two ways to bridge it

### (a) Write the tool — `tools/laproteina/`

Five files, per the `authoring-a-tool` skill:

| file | role |
| --- | --- |
| `spec.py` | the declarative `Tool` descriptor: name, `action: create`, default input column, resources |
| `run_laproteina.py` | builds the task manifest — one task per generation job |
| `collect_laproteina.py` | folds the output PDBs back into a child table |
| `laproteina.sh` | the per-task script that invokes `complexa` |
| `modal_image.py` | the Modal image: LPC env + checkpoints |

`action: create`, because it produces new entities (backbones), so it mints a
child table and links each row to its parent.

Pros: LPC becomes a first-class step. Lineage, reruns, `--force`, status columns
and the `.exit` loop all work. atomium reads its output column directly. This is
what the repo is built for, and `tool-creator` + `authoring-a-tool` exist to do
exactly this.

Cons: it is real work, and it is the piece most likely to need several iterations
against Modal.

### (b) Run LPC standalone, import the backbones

Run `complexa design` by hand (on Modal, in a plain container), then write a
small script that drops the resulting PDBs into a prosapia `run_dir` as a table
with a `laproteina_path` column, so atomium can pick them up.

Pros: much less code; gets to the science sooner.
Cons: no lineage, no rerun safety, no status columns. A one-off bridge that
someone has to maintain, and that will quietly rot.

## Recommendation

**(a), the real tool** — but only after a standalone smoke run proves LPC
generates sensible backbones against target `90_TNFa_HUMAN_xreact`.

Writing a tool wrapper around a program that has never produced a single correct
output is the wrong order. The smoke run tells us the flags, the output layout
and the per-design cost, and all three are inputs the tool wrapper needs anyway.

This is the same staging already agreed: cheap control first, then the real
thing.

## Immediate next steps

1. Decide where the checkpoints live (local vs. straight into a Modal Volume).
   The imec network intercepts TLS to Modal blob storage, so `modal volume put`
   from here is unreliable — pulling them inside a Modal container avoids that.
2. Configure Modal (`modal setup`   — browser login, user action).
3. Smoke run: a handful of designs against `90_TNFa_HUMAN_xreact`, search off,
   one GPU. Goal is a measured cost-per-design and a known output layout, not
   science.
4. Then write `tools/laproteina/`, informed by what step 3 showed.
5. Then the tied-positions atomium step, which is where cross-reactivity is
   actually designed in.

## Open items carried forward

- [ ] Checkpoints: local download vs. direct into a Modal Volume
- [ ] `modal setup` (user action)
- [ ] Confirm access to the Modal workspace holding the `github-token` /
      `github-username` Secrets — atomium's image build needs them, and atomium
      is the cross-reactivity step
- [ ] Smoke run on Modal
- [ ] `tools/laproteina/`
- [ ] Decide how the two-complex structure for tied positions gets built
      (nothing in `tools/` does this yet)
