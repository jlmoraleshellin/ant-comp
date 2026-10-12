---
name: boltz
description: >-
  Boltz structure prediction on Modal: flags, weight cache, timing and confidence columns. Load before composing a boltz run, and before reading any confidence number from a binder complex. Traps: `iptm` is diluted by chain count (the binder-specific numbers are `pair_chains_iptm`, `chains_ptm` and PAE), and forced templates fail silently.
---

# boltz

Predicts structures for designed sequences. **`action: update`** — annotates the table it
reads, in place. It does **not** create a table, and `-t` is required.

There is no dedicated doc page; `sapia run boltz --help` is authoritative.

IMPORTANT: prosapia's bundled boltz is a wrapper of the original. Check the github repo for all information: https://github.com/jwohlwend/boltz/tree/main/docs

## Verified invocation

```bash
sapia run boltz <run_dir> -t table1
```

10 sequences predicted in ~16 min on one task, including the first-time weight download.
Later runs reuse the cache and start much faster.

`default_input_column` is `proteinmpnn_sequence`, so straight after a ProteinMPNN step no
`-i` is needed. From any other source, pass `-i <sequence column>`.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--shard-size` | 10 | YAML inputs per shard directory / task. |
| `--devices` | 1 | GPUs per task; also sets `--gpus-per-task` to match. |
| `--use-msa-server` | off | Adds MSA information **to every chain**. Slower, and calls an external server — don't enable it without being asked. |
| `--msa-empty-chains` | unset | Chain mini-language (`B`, `A,C`, `B:D`). Chains pinned to `msa: empty` **even when `--use-msa-server` is on** — name the de-novo binder chain. Requires `--use-msa-server`; an unknown chain is an error. See the per-chain MSA section. |
| `--chains` | all | Chain mini-language, e.g. `A:D`. Letters beyond the sequence's chain count are dropped. |
| `--positions` | full | Position mini-language: `/` maps groups onto chains, `start:end` inclusive 1-indexed, open ends allowed (`10:`, `:50`), `{expr}` islands resolve up the lineage. |
| `--template-yaml` | none | File spliced verbatim as a `templates:` block into every input YAML. |

Default Modal resources: **A10**, 24 CPU, 64 GiB, 8 h timeout.

## The weight cache

Weights and the CCD download on first use into the Volume named by `SAPIA_MODAL_VOLUME_BOLTZ_CACHE`, mounted at `/boltz_cache`. It holds `boltz2_conf.ckpt`, `boltz2_aff.ckpt`, `mols/` and `mols.tar`. Budget extra time on the very first run of a fresh cache; never delete that Volume to "clean up".

## What it collects

Columns added to the **same** table (leaf-prefixed `boltz_`): `confidence_score`, `ptm`, `iptm`, `ligand_iptm`, `protein_iptm`, `complex_plddt`, `complex_iplddt`, `complex_pde`, `complex_ipde`, plus `_status` and `_path` (the predicted structure).

Reading them: `confidence_score` and `complex_plddt` around 0.9+ is a confident prediction for a de-novo monomer; `ptm` tracks global fold confidence. High confidence means the predictor believes the fold — it is **not** proof the sequence folds to the backbone it was designed for. For that, compare the prediction back to its parent backbone with `usalign`.

##  Important notes

- Forced templates steer; they don't constrain fully (e.g.: `force: true, threshold: 2.0`). Use the flag anyway but always verify independently and exclude rows where the target didn't land.

## Binder complexes: what the collected confidence does and does not mean

**`iptm` is diluted by chain count.** With N protein chains it averages N(N−1)/2 interchain
pairs. In a binder complex all but a handful of those are target–target — and if you supplied a
template, they are pairs *you forced*. A 9-chain complex has 36 pairs of which **35 are
target–target and 1 is binder–target**, so the collected `iptm` is overwhelmingly a statement
about template reproduction.

Measured on a 27-design binder campaign against a trimmed 4-protomer target:

| | |
| --- | --- |
| `confidence_score` vs **target** RMSD | **r = −0.896** — a perfect separator (landed 0.706–0.736, distorted 0.542–0.634, empty band) |
| `confidence_score` vs **binder-interface** ipTM, among rows that landed | **rho = −0.250** |

So: **use `confidence_score` as a gate for "did the target assemble", then stop using it.** It
carries no binder signal once past that gate. The highest-`confidence_score` design in that set
also had the highest target pLDDT and near-bottom binder-interface confidence.

### The binder-specific numbers, and where to get them

Not collected by prosapia. They are in the per-prediction confidence JSON at
`<out_dir>/boltz_results_shard_*/predictions/<name>/confidence_<name>_model_0.json`:

- `chains_ptm` — `{"0": …, "1": …}`, **chain index is positional from the input YAML**, so the
  binder is the last index if `mkcomplex` prepended the target.
- `pair_chains_iptm` — full N×N. **Not symmetric** (observed max asymmetry 0.079); symmetrise
  by averaging both directions. The diagonal just repeats `chains_ptm`.

Also on disk per prediction: `plddt_<name>_model_0.npz` (`plddt`, 0–1, per token) and
`pae_<name>_model_0.npz` (`pae`, n_res × n_res). The cif's `B_iso_or_equiv` carries the same
pLDDT on a 0–100 scale (agrees to 7e-5). There is **no PAE or per-residue data in the JSON**.

Caveat: `pair_chains_iptm` is whole-chain vs whole-chain, so a binder touching a small patch of
a large chain is still diluted. An **interface-restricted PAE** (contacting residue pairs only,
from the npz) is the sharper metric and is not provided by anything here.

### The cheapest honest interface signal

**Compare the binder's mean pLDDT alone vs. in complex.** A real interface raises it. Measured:
all 27 designs *lost* 7.5–23.3 points (mean 13.7) on docking — the signature of a
confidently-folded domain the predictor is confidently unsure where to place. This was more
informative than any ipTM flavour.

Rough reference points: a believed interface is ~0.5–0.6 binder-target ipTM and < 10 Å interface
PAE. The best in that campaign were 0.357 and 19.5 Å.

## MSA policy differs per chain in a binder complex

`--use-msa-server` is off by default and the early runs here used `msa: empty` everywhere. For a
**de-novo binder that is correct** — no MSA exists, a search is meaningless, and it would put an
unpublished design sequence on an external server. For a **natural target it is not**: denying
it an MSA weakens target assembly *and* interface confidence, and the result then reads as bad
designs. Boltz supports the split (an entity with `msa: empty` is skipped, an entity without the
key is searched), and **`--msa-empty-chains` is how you ask for it**:

```bash
# 2-chain complex: A the natural target (searched), B the de-novo binder (no MSA)
sapia run boltz <run_dir> -t table1 -i mkcomplex_sequence \
    --use-msa-server --msa-empty-chains B
```

This is the **only configuration in which the "did the target land" gate means anything** — a
target folded without an MSA may fail that gate for reasons that have nothing to do with the
binder. Say which policy you used when reporting; it is recorded in the run's
`<out_dir>/.meta.json` (`use_msa_server`, `msa_empty_chains`) and visible per design in
`<out_dir>/boltz_inputs/<name>.yml`.

Its contract, all enforced at **submit** time (nothing is queued if it raises):

| Situation | Behaviour |
| --- | --- |
| Flag unset | **Byte-identical to bundled boltz**, with or without `--use-msa-server`. Regression-tested. |
| `--msa-empty-chains` without `--use-msa-server` | **Error.** It cannot be honoured (nothing is searched anyway), and accepting it would hide a forgotten `--use-msa-server`. |
| A named chain the complex does not have | **Error**, naming the chains that *are* present. Never a silent no-op — that would mean the binder quietly got an MSA. |
| A named chain whose sequence is shared with an unnamed chain | **Error.** Identical sequences collapse into one boltz entity, which has one MSA policy. Name all of them or none. |
| A named chain that `--chains` dropped | **Error** (it is no longer in the complex). |
| Every chain named | Same YAML as plain `msa: empty` everywhere. |

Chain letters are **positional in the `/`-chainbreak sequence** (first segment is `A`), *after*
`--chains` has narrowed it — not the chain names in any PDB you started from. After `mkcomplex`
prepended a target, the binder is the last letter.

What it does **not** do: supply a *pre-computed* MSA (boltz's `msa: <path.a3m>` form) — the
choice is only searched vs empty; and it cannot tell you whether a searched chain actually found
homologues. A natural target with no hits behaves like an empty one and no column says so.

### Where this flag lives

It is **not** in prosapia's bundled boltz. `tools/boltz/` in this repo is a user-side **shadow**
(`get_builtin("boltz").with_overrides(...)`, replacing only `add_run_args_fn` and
`build_manifest_fn`); the task script, Modal image and collector are still the bundled ones, and
the tool keeps the name `boltz`, so every column and every downstream step is unchanged.
Consequence: like every other custom tool here, **it is only present when `sapia` is run with
this repo's `tools/` on `$PROSAPIA_TOOLS_DIR`** — i.e. `sapia modal-shell` launched from the
project root, or the synced vib workspace. Launch the workstation from elsewhere and boltz
silently reverts to the bundled tool, where `--msa-empty-chains` is simply an unrecognised flag
(argparse errors — it does not run with the flag ignored). Self-test:
`uv run python tools/boltz/test_boltz.py`.

## Templates

**Do not hand-write the CIF.** `utils/boltz_template.py` builds it and the matching YAML,
and verifies every trap below by re-reading the written file. It is not a sapia tool — run it with
the workstation's python:

```bash
sapia modal-shell --cmd 'python /runs/utils/boltz_template.py \
    --fetch 7OJG --out /runs/inputs/7ojg_tmpl_kabc.cif \
    --chains K,A,B,C --rename-to A,B,C,D --residues K:19-59+107-155 \
    --yaml /runs/inputs/7ojg_tmpl_kabc.yaml'
```

`--chains` picks source chains in order, `--rename-to` renames positionally, `--residues` trims
per **source** chain (`K:19-59+107-155`), `--predict-chains` sets the YAML's `chain_id` when the
prediction's chain IDs differ from the template's. `--check <file> --expect-chains A,B,C` verifies
a CIF from anywhere and exits non-zero. Outputs go under `/runs/` — each `modal-shell` call is a
fresh container. `--fetch` downloads from RCSB inside the container, which is how inputs get
staged here (`modal volume put` is unreliable on this network).

**The script lives at `/runs/utils/boltz_template.py` in the workstation.** It is NOT under
`$PROSAPIA_TOOLS_DIR/_utils/` — that path appeared here until 2026-09-30 and fails with `rc=2`.

### A source chain with UNRESOLVED residues: `--keep-gaps`

**Default mode compacts the declared sequence, and that is wrong for any crystal structure with
an unbuilt loop.** Without the flag the template declares exactly what is MODELLED and numbers it
1..N, so if the source is missing residues 6–8 and a loop at 103–110, position *i* of the
template is no longer position *i* of the query and **every residue after the first gap is
shifted** by the size of that gap. Nothing errors; the prediction is just forced onto the wrong
correspondence. The only sign in the default mode is a `note:` on stderr — *"143 residues for a
range spanning 152"* — which is easy to read as harmless.

`--keep-gaps` (off by default; without it the output is unchanged — a regression suite rebuilds
seven legacy flag combinations with a pre-change copy of the script and compares the md5 of **both**
the CIF and the YAML, opt-in via `BOLTZ_TEMPLATE_REFERENCE`) fixes the mapping for **one
`--residues <Src>:a-b` range per chain**:

- `_entity_poly_seq` declares the **whole** range `a..b`, taking the names of unresolved positions
  from the source's own `_pdbx_poly_seq_scheme`, matched by **author number** (the source entity
  is usually longer than the range, so matching by position in the entity would be wrong);
- `_atom_site.label_seq_id` of a modelled residue is `author number - a + 1`, so it has **gaps**;
- an unresolved position gets **no atoms**. Coordinates are never invented — a filled gap would be
  restrained by the `force`/`threshold` potential like a real one;
- the YAML is unchanged, and each chain is still its own entity (no merging).

Anything it cannot do honestly is an error, not a shorter template: a range running past the
declared sequence, more than one range for a chain, `--renumber` or `--keep-het` at the same time,
a source with no `_pdbx_poly_seq_scheme`, microheterogeneity (two scheme rows at one `seq_id`)
inside the range, a modelled residue whose name disagrees with the declared name at its position,
or an unresolved residue with no author number whose resolved neighbours disagree about the offset.
Only the chains you ask for are read, so an oddity in a chain you do not template cannot block you.

**A failed build writes nothing.** The CIF and YAML go to `<out>.tmp` / `<yaml>.tmp` and are
renamed into place only once the checks pass; on failure both are deleted. So **an absent output
file is the failure signal** — which is the only one that survives `sapia modal-shell --cmd`,
since that always exits 0. Check the file exists before using it; do not trust the printed text.

```bash
sapia modal-shell --cmd 'python /runs/utils/boltz_template.py \
    --in /runs/inputs/8Z8M_fetched.cif --out /runs/inputs/8z8m_tnf_ace.cif \
    --chains A,C,E --rename-to A,B,C --residues A:6-157,C:6-157,E:6-157 \
    --keep-gaps --expect-sequence <152-letter query> --allow-mismatches 1 \
    --yaml /runs/inputs/8z8m_tnf_ace.yaml'
```

The script in the workstation is whatever copy was last staged on the **Volume** at
`/runs/utils/boltz_template.py`; an edit in the repo does not reach it until it is staged, and a
campaign may stage its own copy inside its run_dir. Check the md5 of the copy you are calling.

**Why gaps are safe but names are not** (read from boltz `b1ebfc46`): `parse_polymer` does **not**
use `label_seq_id` — it aligns the modelled residue **names** to `_entity_poly_seq` with
`gemmi.align_sequence_to_polymer`. A declared residue with no atoms becomes a token with
`is_present=False`, zeroed and unrestrained. But a modelled residue whose name differs from the
declared name at its place is treated as absent *without advancing the polymer index*, desyncing
everything after it, and a modelled residue missing from `_entity_poly_seq` can overrun
(boltz issue 669, `IndexError` at `mmcif.py:588`). Both are refused at build time and by `--check`.

**The query→template map is a gapless diagonal.** With `chain_id` and `template_id` both given it
is a local BLOSUM62 alignment with gap penalty −1000 against the **declared** sequence, then
positional, with **no identity threshold** — a wrong-length or shifted declaration is used anyway
and silently. `--expect-sequence <query>` is the check for that: the query is upper-cased and
stripped of whitespace, then compared position by position with the declared sequence. A length
difference **fails**, and so does any number of differing positions above **`--allow-mismatches N`
(default 0)** — a template shifted by k differs almost everywhere, so it fails rather than printing
a note. A genuine point difference is legitimate and stays mapped, so allow exactly as many as you
have read and accepted: **8Z8M's TNF needs `--allow-mismatches 1`**, for a real **L143D** — that is
declared **position 138** of the 6–157 range, i.e. **author 143** (both numbers are printed in
every message; 138 is not a shift).

**The template is also put through the alignment Boltz itself runs.** `verify()` calls
`gemmi.align_sequence_to_polymer` with the default scoring, exactly as `parse_polymer` does, walks
the match string, and **fails unless every modelled residue lands on the declared position equal to
its `label_seq_id`**. This is not a formality: a low-complexity chain (measured on a synthetic
poly-alanine chain with a 12-residue gap) is positionally correct by construction yet gemmi scores
a *gapless* placement higher and puts the C-terminal residues inside the gap — `8M12I`, 5 of 8
residues moved. Such a template is refused, and nothing is written.

### Checking a template

`--check <file> --expect-chains A,B,C` verifies a CIF from anywhere and exits non-zero. A clean
contiguous file still prints exactly the three familiar `OK` lines (`label_asym_id ==
auth_asym_id`, `_entity_poly_seq` count vs modelled, `template chains present`), and a `.` in
`label_seq_id` is now a `FAIL` of its own rather than something to check by hand.

For a `--keep-gaps` file add **`--expect-declared N`** (the range length), or **`A:152,B:140`** when
the chains differ. It then verifies per entity that the declaration is exactly that, that
`label_seq_id` is unique, ascending and inside `1..N`, that no more residues are modelled than
declared, and that each modelled residue's name equals the declared name at its `label_seq_id` —
and it **prints the number and the positions of the unresolved residues per chain, with the
modelled fraction**:

```
OK  chain A: _entity_poly_seq = 152 residues, 143 modelled at label_seq_id 4..152 (keep-gaps: declared complete, atoms only where resolved)
    chain A: 9 unresolved position(s): 1-3,101-106; 143/152 modelled (94%)
OK  chain A: Boltz's alignment places all 143 modelled residues at their label_seq_id
```

A build passes the declared length **per chain**, so chains with ranges of different lengths are
each checked against their own. Without `--expect-declared` the checker recognises a gapped chain
from its own numbering, but **only when the numbering actually has a gap**: a file that declares
more than it models while numbering 1..M contiguously is the old silent-shift failure and still
**FAILS**, so pass `--expect-declared` when a legitimate keep-gaps file happens to have only its
C-terminal tail unresolved. The flip side, on a **foreign** CIF: `--check` now also *passes* a file
whose surplus declarations sit exactly where `label_seq_id` has holes and whose names match there —
if you wanted the old blanket refusal of any over-declaration, that case no longer gives it to you.

**`--expect-modelled-min N`** fails a chain modelling fewer than `N` residues; without it, a
template with 2 of 152 residues modelled passes quietly (the fraction line is the only hint).

A modified residue (MSE and friends) is reported by name rather than normalised, because it
compares as its parent letter in any sequence comparison.

The rest of this section is why that script exists — read it before overriding any of its choices.

`--template-yaml` is spliced **verbatim** at column 0 after the `sequences:` block, so the file
must start at column 0. Schema read from `boltz/data/parse/schema.py`:

```yaml
templates:
  - cif: /runs/inputs/target.cif
    chain_id: [A, B, C, D]      # prediction chains
    template_id: [A, B, C, D]   # template chains
    force: true
    threshold: 2.0              # mandatory when force is true
```

- **Give `chain_id` AND `template_id`, equal length.** Then Boltz zips them positionally — an
  explicit 1:1 map. Omit either and it falls back to `get_template_records_from_search`, which
  scores sequence alignments to pick chains. **With sequence-identical chains (any homo-oligomer)
  that can silently assign the wrong template chain and destroy the geometry.**
- **Boltz names template chains by `label_asym_id`/subchain, not auth chain ID.** A CIF written
  by gemmi defaults to `Axp`, `Bxp`, … and `template_id: [A, …]` then raises *"Template chain A
  is not one of the protein chains"*. Set `label_asym_id == auth_asym_id`.
- **`_entity_poly_seq` must equal the modelled residue count for every entity.** `parse_polymer`
  uses it as the token list, and entries absent from the model still **consume token indices** —
  so a too-long declaration silently shifts every residue index downstream. A common way to get
  this wrong: copied gemmi residues keep their source `subchain` label, `setup_entities()` then
  merges two chains into one entity and takes `full_sequence` from whichever came first. **Count
  the `_entity_poly_seq.mon_id` rows per entity and check them before spending.** The one
  legitimate exception is a `--keep-gaps` file, where the surplus declarations are exactly the
  unresolved positions and are deliberately tokenised as absent — see above.
- **`_atom_site.label_seq_id` must be written, and run 1..N per chain** (under `--keep-gaps`, it is instead the position in the declared range, so it skips the unresolved ones) — a `.` means that residue has no token. gemmi's `assign_label_seq_id()` only numbers residues whose `entity_type` is `Polymer`, and that flag is inherited from the source file: a plain PDB or a bare `_atom_site` mmCIF (most target files here) carries no entity annotation, so hand-built clones arrive as `Unknown` and every residue is skipped without a word. `boltz_template.py` sets it explicitly (fixed 2026-09-30) and raises if any residue is left unnumbered. The symptom, if it ever comes back: `--check` says *"declares N residues but N−k are modelled"* on a file that is otherwise correct, where k is the number of **adjacent identical residue pairs** in the sequence — the checker dedups residues on `(label_seq_id, comp_id)` and `(".", "PRO") == (".", "PRO")`. That is a real bug in the file, never a reason to loosen the check. Self-test: `uv run python utils/test_boltz_template.py`.
- Do not template the binder chain.
- Forced templates steer, they do not constrain. **Always superpose the predicted target chains
  back onto the template and exclude rows where the target did not land.** In one batch only
  9 of 27 landed, and the distorted ones produced the largest interfaces and the best hotspot
  recall in the set.

## Two operational traps

- **`Submitting N designs` counts SHARDS**, not designs (`--shard-size`, default 10). Two designs
  print `Submitting 1 designs`. Count `<out_dir>/boltz_inputs/*.yml` instead.
- **`boltz_shards/` is never cleared** (`mkdir(exist_ok=True)`). A rerun in a run_dir that already
  held a pilot will re-predict rows you meant to skip. Check the staged inputs first.
