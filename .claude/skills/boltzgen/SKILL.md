---
name: boltzgen
description: >-
  How to design binders against a protein target with the custom boltzgen tool on Modal — either
  NANOBODIES (--protocol nanobody-anything, the DEFAULT: new CDRs on a fixed, randomly drawn VHH
  framework) or DE-NOVO MINIBINDERS (--protocol protein-anything --binder-length 80-225). One
  BoltzGen run per target, designs as child rows. Load before composing a boltzgen run, and before
  reading any number it wrote. Three traps it is loaded for: boltzgen_protocol is a column and must
  be read FIRST, because a run that produced nanobodies for a caller who wanted minibinders looks
  identical in every other column; design_iptm is DILUTED BY CHAIN COUNT (a 4-chain average against
  a homotrimer) so design_to_target_iptm is the honest per-pair number; and an ALL-NA column has
  two causes — the metric was not produced, or the name the tool declares is misspelled — which
  only an --all-metrics re-collect tells apart. Covers why
  target geometry (bb_rmsd_target) is the FIRST gate and not the interface, the mapping-trust
  columns, the protocol-dependent step list that exists because `filtering` crashes and cannot be
  turned off, the TWO undocumented identifier frames an epitope depends on (residue numbers are
  label_seq_id not auth_seq_id; chain ids are label_asym_id for mmCIF but the author letter for
  PDB — on 1TUP author E,F,A,B are label A,B,C,D), the per-BATCH binder-length draw, and the fact
  that a binder-alone structure exists under protein-anything but NOT under nanobody-anything.
---

# boltzgen

**What it does.** Runs one of **two** BoltzGen protocols against a protein target (one chain or
many), inverse-folds with SolubleMPNN, refolds the complex with Boltz-2 and analyses it.
`action: create` — one BoltzGen run per target produces N designs, so the row count is known only
at collect time. `-t` is optional: with it, one run per ready row of that table; without it, a
root run that starts `table0`.

| `--protocol` | Binder | Length comes from |
| --- | --- | --- |
| `nanobody-anything` (**default**) | A VHH: framework fixed in sequence **and** structure, drawn **at random per design** from `--scaffolds`; CDR1/2/3 redesigned in sequence, structure *and length*. | the scaffold + CDR sampling |
| `protein-anything` | A **de-novo minibinder**: one chain, no framework, designed from nothing. | `--binder-length` (**required**) |

> **`boltzgen_protocol` is a collected column and is the first thing to read.** It comes from
> the group's own `runs/<group>/.submit.json`, written at submit beside the designs it describes
> — not from a flag at collect time and not from the run-level `.meta.json`, which is merged
> across submits. A run_dir can hold both protocols side by side. This column exists because a
> 60-design nanobody run was once delivered to a caller who had asked for minibinders, and
> nothing else in the table distinguished them.

> **Pre-column run_dirs: `boltzgen_protocol` is NA, and that is expected, not a failure.**
> The column was added 2026-10-09. Two run_dirs predate it and their groups' `.submit.json`
> files carry no `protocol` key:
>
> * `outputs/20261008_134905_boltzgen_verification`
> * `outputs/20261009_071201_boltzgen_pilot_1tnf_d143`
>
> **Both are `nanobody-anything`** — that was the only protocol the tool could run — but the
> collector will not say so, because it reports what the record says and the record is silent.
> It prints `.submit.json records no 'protocol'; it was written before boltzgen_protocol
> existed` and writes NA. **Do not read that NA as a broken attestation or as a mapping
> failure.** Re-collecting them will not fill it in; the fact was never recorded.
>
> **A resubmit into either is refused** (`.submit.json is missing a protocol`): a record that
> never recorded the protocol cannot prove a top-up matches. Use `-l/--dir-label` for a fresh
> leaf. The tables already collected from those runs are untouched.

### Re-collecting a pre-column run_dir: the complete list of what changes

Derived from the code, not estimated. For an **all-`OK` nanobody table** whose group records
carry no `protocol` key, a re-collect changes exactly this much — **everything not listed here
must be bit-identical, and that is the assertion worth testing.**

| Column | Before | After | Why |
| --- | --- | --- | --- |
| `boltzgen_protocol` | *absent* | **NA** | New column. The record never said. **Not** a claim of nanobody. |
| `boltzgen_scaffold` | the scaffold basename | **NA** | Fail-closed: the scaffold check runs only on an explicit `nanobody-anything` record. |
| `boltzgen_scaffold_src` | `in_run_scaffold_list` | **`not_applicable_no_scaffold`** | Same rule. The allow-list is still in the record, but nothing attests the protocol it belonged to. |
| `boltzgen_binder_refold_cif` | *absent* | **NA** | New column from the `ROW_KEYS` seed; nanobody never writes the file (#148). |
| `boltzgen_not_binding_resnames` | *absent*, unless the run used `--not-binding` | **NA** | New from the seed. It used to be written **only** when an anti-epitope existed. |
| `boltzgen_binding_resnames`, `_binding_frame`, `_binding_chain_frame`, `_chain_id_map` | *absent* **if the run used no steering at all**; otherwise present and unchanged | **NA** in that case | Same seed mechanism: these were written only inside `if epitope or anti_epitope`. A steered run keeps its values untouched. |
| any `CORE_METRICS` name absent from that run's CSV | *absent* (dropped) | **NA** | The `pd.NA` fill. **Which names these are is knowable only at collect time** — read the `not in this CSV, so written as NA: …` line the collect prints. |

Three things that do **not** change, and are the core of the non-regression check:

* **`boltzgen_status` stays `OK`.** The fail-closed scaffold rule withholds a name; it does not
  error a row. Nothing here turns a good design into a failure.
* **Every science column is identical** — all the confidence, RMSD, SASA, PLIP, liability and
  composition metrics that *were* in the CSV, plus `boltzgen_path`, `boltzgen_cif_path`,
  `boltzgen_binder_chain`, `boltzgen_sequence`, `boltzgen_cdr_sequence`,
  `boltzgen_design_resnums`, `boltzgen_fixed_resnums`, `boltzgen_n_design_res`,
  `boltzgen_binder_len`, `boltzgen_id`.
* **Every mapping-trust column is identical**: `boltzgen_chain_map`,
  `boltzgen_mapping_seq_identity`. The chain mapping does not depend on the protocol.

The three names removed from the declared set (`boltzgen_plip_hydrophobic_refolded`,
`boltzgen_target_id`, `boltzgen_free_cysteines`) do **not** reappear — they are no longer
declared, so the NA fill does not reach them.

**Input column.** `default_input_column` is **`pdb_path`** — the generic structure column, and
it is what `-i` falls back to. It is rarely the right one: the target usually arrives as
`chainsel_path` (a cropped target), `mkcomplex_path`, or whatever column actually holds the
structure in your table. **Pass `-i` explicitly.** With no `-t` the target comes from
`--target-structure` instead and `-i` is unused.

**Modal only.** There is no `SAPIA_ACTIVATE_BOLTZGEN`; the cluster has no BoltzGen install and
the ~6.4 GB of checkpoints are staged on a Modal Volume. Run it with `-e modal`.

> **Every column this skill documents exists on every row of every table this tool writes.**
> The structural, trust and attestation columns are seeded on each row before anything else, so
> a row that errored before the mapping ran still has `boltzgen_chain_map`,
> `boltzgen_scaffold_src` and the rest — as NA. **The trust filter below therefore cannot raise
> `KeyError`, including on a collect where every single design errored**, which is the shape it
> is most needed for. (Metric columns are a separate, protocol-dependent set; they are NA-filled
> against whichever set applies to the row's protocol. Tables collected **before 2026-10-09** are
> not rewritten and still lack columns — guard those with `if col in df.columns`.)

---

## Premise, and when its numbers are meaningless

### Shared by both protocols

| Assumption | What happens when it does not hold |
| --- | --- |
| **The binder is ONE designed chain.** | The collector refuses the design (`boltzgen_status` = `error: MappingError: expected exactly one designed chain …`; every mapping refusal is prefixed `error: <ExceptionType>: `). It never guesses which chain a sequence belongs to. |
| **The target is a real structure**, in whatever frame you handed it. | BoltzGen renumbers nothing you can rely on; residue indices you supply are `label_seq_id` (see below). |
| **The target stays where you put it.** | **It may not.** Upstream issue #226 (open, unanswered): on a multi-chain target a second protomer can MOVE during design, so the complex has the wrong binding geometry. `boltzgen_bb_rmsd_target` is the only thing that catches it. |
| **Designs are not reproducible.** BoltzGen has no `--seed` flag at all. | Deliberate. Under `nanobody-anything` the framework is `random.choice`d per design (`schema.py:1947-1949`) and `boltzgen_scaffold` records which one; under `protein-anything` the binder length is drawn randomly. Compare distributions, never design identity. |

### `nanobody-anything` only

| Assumption | Consequence |
| --- | --- |
| The framework is fixed; only the CDRs are free. | `boltzgen_cdr_sequence`, `boltzgen_design_resnums` and `boltzgen_fixed_resnums` describe a real CDR/framework split, and `fixed_resnums` is what `proteinmpnn --fixed-positions` consumes. |
| A scaffold allow-list exists. | `boltzgen_scaffold_src` can reach `in_run_scaffold_list`. |
| **There is no binder-alone structure** (issue #148). | Self-consistency needs `chainsel` to cut the binder out of `boltzgen_path`, then `usalign`. |
| `liability_modality=antibody` (`cli/boltzgen.py:99-105`). | `boltzgen_liability_*` are antibody developability rules. |

### `protein-anything` only

| Assumption | Consequence, and when the numbers mean nothing |
| --- | --- |
| **The whole chain is designed.** | `boltzgen_n_design_res == boltzgen_binder_len`, `design_resnums` is `1:N` and `fixed_resnums` is **empty** — correct, not a bug, and it means there is nothing to hand `proteinmpnn --fixed-positions`. `boltzgen_cdr_sequence` is **NA**: there is no CDR/framework split and `boltzgen_sequence` already carries the full sequence. |
| **There is no framework.** | `--scaffolds`/`--scaffold-dir` are **REFUSED at submit time**, `boltzgen_scaffold` is **NA**, and `boltzgen_scaffold_src` is `not_applicable_no_scaffold` — **never** `distinct_from_target`, which would read as a check that had passed. |
| **One binder length is drawn per BATCH, not per design.** Upstream's own `--diffusion_batch_size` help says so: *"all designs generated in the same batch will share the same length"*. Mechanically `multiplicity = ceil(num_designs / diffusion_batch_size)` (`cli/boltzgen.py:1056-1059`), each item re-parses the YAML (`data_from_yaml.py:247-266`), and the length is drawn once per parse (`schema.py:546`). | Upstream picks `diffusion_batch_size=10` for `--num-designs >= 100`, so a 200-design run with `--binder-length 80-225` samples only **20** distinct lengths out of 146. **Pass `--diffusion-batch-size 1` whenever the range matters.** Check it afterwards: `df["boltzgen_binder_len"].nunique()`. Same trap as rfd3's `--length min-max`. |
| **Liability rules are PEPTIDE rules, and they are meaningless here.** `protocol_configs["protein-anything"]` is `{}` (`cli/boltzgen.py:85`) — the base config — and `analysis.yaml` sets `liability_modality: peptide`. | Nothing carries the nanobody protocol's `antibody` rules across (the configs are separate dict entries with no inheritance), but the peptide motif set (`analyze_utils.py:921-941`) is `ProtTryp` on **every K/R**, an N-terminal cyclization motif and a DPP4 cleavage motif. On a folded 80-225 aa protein `boltzgen_liability_score` is close to a lysine+arginine count. **Do not gate on it, and never compare it with a nanobody row's.** Turn it off with `--config-set 'analysis liability_analysis=false'` if the noise is unwanted. |

This tool does **not** design Fabs, peptides or small-molecule binders, does not run affinity
prediction, and does not rank anything.

### The step list is PROTOCOL-DEPENDENT, and that is not an oversight

| `--protocol` | Forced `--steps` |
| --- | --- |
| `nanobody-anything` | `design inverse_folding folding analysis` |
| `protein-anything` | `design inverse_folding folding design_folding analysis` |

There is no `--steps` flag on this tool. Three facts force the omissions, **all re-verified for
both protocols** at the pinned SHA:

* generation never touches ESMFold2 — `task/predict/` contains no reference to it;
* `analysis` is launched with `esmfold2_metrics={not use_affinity}` and `analyze.py:638-652` then
  reads `esmfold2_scores/<id>.input.json` unconditionally → **it crashes** unless that key is
  false. `use_affinity` is `protocol in ["protein-small_molecule"]` (`cli/boltzgen.py:1259`), so
  it is False for both protocols here;
* `filtering` sets `self.use_esmfold2 = not use_affinity` at `filter.py:218` — **hard-wired,
  reachable by no config key** — and `filter.py:421` raises
  `ValueError("Missing ESMFold2 scores; run esmfold2_scoring and analysis before filtering")`.
  **It can only be omitted**, under either protocol.

> **Trap: anyone running BoltzGen with default steps gets a crash at the last step.** The
> signature is that `intermediate_designs_inverse_folded/` is full and
> `final_ranked_designs/` is absent, with that ValueError in the task `.err`.
> `--esmfold2_acceleration off` does **not** help — it selects a numerical backend, it does
> not disable the model.

> **Trap, worse than a crash: `protein-anything` WITHOUT `design_folding` collects NOTHING.**
> `cli/boltzgen.py:1237` sets `do_design_folding = protocol in ["protein-anything",
> "protein-small_molecule"]` and then passes `designfolding_metrics={do_design_folding}` to
> analysis at `:1310`. `filter_steps` (`:1363`) only drops steps from the list — it does not
> rebuild their arguments. So an analysis step with `designfolding_metrics=true` and no
> `fold_out_design_npz/` prints *"Folded path does not exist. Skipping"* and `return None` for
> **every** design (`analyze.py:972-977`), and `aggregate_metrics_analyze.csv` comes back with
> **zero rows** after the whole GPU spend. Not a crash — an empty table. The tool derives the
> step list from `--protocol` precisely so this cannot happen; it is the reason
> `--extra-flags '--protocol …'` is refused.

What `design_folding` **buys** under `protein-anything` is a real new capability: it writes
`refold_design_cif/<id>.cif` — **the binder folded alone, without the target**
(`task/predict/writer.py:94-104`) — collected as `boltzgen_binder_refold_cif`, plus the
`designfolding-*` metrics. So binder self-consistency is in the table directly, as
`boltzgen_designfolding-bb_rmsd_design`, with no `chainsel` + fold round-trip.

Consequences, both protocols: no `final_ranked_designs/`, no `results_overview.pdf`, no
`--budget`, and no `esmfold2_ipsae_min` (upstream's own ranking column). **Ranking happens in
the table with `-f`.** That is the right place for it anyway.

---

## Trust columns — read these before any science column

Three naming systems have to be reconciled before a single sequence can be believed, and
getting it wrong produces a plausible, wrong answer every time:

1. the CSV's sequence columns are keyed by **feature asym id**, not PDB chain label;
2. `source_context[...]["output_chain"]` is BoltzGen's **internal** chain name;
3. the written mmCIF label is `re.sub(r"\d+", "", internal_name)` and, **on a collision, an
   arbitrary letter from a reversed-alphabet pool** (`data/write/mmcif.py:63-80`) — so the
   label in the file need not be either of the other two.

| Column | Read it as | Expected |
| --- | --- | --- |
| `boltzgen_protocol` | **WHICH CLASS OF BINDER THIS ROW IS — read it before everything, including `_status`.** `nanobody-anything` or `protein-anything`, taken from the group's own `runs/<group>/.submit.json`. Written on **every** row including error rows, because a row that cannot say what it is is the problem the column solves. **NA has exactly two causes: the group's record is missing/unreadable (an error row anyway), or it predates this column — see "Pre-column run_dirs" below. NA is NEVER a failed attestation and NEVER a claim that the protocol was nanobody.** It decides the meaning of `scaffold`/`scaffold_src`, `cdr_sequence`, `fixed_resnums` and `liability_*`, and whether `binder_refold_cif` and `designfolding-*` exist at all. Never compare a `liability_*`, `cdr_sequence` or `scaffold` value across two different values of this column. | the protocol you asked for |
| `boltzgen_status` | **The only proof of success.** Anything but `OK` means no chain label, sequence, scaffold or residue list was written for that row — only the metrics, which are true of *some* chain pair but not of a named one. | `OK` |
| `boltzgen_mapping_seq_identity` | Identity between the CSV's `designed_chain_sequence` and the residues of chain `boltzgen_binder_chain` in the structure `boltzgen_path` points at. **This is about the chain MAPPING, not about sequence recovery against the parent framework** — hence the name. On an `OK` row it is always exactly `1.0`, because anything lower is refused; on the **identity-failure** error row the measured value is still written, so the number that condemned the design is visible. Other error kinds (a count mismatch, a corrupt file) leave it NA — they fail before an identity is computed. | `1.0` on `OK`; read it on `error` rows |
| `boltzgen_chain_map` | How the binder chain was tied to its `source_context` entry. **Only `ok` and `ok_renamed` mean the link holds.** `ok` — `output_chain`, digits stripped, equals the written label. `ok_renamed` — `output_chain` is present and differs, i.e. the mmCIF writer re-pooled the label; the link rests on chain **position**, which upstream enforces with `strict=True` (`source_context.py`), so it still holds. `no_source_context` — the whole `source_context` was JSON `null` or unreadable: **provenance is LOST, not benignly absent**; treat the row as untrusted and `boltzgen_scaffold` is NA. A per-chain missing `output_chain` cannot occur upstream (`update_designed` writes it for every entry, `source_context.py:200-205`) and is a hard **error** if it ever does. | `ok` or `ok_renamed` |
| `boltzgen_scaffold_src` | **How far `boltzgen_scaffold` itself was validated. Read `boltzgen_protocol` first.** **The scaffold check is FAIL-CLOSED: it runs only when the record explicitly says `nanobody-anything`.** `protein-anything`, an absent protocol (a pre-2026-10-09 record) and an unrecognised string all get `not_applicable_no_scaffold` + NA — the tool will not infer a protocol in order to publish a provenance verdict. `in_run_scaffold_list` — the name matched a scaffold **this group** was submitted with (read from `runs/<group>/.submit.json`, so a later submit with a different `--scaffolds` cannot re-point it). `distinct_from_target` — no allow-list existed (a `--spec` nanobody run), but the source differs from every target chain's. `not_applicable_no_scaffold` — **nothing was checked**; `boltzgen_scaffold` is NA. A separate verdict on purpose: a de-novo entity's `source_context` entry says `"inline"` (`schema.py:1607`) with no allow-list, so the nanobody branch would report the scaffold NAME `inline` under `distinct_from_target`, which reads as *checked and passed*. A withheld name costs a diversity analysis; a wrong one corrupts it. `not_in_run_scaffold_list` / `is_a_target_chain_source` / `no_source_context` — **NA**; the positional lookup produced something that is not a scaffold of this run, or is the target itself. | `in_run_scaffold_list` (nanobody) / `not_applicable_no_scaffold` (protein) |
| `boltzgen_n_design_res` | Number of designed residues — the CDRs under `nanobody-anything`, **the whole chain** under `protein-anything`. Cross-checked **three ways under `nanobody-anything`, FOUR under `protein-anything`** — the design-mask B-factor of the design mmCIF, the npz `design_mask`, the CSV's own `num_design`, and (protein only) **equality with `boltzgen_binder_len`**, since the whole chain is designed there. **Any disagreement is an error row, never a trimmed fit — and so is an absent `num_design`**, because the three-witness guarantee would otherwise silently become a two-witness one. The first three all count *designed* residues and so agree with each other; the fourth is the only one that catches a mask and a chain describing different objects. **A fourth-witness error row still carries `boltzgen_mapping_seq_identity`, `boltzgen_chain_map` and `boltzgen_scaffold_src`** — they were measured before the check fired, and they are what tell you whether the mapping was sound and only the count wrong. | never 0; `== binder_len` under `protein-anything` |
| `boltzgen_binding_resnames` | The residue **name** found at each `--binding`/`--not-binding` position, in the order you listed them (`A91=LYS,A128=GLU,…`), resolved at **submit** time against the target you passed. Eyeball it against your reference structure. NA when no `--binding` was given. If one **was** given and the record is missing, stale, or at odds with `specs/<group>.yaml`, every design of that group is an **error** — so NA can only ever mean "none asked for". | the residues you meant |
| `boltzgen_binding_frame` | Which RESIDUE frame those positions were read in, per chain (written for **either** steering claim, not just `--binding`): `declared_sequence` (the file declared a polymer sequence — mmCIF `_entity_poly_seq` or PDB `SEQRES`) or `synthesized_from_model` (it did not, so the sequence was synthesized from coordinates, exactly as BoltzGen's `parse_pdb` does). A `+`-joined pair means the epitope spans chains in different frames. NA when no epitope was given. | `declared_sequence` for a deposited structure |
| `boltzgen_binding_chain_frame` | Which CHAIN identifier BoltzGen resolved against: `label_asym_id` (mmCIF input) or `author_letter` (PDB input). See the chain-ID trap above. NA when no epitope was given. | matches your input format |
| `boltzgen_chain_id_map` | The author→BoltzGen chain translation, `H>A,I>B`. **Both sides, so a wrong translation is visible rather than inferred.** Identity map on a PDB. NA when no epitope was given. | identity, or a translation you expected |
| `boltzgen_not_binding_resnames` | The residues `--not-binding` steered the design AWAY from, same `A91=LYS` notation. **A separate column on purpose** — it is the opposite claim, and sharing one column would let a reader checking the epitope read anti-epitope residues as epitope residues. A `--not-binding`-only run is a legitimate invocation and collects normally, with the residue frame and chain map beside it. NA when no `--not-binding`. | the residues you meant to avoid |

### Where the epitope attestation lives, and what it is bound to

`boltzgen_binding_resnames` and its witnesses are a claim about **the designs in
`runs/<group>/`**, so the record lives there — `runs/<group>/.submit.json`, written at submit
time — and not in the run-level `.meta.json`. That matters because group names are
deterministic (`<stem>_bg`, or the parent row name), so re-running a target into the same
out_dir is the ordinary thing to do, and a run-level record keyed by name would have been
merged across submits and never invalidated by anything happening to the outputs.

Three things follow, all of which the collector enforces:

* a group resubmitted **without** `--binding` comes back with **no** epitope, not the previous
  one;
* a resubmit into a group dir that already holds designs is **refused at submit time** unless its
  record *proves* the same spec — BoltzGen is never given `--overwrite`, so those designs would
  survive. An absent, truncated or deleted record is refused for the same reason: no evidence is
  not the same as evidence of a match. Recovery is deleting the **whole** `runs/<group>/`, never
  the record alone;
* at collect, the record's `spec_sha256` and its list of requested claims are cross-checked
  against `specs/<group>.yaml`, which is rewritten on every submit. A disagreement errors every
  design of that group.

`.meta.json` now carries only run-level facts (`runs_root`, `protocol`, `steps`,
`root_designs`), and `root_designs` is rebuilt from the group dirs on disk — `prosapia` rewrites
that file *before* the tool's builder runs (`base_run.py:345` vs `:377`) and silently replaces a
corrupt one, so nothing a design's trust depends on may live there. There is **no fallback** to
an older layout: a group with no record of its own is an error, never a group that inherits an
attestation from somewhere else.

A **hand-written `--spec` YAML carrying `binding_types` cannot collect `OK`.** Its claims are read
back off the YAML, but this tool only resolves residue *names* for an epitope it built itself, so
the attestation is empty and every design of that group errors. That is the same bar the
flag-built path meets — if you want a steered run to collect, use `--binding`/`--not-binding`.

### The chain-ID trap — undocumented upstream, and it relocates whole epitopes

**Which identifier BoltzGen resolves a chain letter against depends on the input file
format**, and neither upstream README nor any doc says so:

| Input | What a `chain: id:` is matched against | Why |
| --- | --- | --- |
| **mmCIF** | **`label_asym_id`** | `parse_mmcif` builds chains from gemmi *subchains* (`mmcif.py:1051-1053,1073,1112`), and a subchain id **is** the `label_asym_id`. Matched at `schema.py:2299-2303`. |
| **PDB** | the **author letter** | `parse_pdb` round-trips through `make_mmcif_block()` (renaming `A` → `Axp`) and then `patch_chain_names` truncates to the **first character** with a collision pool (`pdb_parser.py:102,107-126`). |

> **Measured on 1TUP: author chains `E,F,A,B` are label chains `A,B,C,D`.** Asking for "A"
> meaning the p53 protein would have selected the DNA strand, and the epitope attestation would
> have come back as perfectly plausible nucleotides. The *common* case is milder and just as
> wrong: an author chain split into a polymer subchain plus ligand subchains shifts every later
> chain's `label_asym_id`.

**You always write the AUTHOR letters** — the ones a viewer shows — in `--target-chains`,
`--binding` and `--not-binding`. The tool translates them and records **both sides**:

* `boltzgen_binding_chain_frame` — `label_asym_id` (mmCIF) or `author_letter` (PDB);
* `boltzgen_chain_id_map` — e.g. `H>A,I>B`, author letter on the left. On a PDB it is the
  identity map. **If this shows a translation you did not expect, stop.**

Refused at submit time: an author letter that is not in the file; and, for a PDB, two chains
sharing a first character, where `patch_chain_names` re-letters one from a pool so the id
BoltzGen ends up with is not predictable.

> The upstream README adds to the confusion: `README.md:121` says to use "`label_asym_id`, and
> not the `auth_asym_id` author residue index". Those are **chain** ids; the *residue* index is
> `label_seq_id` / `auth_seq_id`, and the image it links (`assets/label_seq_id.png`) shows a
> residue index. The substance is right for both frames, the field name as written is wrong —
> do not follow it literally.

### The `label_seq_id` trap

**BoltzGen residue indices are 1-based positions in the target's DECLARED polymer sequence**
(mmCIF `_entity_poly_seq`, PDB `SEQRES`) — i.e. `label_seq_id`. The upstream README marks this
"IMPORTANT ⚠️". Two things follow, and the second is the one that bites:

* it is **not** `auth_seq_id` — a target a depositor numbered 18–155 is 1–138 to BoltzGen;
* it is **not** the count of *modelled* residues either. The declared sequence includes
  positions with no coordinates, and `parse_range` (`schema.py:649-667`) indexes straight into
  the array BoltzGen built from it (`mmcif.py:686-706`, one entry per `_entity_poly_seq` row,
  present or not).

> **Measured on 1TUP (p53), fetched 2026-10-08:** chain A declares 219 residues and models 196;
> chains B and C begin at `label_seq_id` 3 and 2. Counting modelled residues — which an early
> version of this tool did — reports a *different, perfectly plausible* residue for every
> position after the first gap, and `boltzgen_binding_resnames` would then confirm the wrong
> frame. This tool now replicates BoltzGen's own parser instead (including `SEQRES` alignment
> and the synthesis step for a PDB with no `SEQRES`), and reports which frame it used.

The defences are all mandatory and all present:

* a `--binding` position **absent** from the declared sequence raises at submit time, naming
  the chain and its real range;
* a position that is declared but has **no coordinates** also raises. BoltzGen would accept it
  (the residue is parsed with `is_present=False`) and it would steer nothing, while the table
  asserted an epitope;
* a `--binding` chain that `--target-chains` does not include raises — otherwise the run either
  dies inside BoltzGen or proceeds with **no epitope steering at all**;
* `boltzgen_binding_resnames` and `boltzgen_binding_frame` put the answer in the table;
* `boltzgen check` runs on the written spec **before** any GPU work (see below).

> Also note upstream issue #215: `binding:` is a **soft preference, not a constraint** — a
> listed hotspot can be ignored. And issue #217: `binding` and `not_binding` may not work on
> the **same chain**; put them on different chains or verify on a small run first.

---

## The questions it answers, in the order the gates must be applied

> **This section gives no numeric cutoffs, deliberately.** The tool records scores; which
> value separates a design worth keeping from one worth discarding is a scientific decision
> that belongs to whoever is running the campaign, with the distribution of their own run in
> front of them. **No observed range is given for any column.** A 60-design `nanobody-anything`
> run exists (2026-10-09) but its distributions have not been written down here, and
> `protein-anything` has never been run at all. A threshold written into a skill is copy-pasted
> long after it has gone stale —
> `polyv`'s first skill shipped `sc >= 0.55` when that tool's measured range topped out at
> 0.40, which would have selected zero rows forever. What *is* given below is the **order** the
> gates must be applied in, and which direction is good; take the numbers from your own run and
> record them in the campaign log.

**Gate 1 — did the target hold still? `boltzgen_bb_rmsd_target` (Å, lower is better).** This is
the first gate and it is not optional. On a multi-chain target, issue #226 says a protomer can
move; a complex whose target moved has the wrong geometry, and **a large interface measured on
it is a fiction**. Read this before `delta_sasa_refolded`, before `design_to_target_iptm`,
before anything. It is a **floor for discarding**, not a ranker: a design that fails it is
unusable, but a lower value among the survivors does not make a design better. What counts as
"moved" depends on how many chains the target has and how much of it was cropped — plot the
distribution of your own run, pick the bar, and write it in the campaign log.

> **MEASURED 2026-10-09 — under `protein-anything`, `design_iptm` is NOT an independent number.**
> In the raw `aggregate_metrics_analyze.csv` of a minibinder run, `design_iptm`,
> `design_to_target_iptm` and `design_residue_iptm` are **identical in every row** (10/10
> checked; each has 10 distinct values, so they are not constant). Under `nanobody-anything` the
> same CSV has `design_iptm` differing from the other two in **10/10** rows. **This is upstream,
> not our collector** — the collector copies the CSV names verbatim and was verified to alter
> nothing.
>
> **Which column is still honest?** `design_to_target_iptm` ≡ `design_residue_iptm` in **both**
> protocols — that pairing is stable, and it is `design_iptm` (the average) that moved onto them
> under `protein-anything`, not the reverse. The parsimonious reading is therefore that
> **`design_to_target_iptm` remains the per-pair number and is still the one to gate on**; what
> you lose under `protein-anything` is `design_iptm` as a *separate* diluted figure, which now
> carries no additional information. **That is an inference from column identity, not proof** —
> settling it means reading how upstream computes `design_iptm` per protocol. Until then, do not
> quote `design_iptm` as a chain-diluted value on a minibinder run, and do not treat the two as
> two pieces of evidence.
>
> Also measured on the same run: **every `designfolding-*iptm` column is constant 0**, because
> the design-folding step folds the binder **alone with no target** — there is no interface to
> score. Expected, and a useful confirmation that the design-folding block really is
> target-free. For the same reason `designfolding-bb_rmsd_design` ≡ `designfolding-bb_rmsd`
> (upstream, 10/10) — with no target, "the design" and "everything" are the same atoms.

**Gate 2 — is the pose confident? `boltzgen_design_to_target_iptm` (higher is better),
`boltzgen_min_interaction_pae` (Å, lower is better).**

> Note `no_source_context` is deliberately **not** in that list: a null `source_context` is a
> loss of provenance, not a benign absence.

> **Trap, and the single easiest misreading of this tool's output: `design_iptm` is diluted by
> chain count.** It is an iPTM averaged over every chain pair in the complex. Against a
> homotrimer the complex has **four** chains, so a nanobody that binds one protomer perfectly
> is averaged against pairs it was never meant to form, and a perfectly good design reads low.
> **`design_to_target_iptm` is the per-pair number and the honest one.** Gate on that.
> (Same arithmetic as boltz's `iptm` vs `pair_chains_iptm`; the boltz skill says the same
> thing about the same quantity.) The dilution also means **no iPTM bar you have seen quoted
> for a two-chain complex transfers here** — the number is not on the same scale.

**Gate 3 — does the design fold back to what was designed? `boltzgen_bb_rmsd_design`,
`boltzgen_bb_target_aligned_rmsd_design`.** The second is the stricter one: it aligns on the
**target** and then measures the design, so it tests the binding pose, not just the fold.

> Upstream issue #225: on a fixed VHH framework, **every** design failed the default
> `refolding_rmsd_threshold`, even at 10 Å, and the reporter suspects `filter_rmsd_design` is
> the right metric instead. We never run `filtering`, so that default cannot silently return
> zero designs here — but **a conventional RMSD bar read off a `protein-anything` paper or
> README does not transfer to nanobodies**, and the reverse is equally untested.

**Gate 3b — `protein-anything` only: does the binder fold on its own?
`boltzgen_designfolding-bb_rmsd_design` (Å, lower is better).** The binder sequence folded
**without the target**, compared back to the designed backbone. A minibinder that only adopts
its shape when the target is holding it there will pass gate 3 and fail this one. There is no
nanobody equivalent — a VHH framework is a folded protein by construction, which is part of why
upstream disables the step for it. **Still Boltz-2 checking Boltz-2**; an independent predictor
is a separate step.

**Gate 4 — is the interface real? `boltzgen_delta_sasa_refolded`, `boltzgen_plip_*_refolded`.**
Only after gates 1–3.

> **What `delta_sasa_refolded` actually is, because this workspace has already been bitten by
> a factor of two here (`polyv`).** `analyze_utils.py:662-716`, `get_delta_sasa`, docstring at
> `:663` — "Return target SASA lost on binding, unbound SASA and bound SASA". Both partners
> enter the bound surface (`bound_mask = atom_design_mask | atom_target_mask`, `:666`) but
> **only target atoms are summed** (`:699`), the unbound reference is the target alone
> (`:668`, `:709`), and `delta = target_area.sum() - target_bound` (`:715`). Probe 1.4 Å,
> `point_number=960`, custom per-atom radii (`:686-696`, `:709-714`). Masks are not swapped at
> either call site (`analyze.py:774-778`, `:1168-1172`, both passed by keyword).
>
> So it is **one-sided BSA, NOT summed over both partners**. `polyv_dsasa` IS summed over both,
> so it is roughly **2×** this — before probe, radii and structure differences. Do not compare
> them, and do not compare either with PyRosetta's `if_dSASA` without checking that tool's
> convention first.
>
> **The `design_sasa_*` columns are misnamed upstream, in all four variants.**
> `get_delta_sasa` returns `(delta, target_area.sum(), target_bound)` (`:716`) — all three
> target-side — and `analyze.py:1174-1176` (refolded) and `:779-781` (original) unpack the last
> two into `design_sasa_unbound_*` / `design_sasa_bound_*`. They are the **target's** bound and
> unbound SASA. No upstream issue reports it. Collected under upstream's spelling so they match
> the raw CSV; the meaning is here.
>
> One more: the "design" side of the mask is `resolved_atoms_chain_mask`, built from
> `chain_design_mask` — the **whole designed chain including its fixed framework**
> (`analyze.py:695-696`), not just the CDRs.

**Gate 5 — is it developable? Protocol-dependent, and the two are NOT comparable.**

* `nanobody-anything`: `boltzgen_liability_score`, `boltzgen_liability_num_violations`,
  `boltzgen_design_hydrophobicity`. The pinned SHA sets `liability_modality=antibody` for this
  protocol — **PyPI 0.3.2 does not** (issue #228, a PR was made 2026-05-19 and never merged),
  which is one of the reasons this tool pins a git SHA.
* `protein-anything`: `boltzgen_design_largest_hydrophobic_patch_refolded` and
  `boltzgen_design_hydrophobicity`. **`liability_*` is a PEPTIDE motif panel here and should
  not gate anything** — see the premise table. There is no protein-appropriate liability
  modality in BoltzGen at all (`get_motif_set` asserts `modality in {"antibody","peptide"}`,
  `analyze_utils.py:942`).

> **`boltzgen_free_cysteines` is NOT collected, under either protocol.** This skill used to
> list it. Upstream guards it with `self.free_cys` (`analyze.py:920`), whose constructor default
> is `False` (`:124`) and which no protocol config and no key in `analysis.yaml` ever sets —
> so it never reaches the CSV. Confirmed on the real 60-design run: absent from the table under
> every spelling. Count unpaired cysteines from `boltzgen_sequence` yourself if you need it.

---

## Invocation and flags

```bash
# NANOBODIES, root run: one BoltzGen run against a target not in any table yet
sapia run boltzgen <run_dir> -e modal \
    --target-structure inputs/trimer.cif --target-chains A,B,C \
    --binding A91,A128,B91,B128 --num-designs 50

# DE-NOVO MINIBINDERS, 80-225 aa, against the same target.
# --diffusion-batch-size 1 is NOT decoration: without it, upstream picks 10 for
# --num-designs >= 100 and a 200-design run samples only 20 of the 146 lengths.
sapia run boltzgen <run_dir> -e modal \
    --protocol protein-anything --binder-length 80-225 \
    --diffusion-batch-size 1 \
    --target-structure inputs/trimer.cif --target-chains A,B,C \
    --binding A91,A128,B91,B128 --num-designs 200

# child run: one BoltzGen run per target already in a table.
# -i names the column holding the TARGET structure. default_input_column is
# 'pdb_path' and is rarely right -- pass -i explicitly.
sapia run boltzgen <run_dir> -t table0 -i chainsel_path -e modal \
    --binding 'A{epitope_start}-{epitope_end}' --num-designs 200

sapia collect boltzgen <run_dir> -t <the table the run reserved>
```

| Flag | Default | Notes |
| --- | --- | --- |
| `--protocol` | **`nanobody-anything`** | `nanobody-anything` or `protein-anything`. **The default is deliberate and unchanged**: every invocation written before 2026-10-09 behaves identically. `protein-anything` is strictly opt-in. It decides the step list, the spec shape, and the meaning of five columns. Recorded per group and collected as `boltzgen_protocol`. |
| `--binder-length` | none | **`protein-anything` ONLY, and REQUIRED there** (there is no default length). `N`, `MIN-MAX`, or BoltzGen's own `MIN..MAX`; inclusive at both ends. `{expr}` islands resolve per design. **REFUSED under `nanobody-anything`** — the framework fixes the chain and the CDR lengths are sampled from the scaffold. **The `MIN-MAX` form is TRANSLATED, not forwarded:** BoltzGen tokenises a sequence spec with `re.findall(r"\d+\.\.\d+\|\d+\|[a-zA-Z]", …)` (`schema.py:947`), so the literal `80-225` would be two design blocks of 80 and 225 — a **305-residue** binder. Anything richer than a number or a range (`3..5C6C3`, comma-separated blocks) is refused; use `--spec` for those. **One length per BATCH — see the premise table.** |
| `--binder-chain` | auto (first free letter of `ZYX…`) | **`protein-anything` ONLY, and ADVISORY — the id you ask for does NOT appear in the output.** See the measurement below. The flag still writes the `protein:` entity's `id:` and is preserved as `source_context`'s `source_chain`, and an explicit letter the target already uses is still **REFUSED** rather than re-pointed (`Structure.concatenate`, `schema.py:1633-1638`, would rename one of the two and only log it) — kept as cheap insurance in case upstream changes. But **never assume the requested letter is the output letter: read `boltzgen_binder_chain`.** **REFUSED under `nanobody-anything`.** |
| `--target-structure` | none | Root run only. With `-t` the target comes from `--input-column`. Group name is `<stem>_bg`. |
| `--target-chains` | `all` | prosapia chain mini-language (`,` separates, `:` is an inclusive letter range): `A`, `A,B,C`, `A:C`. **AUTHOR chain letters**, matched **case-sensitively** — the ones a viewer shows; the tool translates them (see the chain-ID trap). **An explicit list writes `include:` with the POLYMER subchains only, so any cofactor, metal or glycan on those chains is DROPPED — the default `all` keeps them.** The submit prints what it drops; the two modes design against materially different targets. |
| `--binding` | none | Epitope, chain-prefixed: **AUTHOR chain letter** + **1-based position in the DECLARED polymer sequence** (`label_seq_id`), e.g. `A91,A128,B91-95`. `{expr}` islands resolved per design up the lineage. Soft preference, not a constraint (#215). |
| `--not-binding` | none | Same syntax and frame. See #217 before putting both on one chain. |
| `--scaffolds` | `7eow,7xl0,gontivimab,isecarosmab,sonelokimab` | **`nanobody-anything` ONLY — REFUSED under `protein-anything`, not ignored.** The VHH frameworks BoltzGen draws from, **at random, one per design**. Bare names resolve against `--scaffold-dir`; a name with `/` is a path. `8coh`, `8z8v` and `gefurulimab` also ship. One name pins the framework and kills the diversity. |
| `--scaffold-dir` | `$BOLTZGEN_SCAFFOLD_DIR` or `/opt/boltzgen/example/nanobody_scaffolds` | **`nanobody-anything` ONLY — REFUSED under `protein-anything`.** Inside the task container (Modal only). **Not checked at submit time** — it does not exist on the workstation. |
| `--num-designs` | **50** | **Upstream's default is 10000 and its advice for real work is 10,000–60,000.** 50 is upstream's own smoke-test size and is the default here *because per-design runtime on this image is unmeasured*. Raise it once a verification run has timed it. Issue #139: values above 10000 are rejected. |
| `--diffusion-batch-size` | upstream auto (**1 below 100 designs, else 10**) | **Every design in one batch shares the same sampled length** — the CDR lengths under `nanobody-anything`, the whole binder length under `protein-anything`. A large batch relative to `--num-designs` under-samples lengths. **The default flips to 10 at `--num-designs 100`, which is exactly where a real campaign starts**, so a 200-design minibinder run gets 20 lengths unless you pass `1`. Same trap as rfd3's `--length min-max`. |
| `--num-workers` | `1` | → `--num_workers`. |
| `--spec` | none | A BoltzGen design YAML used verbatim (still `{expr}`-resolved). **`boltzgen_protocol` then attests THE FLAG, NOT THE YAML**: your spec decides what is actually designed, while `--protocol` only selects BoltzGen's step list and config defaults. The two can disagree and nothing can detect it — the column records what was requested, which is all this tool knows. (The resubmit guard compares the recorded protocol explicitly for this reason: the YAML hashes identically under either protocol.) Mutually exclusive with the spec-building flags — passing both is refused, not silently ignored. **Use absolute paths inside it**: BoltzGen resolves relative paths against the *original* YAML's directory, not the resolved copy's. **No chain translation happens on a `--spec` YAML** — you write the ids BoltzGen resolves (`label_asym_id` for an mmCIF target), and `boltzgen_chain_id_map` is NA, so nothing in the table attests them. |
| `--extra-spec` | none | YAML/JSON mapping merged into every written spec. A key this tool also writes (`entities`) is an error. **`{expr}` is resolved but chain ids are NOT translated** — any chain id inside it reaches the YAML verbatim and is not covered by `boltzgen_chain_id_map`, so write the ids BoltzGen resolves (`label_asym_id` for an mmCIF target). |
| `--config-set` | none | Repeatable, `'<step> <key>=<value>'` → `--config <step> <key>=<value>`. Appended after the forced `analysis esmfold2_metrics=false`; **setting that key to anything but false is refused at submit time**. |
| `--extra-flags` | `""` | Verbatim `boltzgen run` tokens. **Refused at submit time: `--protocol`, `--steps`, `--output`, `--num_designs`, `--config`, `--cache`, `--overwrite`, `--reuse`, `--diffusion_batch_size`** — argparse is last-wins, so each would silently override the tool's own value. The last three were added 2026-10-09: `--overwrite`/`--reuse` falsify the resubmit guard's premise, and `--diffusion_batch_size` has a dedicated flag whose value the submit log reports. |
| `--metrics` (collect) | `""` | Extra CSV columns on top of the core set. |
| `--all-metrics` (collect) | off | Every column (~100). Names containing `<`/`>` (upstream's threshold booleans) are still skipped — a `--filter` cannot address them. |

**Resources (Modal only).** `RESOURCES = {"gpu": "H100", "cpu": 8, "memory": "64G", "timeout":
"12:00:00"}`. GPU is **not** forced to 0 — this is a GPU tool. Override the GPU per run with
the base flag `--modal-gpu` (`--gpu-type` does not exist). 64 GB of host RAM is not a guess:
the official SLURM example asks for it and issue #208 reports the `analysis` step OOM-ing at
16 GB with only 100 designs.

**Weights (Modal only).** ~6.4 GB on the Volume `SAPIA_MODAL_VOLUME_BOLTZGEN_CACHE`
(default `boltzgen-cache`), mounted **read-write** at `/boltzgen_cache` and passed as
`--cache`. Read-only reuse is undocumented upstream, and its own Dockerfile chowns the cache
to the runtime user. **Warm it with one small run before fanning out**, or several first
containers each pull the same 6.4 GB.

**First build is long.** The image clones BoltzGen at a pinned SHA and installs ~100 exactly
pinned wheels including torch 2.10.0 and the CUDA-12 stack. Allow 15+ minutes before assuming
a first submit is stuck.

---

## Gotchas, each with its signature

| Gotcha | What you will see |
| --- | --- |
| **`filtering` crashes and cannot be configured off.** | Not reachable through this tool — but if anyone runs `boltzgen` by hand with default steps: `ValueError: Missing ESMFold2 scores; run esmfold2_scoring and analysis before filtering` in the `.err`, with `intermediate_designs_inverse_folded/` full and `final_ranked_designs/` absent. |
| **`design_iptm` read as the binder's iPTM.** | Good designs cluster at a suspiciously low `design_iptm` on a multi-chain target while `design_to_target_iptm` is high. The gap *is* the chain-count dilution. |
| **Epitope numbers copied off a crystal structure.** | At submit time: `--binding: label_seq_id [91, 128] absent from chain A of trimer.cif, whose polymer is 1..138`. If the numbers happen to exist but mean something else, `boltzgen_binding_resnames` shows the wrong amino acids. |
| **`rmsd` / `rmsd_design` / `plip_hbonds` / `delta_sasa_original` are missing.** | They are **not collected**, deliberately. Upstream passes `allatom_fold_metrics` / `delta_sasa_original` / `noncovalents_original` `= skip_inverse_folding`, which is **False** here, so it never computes them. Use `bb_rmsd*`, `delta_sasa_refolded`, `plip_*_refolded`. |
| **A declared metric is absent from the CSV.** *(Behaviour CHANGED 2026-10-09 — see below.)* | One printed line per collect, naming every absent column: `not in this CSV, so written as NA (the column EXISTS and a -f filter naming it will not raise; …): …`. **Since 2026-10-09 the column is written as `pd.NA` rather than dropped**, so the table's shape is a function of the tool version rather than of what one CSV happened to contain, and a `-f` filter naming any documented column is safe. **NA is the absence of a measurement, never 0.** |
| **A `boltzgen_` column is entirely NA across every row.** | **TWO possible causes, and NA alone cannot tell them apart.** Either (a) the metric genuinely was not produced for this step set / protocol, or (b) **the name this tool declares is misspelled** and the real CSV column is called something else. Since the `pd.NA` change there is no `KeyError` to distinguish them — a typo now presents as a plausible all-NA column. **The way to tell: re-collect with `--all-metrics`, which takes the CSV's own column names, and compare.** This matters most for the five `protein-anything` names (`boltzgen_designfolding-*`, `boltzgen_design_largest_hydrophobic_patch_refolded`), whose spelling is read off `analyze.py` and **has never been observed in a real CSV**. Also read the collect's printed "written as NA" line — it names exactly which declared columns were missing, and it is now the only live detector. |
| **A `-f` filter dies with `KeyError: 'boltzgen_<something>'`.** | On a table collected **before 2026-10-09** this is the old dropped-column behaviour: `table0_a100s1.tsv` and `table0_h100.tsv` of `20261009_071201_boltzgen_pilot_1tnf_d143` have 57 columns each and carry none of `boltzgen_free_cysteines`, `boltzgen_plip_hydrophobic_refolded`, `boltzgen_target_id` under any spelling. **Old tables are not rewritten**, so guard with `if col in df.columns` when filtering one. On a table collected after, a `KeyError` means a name this skill does not document — check `[c for c in df.columns if c.startswith("boltzgen_")]`. |
| **`boltzgen_design_ipsae_min` / `_design_to_target_ipsae` / `_target_to_design_ipsae` are all `0.0`.** | **Measured 2026-10-09: identically `0.0` in all 60 rows of a real run — present, but carrying no information.** They are in the CSV and are collected, so they do not trip the dropped-column trap; they simply have zero variance. **Do not rank on them.** Upstream's own ipSAE ranking column is `esmfold2_ipsae_min`, which needs the ESMFold2 step this tool never runs. Why these three come back zero is **unknown**. |
| **Sibling boltzgen tables cannot be merged on `name`.** | BoltzGen derives a design's `id` from the target stem plus an index, so **every shard of the same target produces the same 10 names**. Measured 2026-10-09: six sibling tables of one pilot gave 10 distinct names, each repeated 6×. A `pd.merge(..., on="name")` across them silently fans out to a cross product. Merge on `(parent_table, name)` or `(table, name)`, or re-key before merging — and prefer one table per shard with `--table-label`. |
| **`protein-anything` with no `design_folding`.** | `aggregate_metrics_analyze.csv` comes back with **zero rows** and the `.out` is full of `Folded path does not exist. Skipping`. No crash, no error — just a header-only CSV after the whole GPU spend. Not reachable through this tool (the step list comes from `--protocol`), but the collect names it anyway: a zero-row group prints **`NO DESIGNS -- aggregate_metrics_analyze.csv exists but has zero rows. This is not an empty success…`** and specifically **never prints `OK`**, because "OK (0 design(s))" is how this failure used to read. |
| **`--scaffolds` or `--binder-length` under the wrong protocol.** | At submit time: `--protocol protein-anything does not use ['--scaffolds']. They would be silently ignored, so they are refused instead…`. The refusal names **every** mismatched flag, not just the first. `protein-anything` with no `--binder-length` is refused the same way: `needs --binder-length`. |
| **The same `--spec` YAML resubmitted under a different `--protocol`.** | At submit time: `already holds designs made under a DIFFERENT PROTOCOL (recorded … this run …). The spec hash matches, which under --spec it will for ANY protocol`. **A `--spec` YAML is byte-identical whichever protocol you pass**, so the spec hash alone cannot tell the two submits apart — the guard compares the recorded protocol and step list explicitly. Without it the second submit would rewrite `.submit.json` and stamp its own protocol onto the first's designs, which is exactly what `boltzgen_protocol` exists to prevent. Use `-l/--dir-label`. |
| **A resubmit into a pre-2026-10-09 group.** | `already holds designs but its .submit.json is missing a protocol`. Intended: a record that never recorded the protocol cannot prove a top-up matches, and the guard requires a readable **match**, never merely the absence of a known mismatch. Use `-l/--dir-label` for a fresh leaf. |
| **`--extra-flags '--overwrite'` / `'--reuse'` / `'--diffusion_batch_size 10'`.** | Refused at submit time. `--overwrite` and `--reuse` **falsify the resubmit guard's own stated premise** ("BoltzGen is never run with `--overwrite`, so those designs would survive") — with them the designs the guard protects are destroyed instead. `--diffusion_batch_size` has a dedicated flag and argparse is last-wins, so passing it here silently beats `--diffusion-batch-size 1` **while the submit log still prints the per-batch warning computed from the tool flag** — you are told one thing and get another, on the parameter that decides whether a length range is sampled per design or per batch of ten. |
| **A `protein-anything` design whose mask does not cover the whole chain.** | An error row: `under protein-anything the whole binder chain is designed, so the designed-residue count must equal the chain length -- but N residues are masked as designed in a chain of M`. The three count witnesses all count *designed* residues and agree with each other; this is the only check that catches a mask and a chain describing different objects, which would displace every residue-indexed column. **The row keeps `boltzgen_mapping_seq_identity`, `boltzgen_chain_map` and `boltzgen_scaffold_src`** — `map_binder` had already returned when the check fired, so those are real measurements. `mapping_seq_identity == 1.0` beside this error means the chain mapping was sound and *only* the mask count is wrong, which is a very different diagnosis from a mapping that is adrift. Everything residue-indexed or count-derived is NA. |
| **`--binder-length 80-225` producing 305-residue binders.** | It cannot, through this tool: the hyphen form is translated to `80..225` and anything else is refused. If you ever write `sequence: 80-225` into a `--spec` YAML by hand, BoltzGen's tokeniser (`schema.py:947`) reads it as two design blocks and you get 305 residues with nothing to say so. Check `boltzgen_binder_len` against what you asked for. |
| **Every `boltzgen_binder_len` is one of a handful of values.** | The per-BATCH length draw — **MEASURED, not inferred, see below.** Count them: `df["boltzgen_binder_len"].nunique()` should be about `num_designs / diffusion_batch_size`. Re-run with `--diffusion-batch-size 1`. |

### MEASURED 2026-10-09: the per-batch length draw is real, and it buys throughput

Two 10-design `protein-anything` runs against the same 3-chain target, `--binder-length 80-120`,
H100, identical in every respect except `--diffusion-batch-size`
(`outputs/20261009_111006_boltzgen_protein_verification`, `table0` and `table0_batch10`):

| | `--diffusion-batch-size 1` | `--diffusion-batch-size 10` |
| --- | --- | --- |
| `binder_len.nunique()` over 10 designs | **9** | **1** (all ten were 112) |
| total wall clock | 11 min 42 s | 9 min 14 s |
| startup (task start → first design file) | ~92 s | ~112 s |
| **per design, post-startup** | **~61 s** | **~44 s** |
| design-stage file span | 2 min 50 s | **7 s** |

**Both halves of this matter and they pull in opposite directions.** Batching really does
collapse length diversity to one draw per batch — upstream's `--diffusion_batch_size` help text
is correct. It is also **~28 % faster per design**, and the gain is concentrated entirely in the
diffusion stage (a 2 min 50 s file span becomes 7 s); the refold stages are per-design and do not
batch, which is why the end-to-end gain is smaller than the design-stage gain.

**The trade-off, priced.** At 200 designs over an 80–120 range:
`--diffusion-batch-size 1` ≈ 3.4 GPU-h and ~200 length draws;
`--diffusion-batch-size 10` (upstream's auto-choice above 100 designs) ≈ 2.4 GPU-h and **20**
length draws. **Full length diversity costs roughly 40 % more GPU time.** Decide which you are
buying and write it in the campaign log — do not let upstream's auto-switch decide it for you.

Order-of-magnitude for a large campaign, from the same rates: 10,000 designs ≈ **122 GPU-h**
batched, ≈ **169 GPU-h** unbatched. Measured at 10 designs, so treat as a planning figure.
| **A design that cannot be mapped.** | `boltzgen_status` starts with `error:`, `boltzgen_path` is **empty**, `boltzgen_cif_path` and (protein) `boltzgen_binder_refold_cif` still point at the written files so you can look, and every chain/sequence/scaffold/residue column is NA. The metrics are still there. |
| **`-i boltzgen_cif_path` or `-i boltzgen_binder_refold_cif` picks up REFUSED designs.** | Both are deliberately left populated on a mapping-error row when the file exists — the policy is "a failure must be inspectable". But `filter_ready` keys on a column being **non-empty**, not on `<leaf>_status`, so a downstream step reading either column gets the refused rows too. **`boltzgen_path` is the safe input column**: it is empty on every error row by construction. Otherwise filter `boltzgen_status == "OK"` first. |
| **Two runs of identical commands give different designs.** | Expected. No seed exists; the scaffold is a random draw. Compare `boltzgen_scaffold` distributions, not design identity. |
| **`mols.zip` download hangs at 0%.** | Upstream issue #216: a TLS-intercepting proxy blocking `cas-bridge.xethub.hf.co` — the same failure mode this network has with Modal blob storage. Stage the Volume from a container; do not weaken TLS. |
| **A spec YAML passed with `--spec` references files relatively.** | BoltzGen resolves them against the *original* YAML's directory. The resolved copy lives in the run_dir, so a relative path points at nothing. Use absolute paths. |
| **`boltzgen check` fails and the task stops before any GPU work.** | Intended. The task script validates the written spec first and `set -e` kills the task on a non-zero exit; the error is in the `.err` and a colour-coded mmCIF of the binding site is under `<out_dir>/runs/<group>/spec_check/`. This is the only independent test that the YAML nesting is what BoltzGen's schema wants — without it a schema parser that ignores an unknown key would turn a shape error into "no epitope was requested" while `boltzgen_binding_resnames` still reported one. |
| **`--extra-flags '--protocol protein-anything'` (or `--steps`, `--output`, `--num_designs`, `--config`, `--cache`).** | Refused at submit time with `--extra-flags may not contain [...]`. These are interpolated after the tool's own flags and argparse is last-wins, so they would silently change what the tool IS while every column, docstring and runtime guard still said nanobody — and the single-designed-chain guard would still pass. `--config-set` likewise refuses to set `esmfold2_metrics` to anything but false. |
| **An epitope residue with no coordinates.** | At submit time: `label_seq_id [128] on chain A ... are declared but have NO COORDINATES`. BoltzGen would accept it and steer nothing. Pick a modelled residue or crop the target. |
| **An epitope on a chain `--target-chains` excludes.** | At submit time: `--binding: chain(s) ['C'] are not in --target-chains`. Without the guard the run either dies inside BoltzGen or proceeds with no steering while the table asserts an epitope. |
| **A duplicated `id` in `aggregate_metrics_analyze.csv`.** | One `error:` row naming it, never a silent last-wins collapse of two designs into one. The by-name join (`file_name == id + ".cif"`) is ambiguous for every artifact of that design. |
| **One corrupt `.npz` in a 12-hour run (upstream issue #223, intermittent `BadZipFile`).** | That design becomes an `error:` row naming the exception type; the rest of the collect completes. It does **not** abort the collect — which it used to, losing the whole run's table and collecting nothing on re-run until the file was deleted by hand. |
| **Every `boltzgen_scaffold` is NA and `scaffold_src` is `no_source_context`.** | `source_context` was null or unreadable for every design; the collector prints the exception type per file. The sequences and chain labels are unaffected — only provenance is lost. |
| **Every `boltzgen_scaffold` is NA with a mix of other verdicts.** | The collector prints `NO design got a validated scaffold` plus the verdict counts. `source_context["chains"]` holds an entry per **non-water chain**, polymer or not (`source_context.py:14-30` has no type filter), zipped `strict=True` against the written chains — so the correspondence is enforced, but only if both sides count the same chains. A target carrying a ligand, ion or glycan chain is where to look. |
| **A target with microheterogeneity.** | Refused at submit time: `chain A has microheterogeneity in its declared sequence`. Upstream appends one array slot per `_entity_poly_seq` row with no dedup (`mmcif.py:998`) and then trips its own `assert len(set(res_indices)) == len(res_indices)` (`:803-806`) — so today it is a hard parse failure *inside a GPU task*. Refusing at submit time costs a submit instead. **If that assert is ever relaxed the silent off-by-N returns**, which is why this guard exists rather than relying on upstream. (`.pdb` inputs are pre-cleaned by `remove_alternative_conformations()`, `pdb_parser.py:36-38`, and no shipped example exercises the path.) |
| **A run refused most of its targets.** | `Submitting N of M target(s); K REFUSED` plus one line per refusal. **A skipped target leaves no status anywhere** — a create tool cannot stamp a parent it never submitted — so this summary is the only record. All targets refused raises instead of exiting 0. |
| **A second `sapia run boltzgen` into the same out_dir, after the first has finished.** | Supported: each group's attestation lives in `runs/<group>/.submit.json`, **beside the designs it describes**, so a resubmit replaces that group's record and touches no other. A group resubmitted *without* `--binding` correctly comes back with no epitope — it does not inherit the previous one. |
| **A second submit into the same leaf while the first is still RUNNING.** | **Not safe. Wait, or use `-l/--dir-label`.** The two submits share `specs/<group>.yaml` and the manifest file, and the staleness check reads directory state at submit time while `intermediate_designs_inverse_folded/` does not appear until BoltzGen's inverse-folding step runs minutes to hours later — so two submits close together can both pass it before either has designs on disk. |
| **`already holds designs but its .submit.json is …` / `… made from a DIFFERENT spec`.** | A resubmit into a group dir that already has outputs, where the record does not *prove* the same spec. BoltzGen is never given `--overwrite`, so those designs would survive and collect under this run's epitope. **Delete the WHOLE `runs/<group>/` directory, or use `-l/--dir-label`.** Deleting `.submit.json` alone leaves the old designs in place and is exactly what the check exists to catch — the refusal requires a readable *match*, not merely the absence of a known mismatch. |
| **`no .submit.json beside runs/<group>`.** | Every design of that group is an `error:` row. **The usual cause is a `-l/--dir-label` or `-t/--table` that resolves to a different out_dir than the run used** — check those before suspecting the outputs. |
| **`.submit.json was written for a different spec than specs/<group>.yaml`, or the two disagree about which steering was requested.** | The record and the spec BoltzGen consumed came from different runs. Every design of that group errors rather than publishing an attestation that does not describe them. |
| **The task dies with empty `.out` and `.err`.** | Container never started — almost always the Volume refusing to mount on a non-empty path. The image has a build-time guard against exactly this; if it fires you will see it in the *build* log instead, which is the point. |

---

## What it collects

Columns are leaf-prefixed `boltzgen_` (or `boltzgen<label>_` with `-l`). Trust columns first.

**P** marks a column that exists **only under `protein-anything`**; **N** one that is NA there.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `boltzgen_protocol` | `nanobody-anything` / `protein-anything`, from the group's `.submit.json` | **TRUST COLUMN, read FIRST.** On every row, including error rows. **NA = the record is missing/unreadable, or predates this column (the 2026-10-08/09 pilots) — never a claim that it was nanobody.** It decides the meaning of `scaffold`, `scaffold_src`, `cdr_sequence`, `fixed_resnums` and `liability_*`. Never compare those across two values of this column. |
| `boltzgen_status` | Per-design outcome | **`OK` is the only proof of success.** `error: …` carries its own diagnosis. |
| `boltzgen_mapping_seq_identity` | CSV sequence vs the binder chain of the structure `boltzgen_path` points at | `1.0` on every `OK` row by construction; the **measured** value is written on the identity-failure error row so you can see what condemned it, and is NA on error kinds that fail earlier. About the chain MAPPING — **not** sequence recovery against the parent framework. |
| `boltzgen_chain_map` | `ok` / `ok_renamed` / `no_source_context` | **Only the two `ok*` hold.** `no_source_context` means the whole `source_context` was null or unreadable — provenance LOST; `boltzgen_scaffold` is NA and the row should not be trusted. A per-chain missing `output_chain` cannot occur upstream and errors the row if it does. |
| `boltzgen_scaffold_src` | `in_run_scaffold_list` / `distinct_from_target` / **`not_applicable_no_scaffold`** / `not_in_run_scaffold_list` / `is_a_target_chain_source` / `no_source_context` | **Read this before `boltzgen_scaffold`, and read `boltzgen_protocol` before this.** Only the first two carry a name; the last four mean `boltzgen_scaffold` is NA. `distinct_from_target` = no allow-list to check against (a `--spec` nanobody run). **`not_applicable_no_scaffold` = `protein-anything`: there is no framework, so nothing was checked** — it is a distinct verdict precisely so it cannot be read as `distinct_from_target`, which the un-widened code would have emitted beside the scaffold "name" `inline`. |
| `boltzgen_n_design_res` | Designed residue count — CDRs (nanobody) or the **whole chain** (protein) | Agreed by three independent witnesses or the row errors — and an **absent** third witness errors too. **Under `protein-anything` a FOURTH witness is enforced: it must equal `boltzgen_binder_len`, and an inequality is an error row.** The other three all count designed residues, so they agree with each other and would miss a mask that does not describe the same chain; this one is the only check that catches it, and it is free. A partial mask is of course normal and correct under `nanobody-anything` — that is what a CDR/framework split *is* — so the check does not run there. **A fourth-witness error row keeps `boltzgen_mapping_seq_identity`, `boltzgen_chain_map` and `boltzgen_scaffold_src`**, which were measured before the check fired: read them first to tell "the mapping was sound and only the count disagreed" from "the whole mapping is adrift". |
| `boltzgen_binding_resnames` | `A91=LYS,A128=GLU,…` at submit time | Eyeball against your reference. NA when no `--binding`; if one was given and the record is missing, the designs are **errors**, so NA always means "none asked for". |
| `boltzgen_binding_frame` | `declared_sequence` / `synthesized_from_model` (per chain, `+`-joined if mixed) | Which RESIDUE frame the epitope was read in. `synthesized_from_model` means the file declared no polymer sequence, so positions are the modelled ordinals — fine for a predictor's output, a smell for a deposited structure. NA when no `--binding`. |
| `boltzgen_binding_chain_frame` | `label_asym_id` (mmCIF) / `author_letter` (PDB) | Which CHAIN identifier BoltzGen resolved against. NA when no `--binding`. |
| `boltzgen_chain_id_map` | `H>A,I>B` — author letter > the id written into the YAML | **Read it.** On an mmCIF the two routinely differ, and a wrong translation moves the whole epitope. Identity on a PDB. NA when no `--binding`. |
| `boltzgen_not_binding_resnames` | The ANTI-epitope, same `A91=LYS` notation | A separate column from `binding_resnames` on purpose: opposite claims must not share one. NA when no `--not-binding`. |
| `boltzgen_path` | The **refolded complex**, as PDB | Target **+** binder. Empty on an error row. Under `nanobody-anything` this is the **only** structure — a binder-alone one needs `chainsel`. |
| `boltzgen_cif_path` | The same complex as the original mmCIF | Present on an error row **when the file exists**, so a failure is inspectable; empty when it does not, so a downstream `-i boltzgen_cif_path` cannot pick up a path that is not there. |
| **P** `boltzgen_binder_refold_cif` | The **binder folded ALONE**, no target — `refold_design_cif/<id>.cif`, as mmCIF | **Carries a path only under `protein-anything`** (the `design_folding` step; issue #148 means the file is never written for nanobodies, so the column is **NA** there — present, like every documented column, but empty). B-factors are per-atom pLDDT. Feed it to `usalign` against the designed backbone, or to `pyrosetta` for a monomer score. Its absence on a `protein-anything` row is an **error row**, never a quietly missing column. **Same policy as `boltzgen_cif_path`, and the same trap: it stays populated on a MAPPING-error row when the file exists, so a failure is inspectable — while `boltzgen_path` goes empty.** `filter_ready` keys on a column being non-empty rather than on `<leaf>_status`, so **`-i boltzgen_binder_refold_cif` will pick up refused designs.** Filter on `boltzgen_status == "OK"` first, or use `boltzgen_path`, which is empty on every error row by design. |
| `boltzgen_binder_chain` | PDB chain label of the binder in `boltzgen_path` | What you pass to `chainsel --chains` / `cms --binder-chains` / `epitope --design-chains`. **Measured off the written structure's design mask, not off the requested id — and the two genuinely differ, so this column is the only safe source.** See the chain-id measurement below. |

### MEASURED 2026-10-09: BoltzGen discards the requested binder chain id

Verified on `outputs/20261009_111006_boltzgen_protein_verification`, a `protein-anything` run
against a 3-chain target (A/B/C) where the tool auto-picked and requested **`Z`**:

| Layer | Value |
| --- | --- |
| written spec, `protein:` entity `id:` | **`Z`** — we did ask for it |
| `source_context` designed entry, `source_chain` | **`Z`** — the request survives here |
| `source_context` designed entry, `output_chain` | **`D1`** |
| design mmCIF and refolded complex, chain labels | **A, B, C, D** — no `Z` anywhere |
| `boltzgen_binder_chain` (collected) | **`D`** ×10 |
| `boltzgen_chain_map` | **`ok`** ×10 — the link held; this is not a mapping failure |

The four `output_chain` values were `A1, B1, C1, D1` — **sequential by entity order**, not
collision-avoidance. No `Renaming`/`conflict` message appeared in `.out`, `.err` or
`spec_check` (unlike a nanobody run, which does print `Chain ids conflict … Renaming with
{'A': 'D'}` when the scaffold's chain genuinely collides with the target's).

**What this means.** BoltzGen appears to assign output chain ids positionally at parse time and
discard the requested id. The binder is the last entity, so on a 3-chain target it lands on `D`
— which is *also* where the nanobody protocol lands from a scaffold chain `A`, by a different
route. **Not proven**: upstream source was not read, and a `--binder-chain Q` run or a 4-chain
target would test the positional rule. What *is* measured is that the requested id does not
reach the output.

**Consequences, in order of how much they will cost you:**

1. **Never hard-code or predict the binder chain.** `boltzgen_binder_chain` is the only source.
   Everything downstream — `chainsel --chains`, `cms --binder-chains`, `epitope --design-chains`
   — must take it from the column.
2. The collision refusal on `--binder-chain` guards a collision that upstream may resolve on its
   own. It is retained because it costs nothing and upstream behaviour is unverified.
3. `chain_map == ok` (not `ok_renamed`) on every row: the `D` is already present at the
   `source_context` layer, so this is **not** the mmCIF writer's re-pooling described above —
   it happens earlier.
| `boltzgen_sequence` | The designed chain's **full** sequence | What a structure predictor should fold. Under `protein-anything` the whole chain is designed, so this is also the designed sequence. |
| **N** `boltzgen_cdr_sequence` | Only the designed residues, concatenated | Diversity bookkeeping; not a foldable sequence. **NA under `protein-anything`**: the whole chain is designed there, so upstream's value equals `boltzgen_sequence` and a column named `cdr_sequence` would assert a CDR/framework split that does not exist. The column still exists (all-NA), so a filter naming it cannot `KeyError`. |
| **N** `boltzgen_scaffold` | Basename of the VHH scaffold `.cif` drawn for this design | The diversity axis. **Read `boltzgen_scaffold_src` first** — NA whenever the name could not be validated, and **always NA under `protein-anything`**. Group by it before concluding anything about frameworks. |
| `boltzgen_design_resnums` | Designed positions, prosapia syntax (`26:33,51:58,97:115`), **1-indexed within the binder chain** | The CDRs under `nanobody-anything`; **`1:binder_len` under `protein-anything`**, where everything is designed. |
| `boltzgen_fixed_resnums` | The complement, same syntax | **What `proteinmpnn --fixed-positions` consumes** for a CDR-only redesign. **EMPTY under `protein-anything`** — correct, not a bug: nothing is fixed, so there is nothing to hand proteinmpnn. |
| `boltzgen_binder_len` | Binder chain length | **Under `nanobody-anything`**: framework + sampled CDR lengths. **Under `protein-anything`**: the drawn length, within `--binder-length`. **Check `.nunique()` against `num_designs / diffusion_batch_size`** — one length is drawn per batch. |
| `boltzgen_id` | BoltzGen's own `id` | The join key into the raw output on disk. |
| `boltzgen_design_ptm` | pTM of the designed chain | Fold confidence of the nanobody alone. |
| `boltzgen_design_iptm` | iPTM **averaged over every chain pair** | **Diluted by chain count — do not gate on it.** |
| `boltzgen_design_to_target_iptm` | Design↔target iPTM, per pair | **The honest interface confidence.** |
| `boltzgen_design_iiptm` | Intra-interface iPTM | Secondary. |
| `boltzgen_design_residue_iptm`, `boltzgen_target_ptm`, `boltzgen_ptm`, `boltzgen_iptm`, `boltzgen_protein_iptm` | Boltz-2's other confidence terms | Diagnostic; same dilution caveat on the complex-wide ones. |
| `boltzgen_min_design_to_target_pae`, `boltzgen_min_interaction_pae` | Minimum interface PAE, Å | Lower is better. No bar is given here: the observed range on this tool is unknown until it has been run. |
| `boltzgen_design_ipsae_min`, `boltzgen_design_to_target_ipsae`, `boltzgen_target_to_design_ipsae` | ipSAE variants | **Measured 2026-10-09: identically `0.0` in all 60 rows of a real run. Present but carrying NO INFORMATION — zero variance. Do not rank or filter on them.** They are in the CSV (so they do not trip the dropped-column trap) and they are collected, but they are not computed to anything useful. Why is **unknown**. Upstream's own ipSAE ranking column is `esmfold2_ipsae_min`, which needs the ESMFold2 step this tool never runs. |
| `boltzgen_complex_plddt`, `boltzgen_complex_iplddt`, `boltzgen_complex_pde`, `boltzgen_complex_ipde` | Complex-level confidence | Note: upstream issue #202 claims pLDDT is missing from BoltzGen metrics; on the pinned SHA `complex_plddt` **is** in `eval_keys_confidence` and is collected. The refolded mmCIF also carries per-atom pLDDT in its B-factor column. |
| `boltzgen_bb_rmsd_target` | Backbone RMSD of the **target** between design and refold, Å | **Gate 1.** The column that catches issue #226. |
| `boltzgen_bb_rmsd_design` | Backbone RMSD of the **design** | Self-consistency of the binder's own fold. |
| `boltzgen_bb_target_aligned_rmsd_design` | Design RMSD after aligning on the target, Å | The **pose** test, stricter than `bb_rmsd_design`. |
| `boltzgen_bb_rmsd`, `boltzgen_bb_rmsd_design_target` | Whole-complex and design+target backbone RMSDs | Context for the two above. |
| `boltzgen_delta_sasa_refolded` | The **target's** buried area on the refolded complex, Å², 1.4 Å probe | **ONE-SIDED BSA, not summed over both partners** (`get_delta_sasa`: `target_area - target_bound`). Roughly half a both-sides number like `polyv_dsasa`; not comparable with it or with pyrosetta's `if_dSASA`. Meaningless until gate 1 clears. |
| `boltzgen_design_sasa_bound_refolded`, `boltzgen_design_sasa_unbound_refolded` | The two terms above, **of the TARGET** | Upstream's `design_` prefix is a misnomer — these are the target's bound and unbound SASA. Collected under upstream's spelling so they match the raw CSV. |
| `boltzgen_plip_hbonds_refolded`, `boltzgen_plip_saltbridge_refolded` | Interface contact counts | Counts, not energies. **There is no `plip_hydrophobic_refolded`** — `count_noncovalents` (`analyze_utils.py:183-225`) returns exactly these two. This skill used to list a third; it never existed upstream and was removed from the declared set 2026-10-09. |
| `boltzgen_design_hydrophobicity`, `boltzgen_design_chain_hydrophobicity` | Hydrophobicity of the designed residues / the whole chain | Developability. Under `protein-anything` the two are the same quantity (everything is designed). |
| **P** `boltzgen_design_largest_hydrophobic_patch_refolded` | Largest contiguous hydrophobic surface patch on the refolded binder, Å² | **`protein-anything` only.** `nanobody-anything` sets `largest_hydrophobic_refolded=false` (`cli/boltzgen.py:99-105`); `protein-anything`'s config is `{}` so `analysis.yaml`'s `true` stands. Higher = more aggregation-prone. **No observed range — this has never been run.** |
| `boltzgen_liability_score`, `boltzgen_liability_num_violations`, `boltzgen_liability_details` | Sequence-motif developability liabilities | **The rule set DEPENDS ON THE PROTOCOL and the two are NOT comparable. Read `boltzgen_protocol` first.** `nanobody-anything` → **antibody** rules (`cli/boltzgen.py:99-105`; PyPI 0.3.2 would silently apply peptide rules, #228, which is one reason the image pins a git SHA). `protein-anything` → the base default **`liability_modality: peptide`** (`analysis.yaml`), whose motif set is `ProtTryp` on **every K/R**, an N-terminal cyclization motif and a DPP4 cleavage motif — **rules for a short synthetic peptide. On a folded 80-225 aa minibinder this is close to a lysine+arginine count and should not gate anything.** Turn it off with `--config-set 'analysis liability_analysis=false'` if you want the noise gone. |
| `boltzgen_num_design`, `boltzgen_num_tokens`, `boltzgen_num_prot_tokens` | Upstream's own counts | `num_design` is the third witness behind `boltzgen_n_design_res`. |
| `boltzgen_helix`, `boltzgen_sheet`, `boltzgen_loop` | Secondary-structure fractions | A VHH is β-rich; a helix-heavy nanobody is suspicious. Under `protein-anything` there is no such expectation — a de-novo minibinder may legitimately be all-helical. |
| **P** `boltzgen_designfolding-bb_rmsd_design` | Backbone RMSD between the designed binder and **the same sequence folded ALONE**, Å | **`protein-anything` only, and it is the self-consistency number** — what `nanobody-anything` needs `chainsel` + `usalign` to get. Lower = the binder folds to its designed shape without the target holding it there. **Note the HYPHEN** — it is upstream's own key (`analyze.py:1037-1039`), kept so the column matches the raw CSV. Address it as `df["boltzgen_designfolding-bb_rmsd_design"]`; it **cannot** be used in a `df.query()` string. |
| **P** `boltzgen_designfolding-bb_rmsd` | Whole-structure backbone RMSD of the same comparison | The binder alone *is* the whole structure here, so this should track `-bb_rmsd_design`. A disagreement between the two is a reason to look at the files. |
| **P** `boltzgen_designfolding-design_ptm`, `boltzgen_designfolding-complex_plddt` | Confidence of the binder folded **alone** | A monomer prediction: no target, no interface. Low here with high `design_to_target_iptm` means a binder that only folds when bound. |

> **The `designfolding-*` and `design_largest_hydrophobic_patch_refolded` names are read off
> `analyze.py` at the pinned SHA and have NEVER been observed in a real CSV.** If a name is
> wrong the column is **dropped**, not NA. The verification run must check for them by name;
> `--all-metrics` will show what is really there.

**Spec validation.** Every task runs `boltzgen check` on its written spec before `boltzgen
run`, and a non-zero exit kills the task. The colour-coded mmCIF it produces lands in
`<out_dir>/runs/<group>/spec_check/` — open it to see the binding site BoltzGen understood.
Protocol-independent.

**On disk**, per target group:
`<out_dir>/runs/<group>/intermediate_designs_inverse_folded/` holding `<id>.cif` (the design,
with the design mask in its B-factors and **zeroed side chains on designed residues**),
`<id>.npz` (`design_mask`, `source_context`), `refold_cif/<id>.cif` (the refolded complex),
`refold_design_cif/<id>.cif` (**the binder alone — `protein-anything` only**) and
`aggregate_metrics_analyze.csv`. The written spec YAMLs are under `<out_dir>/specs/`, and each
group's `.submit.json` (epitope attestation, protocol, step list, scaffold allow-list) sits at
`<out_dir>/runs/<group>/.submit.json`.
The join is by name — `file_name == id + ".cif"` — never by row order.

> **`name` COLLIDES ACROSS SIBLING BOLTZGEN TABLES.** BoltzGen derives each design's `id` from
> the target stem plus an index, so two shards of the same target produce the **same** names.
> Measured 2026-10-09: six sibling tables of one pilot had 10 distinct names each, every one
> repeated 6×. **Never `pd.merge(..., on="name")` across sibling tables** — it fans out to a
> cross product. Merge on `(table, name)`, or give each shard its own `--table-label` and keep
> them apart.

---

## Filters a caller can write

**No numeric cutoff appears below, on purpose.** The only thing filtered here is trust — which
is not a scientific choice but a statement about whether a row's chain labels and residue
numbers mean anything. The scientific gating is left to you, with your own run's distribution
in front of you; see the gate section above for the order and the direction.

```python
# -f module: keep only rows whose mapping was verified. Gates nothing scientific.
# This is the filter to apply BEFORE looking at any number, every time.
#
# PROTOCOL-AWARE. Under protein-anything scaffold_src is ALWAYS
# "not_applicable_no_scaffold" (there is no framework to validate), so the
# nanobody-only list would reject every minibinder row.
def apply_filter(df):
    # .isin(), not .str.startswith(): an error row leaves these NA, .str.* propagates
    # NA, and `&` then raises "Cannot mask with non-boolean array containing NA / NaN
    # values" -- so a startswith filter breaks on exactly the run you most need it for.
    scaffold_ok = [
        "in_run_scaffold_list",       # nanobody, validated against the run's list
        "distinct_from_target",       # nanobody --spec run, no list to check
        "not_applicable_no_scaffold", # protein-anything: nothing to check
    ]
    return df[
        (df["boltzgen_status"] == "OK")
        & (df["boltzgen_chain_map"].isin(["ok", "ok_renamed"]))
        & (df["boltzgen_scaffold_src"].isin(scaffold_ok))
    ]
```

```python
# Say WHICH binder you have, and refuse to mix. Run this before anything else on a
# run_dir that could hold both protocols.
def apply_filter(df):
    print(df["boltzgen_protocol"].value_counts(dropna=False))
    # NA here means the group's .submit.json is missing or predates the column --
    # it is NOT a claim that the protocol was nanobody.
    return df[df["boltzgen_protocol"] == "protein-anything"]
```

```python
# Since 2026-10-09 every documented column EXISTS in a freshly collected table, so a
# filter can name one without guarding. What it must still guard is an ALL-NA column:
# NA is not a number and `<` on it is False, so a naive comparison silently selects
# zero rows and looks like a failed gate rather than a missing metric.
def apply_filter(df):
    col = "boltzgen_designfolding-bb_rmsd_design"   # protein-anything only
    if col not in df.columns:
        # A table collected BEFORE 2026-10-09, or a column this skill does not
        # document. Old tables are not rewritten.
        raise SystemExit(f"{col} absent -- this table predates the pd.NA change")
    if df[col].isna().all():
        # Two causes, and you cannot tell them apart from here: the metric was not
        # produced, or the name declared by the tool is wrong. Re-collect with
        # --all-metrics and compare against the CSV's real column names.
        raise SystemExit(f"{col} is ALL-NA -- re-collect with --all-metrics")
    return df[(df["boltzgen_status"] == "OK") & (df[col] < 2.0)]
```

```python
# Then look at the distributions and pick YOUR cutoffs, in gate order.
# Nothing here has been run on real data, so no expected range is known.
def apply_filter(df):
    ok = df[df["boltzgen_status"] == "OK"]
    for col in ("boltzgen_bb_rmsd_target",            # gate 1, lower better
                "boltzgen_design_to_target_iptm",     # gate 2, higher better (NOT design_iptm)
                "boltzgen_min_interaction_pae",       # gate 2, lower better
                "boltzgen_bb_target_aligned_rmsd_design",  # gate 3, lower better
                "boltzgen_delta_sasa_refolded"):      # gate 4, higher better, ONE-SIDED
        print(col, ok[col].describe())
    return ok
```

Write the cutoffs you choose into the campaign log, with the run they were read off. They are a
property of your target and your batch, not of this tool.

```python
# NANOBODY diversity bookkeeping: did the random scaffold draw actually spread?
# scaffold_src first -- an unvalidated scaffold name is NA, and a value_counts()
# over a column full of NA would report a spread that is not there.
# Meaningless under protein-anything; scaffold is NA for every row there.
def apply_filter(df):
    ok = df[(df["boltzgen_status"] == "OK")
            & (df["boltzgen_scaffold_src"] == "in_run_scaffold_list")]
    print(ok["boltzgen_scaffold"].value_counts())
    return ok
```

```python
# MINIBINDER length bookkeeping: did the length range actually get sampled, or did
# the per-BATCH draw collapse it? Expect about num_designs / diffusion_batch_size
# distinct values; with --diffusion-batch-size 1, about num_designs.
# Record the number you see in the campaign log -- it is a property of YOUR run.
def apply_filter(df):
    ok = df[(df["boltzgen_status"] == "OK")
            & (df["boltzgen_protocol"] == "protein-anything")]
    print(ok["boltzgen_binder_len"].describe())
    print("distinct lengths:", ok["boltzgen_binder_len"].nunique(), "of", len(ok))
    return ok
```

### Handing CDRs to ProteinMPNN (nanobody-anything only)

`boltzgen_fixed_resnums` is written in exactly the syntax `proteinmpnn --fixed-positions`
parses (`,` fragments, `start:end` inclusive, 1-indexed within the chain), so a CDR-only
redesign is `--chains-to-design <binder chain> --fixed-positions <that string>`.

> **Trap: `{expr}` cannot inject a whole string.** prosapia's expression language resolves
> only **integers** (`expr.py` raises on a non-integral column), so
> `--fixed-positions '{boltzgen_fixed_resnums}'` will not work. CDR boundaries vary per design
> because CDR *lengths* are sampled. Either group designs that share a `fixed_resnums` string
> and run proteinmpnn once per group, or fix the whole target chain instead and let the binder
> be fully redesigned.

**Under `protein-anything` this section does not apply**: `boltzgen_fixed_resnums` is empty
because the whole chain is designed, so there is nothing to fix. If you want to re-design a
minibinder's sequence on its backbone, feed `boltzgen_binder_refold_cif` (or a `chainsel` cut
of `boltzgen_path`) to `proteinmpnn`/`atomium` with no `--fixed-positions` at all.

---

## What it is blind to

* **Whether the complex is real.** Every number here is Boltz-2's opinion of its own output.
  Fold the designed sequence with an *independent* predictor (`boltz`, `alphafold3`) and
  compare back — that is self-consistency against a *different* model, and this tool cannot
  provide it under either protocol. `designfolding-bb_rmsd_design` is still Boltz-2 checking
  Boltz-2.
* **A binder-alone structure exists under `protein-anything` ONLY.** `design_folding` runs for
  `protein-anything` and `protein-small_molecule` only (`cli/boltzgen.py:1237`), so under
  `nanobody-anything` `refold_design_cif/` is **never written** (issue #148) and a binder-alone
  check needs **`chainsel`** to cut `boltzgen_binder_chain` out of `boltzgen_path` first, then
  `usalign`. Under `protein-anything` use `boltzgen_binder_refold_cif` directly.
* **Interface energy, packing, shape complementarity.** `delta_sasa_refolded` is burial, not
  energy. Use **`pyrosetta`** (`if_dG`, `if_delta_unsat`, `packstat`) and **`cms`** (CMS, SC).
* **Where it actually bound.** `binding:` is a soft preference (#215), so a design can land
  anywhere. Use **`epitope`** on `boltzgen_path` with `--design-chains boltzgen_binder_chain`
  to measure the real footprint and hotspot recall.
* **Whether the target protomers are correctly arranged relative to the full assembly.**
  `bb_rmsd_target` compares the design to its own refold, not to your reference assembly. For
  an oligomeric target use **`usalign`** against the reference, and **`ringfit`** if the binder
  is meant to straddle two protomers.
* **Expression, aggregation, immunogenicity.** `liability_*` is a sequence-motif panel, not an
  assay — and under `protein-anything` it is a *peptide* motif panel applied to a folded
  protein, so it is not even the right panel. `design_largest_hydrophobic_patch_refolded`
  (protein-anything only) is the closer aggregation proxy, and it is still a surface area.
* **Whether a minibinder is well-behaved as a protein.** `protein-anything` produces de-novo
  chains with no framework guaranteeing solubility, expression or a folded core. Score them
  with **`pyrosetta`** (monomer `total_score`, `score_per_res`, `packstat`, `buried_unsat`) on
  `boltzgen_binder_refold_cif`, and fold the sequence independently with **`boltz`**.
* **Anything about a binder that is not ONE designed chain.** Both protocols assume exactly
  one; nothing here generalises to Fabs, peptides, small-molecule binders or multi-chain
  binders. Those are different BoltzGen protocols this tool does not wrap.

---

## Changes — 2026-10-09, the `protein-anything` widening

The tool was widened **in place**, keeping `name = "boltzgen"` and the `boltzgen_` column
prefix. Existing tables stay readable and the nanobody path is unchanged.

**Nothing about `nanobody-anything` changed.** The default protocol, its step list, its spec
shape and every column it writes are what they were. An invocation written before this date
behaves identically.

**New columns.** `boltzgen_protocol` (all rows, both protocols);
`boltzgen_binder_refold_cif`, `boltzgen_designfolding-bb_rmsd_design`,
`boltzgen_designfolding-bb_rmsd`, `boltzgen_designfolding-design_ptm`,
`boltzgen_designfolding-complex_plddt`, `boltzgen_design_largest_hydrophobic_patch_refolded`
(`protein-anything` only).

**Columns whose MEANING now depends on `boltzgen_protocol`** — rows collected under one value
are **not comparable** with rows collected under the other: `boltzgen_scaffold`,
`boltzgen_scaffold_src`, `boltzgen_cdr_sequence`, `boltzgen_design_resnums`,
`boltzgen_fixed_resnums`, `boltzgen_n_design_res`, `boltzgen_liability_score`,
`boltzgen_liability_num_violations`, `boltzgen_liability_details`.

**`boltzgen_scaffold_src` has a new value**, `not_applicable_no_scaffold`. A filter written
before this date that allow-lists `["in_run_scaffold_list", "distinct_from_target"]` will
reject **every** `protein-anything` row. Update it.

**The scaffold check is now FAIL-CLOSED.** It runs only when the record explicitly says
`nanobody-anything`; `protein-anything`, an absent protocol and an unrecognised string all get
`not_applicable_no_scaffold` + NA. The previous form tested `protocol != protein-anything`, so
an absent or unknown protocol fell into the branch that *publishes* a provenance verdict —
a de-novo design would have come back with scaffold `"inline"` and `distinct_from_target`,
which reads as *checked and passed*. An unrecognised (not merely absent) protocol string now
errors every design of its group rather than being read under the nearest-looking rules.

**The resubmit guard now compares the protocol and step list, not just the spec hash.** Under
`--spec` the YAML is byte-identical whichever protocol was passed, so the hash alone let a
second submit relabel the first's designs. A record with no `protocol`/`steps` — anything
written before this date — can no longer be topped up; use `-l/--dir-label`.

**`--extra-flags` additionally refuses `--overwrite`, `--reuse` and `--diffusion_batch_size`.**

**`protein-anything` gains a fourth count witness**: `n_design_res` must equal `binder_len`,
and an inequality is an error row — one that **keeps the mapping measurements**
(`mapping_seq_identity`, `chain_map`, `scaffold_src`) so the failure is diagnosable from the
table without a re-run.

**The structural, trust and attestation columns are now seeded on every row**, so every column
this skill documents exists on every row of every table the tool writes — including a collect
where every design errored. That is the shape the mandatory trust filter is most needed for,
and the one it used to raise `KeyError` on.

**A zero-row metrics CSV no longer prints `OK`.** It prints `NO DESIGNS` and names
`design_folding` as the `protein-anything` cause.

**Three declared columns were removed**, because upstream never writes any of them (see the
column table for the line numbers): `boltzgen_plip_hydrophobic_refolded`, `boltzgen_target_id`,
`boltzgen_free_cysteines`. No table ever carried them. This matters more after the `pd.NA`
change than before it: left declared, they would now appear in every table as permanently
all-NA columns, which is exactly the signature of a misspelled name.

**`boltzgen_protocol` is NA on anything collected before this change** — the two pre-column
run_dirs are named at the top of this file. That NA is expected, not a failed attestation.

**Defect FIXED: a declared-but-absent metric is now `pd.NA`, not dropped.** It used to be
omitted from the row entirely, so `prosapia` never minted the column and a `-f` filter naming
one raised `KeyError` instead of filtering. Now every column this skill documents exists in
every freshly collected table — **the table's shape is a function of the tool version, not of
what one CSV happened to contain** — which is the precondition for a filter being writable at
all. Additive: existing tables are untouched and are **not** rewritten, so a table collected
before this date still lacks those columns and still needs `if col in df.columns`.

**The sharp edge that fix introduces, and what replaces the lost signal.** The `KeyError` was
loud; an all-NA column is not. A **misspelled declared name** now presents as a plausible
all-NA column rather than an error — and the five `protein-anything` names have never been
seen in a real CSV. Two defences:

* the collect prints one line **naming every declared column that came back absent**
  (`not in this CSV, so written as NA …`). Read it. It is the only live detector;
* an all-NA column has **two** causes — not produced, or misspelled — and **`--all-metrics`,
  which takes the CSV's own column names, is what tells them apart.** Never conclude "not
  measured" from NA alone.

**Not yet verified on real data.** `protein-anything` has never been run. Until it has, treat
as unknown: whether `designfolding-*` and `design_largest_hydrophobic_patch_refolded` appear
in the CSV under those exact names; whether `refold_design_cif/<id>.cif` is written with the
design `id` as its stem; whether a de-novo `protein:` entity's `source_context` entry really
carries `source: "inline"` and an `output_chain`; whether `boltzgen check` accepts the written
`protein:` entity; and **whether one length really is drawn per batch** (count
`boltzgen_binder_len.nunique()` on the first run with a known `--diffusion-batch-size`).
