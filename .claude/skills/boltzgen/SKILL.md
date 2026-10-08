---
name: boltzgen
description: >-
  How to design NANOBODY (VHH) binders against a protein target with the custom boltzgen tool on
  Modal — new CDRs grafted onto a fixed, randomly drawn VHH framework, one BoltzGen run per target,
  designs as child rows. Load before composing a boltzgen run, and before reading any number it
  wrote: design_iptm is DILUTED BY CHAIN COUNT (against a homotrimer it is a 4-chain average and
  reads low for a perfectly good binder) — design_to_target_iptm is the honest per-pair number.
  Covers why target geometry (bb_rmsd_target) is the FIRST gate and not the interface, the
  mapping-trust columns that must be read before any sequence or chain label, the forced step list
  that exists because `filtering` crashes and cannot be turned off, the TWO undocumented identifier
  frames an epitope depends on (residue numbers are label_seq_id not auth_seq_id; chain ids are
  label_asym_id for mmCIF but the author letter for PDB — on 1TUP author E,F,A,B are label A,B,C,D),
  and the fact that there is NO binder-alone structure so self-consistency needs chainsel first.
---

# boltzgen

**What it does.** Runs BoltzGen's `nanobody-anything` protocol: it takes a protein target
(one chain or many), picks one of several VHH frameworks **at random per design**, redesigns
CDR1/2/3 in sequence, structure *and length*, inverse-folds the result with SolubleMPNN,
refolds the complex with Boltz-2, and analyses it. `action: create` — one BoltzGen run per
target produces N designs, so the row count is known only at collect time. `-t` is optional:
with it, one run per ready row of that table; without it, a root run that starts `table0`.

**Modal only.** There is no `SAPIA_ACTIVATE_BOLTZGEN`; the cluster has no BoltzGen install and
the ~6.4 GB of checkpoints are staged on a Modal Volume. Run it with `-e modal`.

---

## Premise, and when its numbers are meaningless

| Assumption | What happens when it does not hold |
| --- | --- |
| **The binder is ONE designed chain** — a nanobody, framework fixed in sequence *and* structure, only the CDRs free. | The collector refuses the design (`boltzgen_status` = `error: MappingError: expected exactly one designed chain …`; every mapping refusal is prefixed `error: <ExceptionType>: `). It never guesses which chain a sequence belongs to. |
| **The target is a real structure**, in whatever frame you handed it. | BoltzGen renumbers nothing you can rely on; residue indices you supply are `label_seq_id` (see below). |
| **The target stays where you put it.** | **It may not.** Upstream issue #226 (open, unanswered): on a multi-chain target a second protomer can MOVE during design, so the complex has the wrong binding geometry. `boltzgen_bb_rmsd_target` is the only thing that catches it. |
| **Designs are not reproducible.** BoltzGen has no `--seed` flag at all, and the framework is `random.choice`d per design (`data/parse/schema.py:1947-1949`). | This is **deliberate** — scaffold diversity is the point of the protocol, and `boltzgen_scaffold` records what each design actually used, which is a stronger guarantee than a seed. Do not "fix" it by pinning `--scaffolds` to one name unless you specifically want one framework. |

This tool does **not** do de-novo backbone design (that is `rfdiffusion3` / `bindcraft2`), does
not design Fabs or peptides, does not run affinity prediction, and does not rank anything.

### The step list is fixed, and that is not an oversight

Every run is `--steps design inverse_folding folding analysis` plus
`--config analysis esmfold2_metrics=false`. There is no `--steps` flag on this tool. Three
facts, all verified against the pinned SHA, force exactly that list:

* generation never touches ESMFold2 — `task/predict/` contains no reference to it;
* `analysis` is launched with `esmfold2_metrics=true` and `analyze.py:638-652` then reads
  `esmfold2_scores/<id>.input.json` unconditionally → **it crashes** unless that key is false;
* `filtering` sets `self.use_esmfold2 = not use_affinity` at `filter.py:218` — **hard-wired to
  the protocol, reachable by no config key** — and `filter.py:421` raises
  `ValueError("Missing ESMFold2 scores; run esmfold2_scoring and analysis before filtering")`.
  **It can only be omitted.**

> **Trap: anyone running BoltzGen with default steps gets a crash at the last step.** The
> signature is that `intermediate_designs_inverse_folded/` is full and
> `final_ranked_designs/` is absent, with that ValueError in the task `.err`.
> `--esmfold2_acceleration off` does **not** help — it selects a numerical backend, it does
> not disable the model.

Consequences: no `final_ranked_designs/`, no `results_overview.pdf`, no `--budget`, and no
`esmfold2_ipsae_min` (upstream's own ranking column for polymer protocols). **Ranking happens
in the table with `-f`.** That is the right place for it anyway.

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
| `boltzgen_status` | **The only proof of success.** Anything but `OK` means no chain label, sequence, scaffold or residue list was written for that row — only the metrics, which are true of *some* chain pair but not of a named one. | `OK` |
| `boltzgen_mapping_seq_identity` | Identity between the CSV's `designed_chain_sequence` and the residues of chain `boltzgen_binder_chain` in the structure `boltzgen_path` points at. **This is about the chain MAPPING, not about sequence recovery against the parent framework** — hence the name. On an `OK` row it is always exactly `1.0`, because anything lower is refused; on the **identity-failure** error row the measured value is still written, so the number that condemned the design is visible. Other error kinds (a count mismatch, a corrupt file) leave it NA — they fail before an identity is computed. | `1.0` on `OK`; read it on `error` rows |
| `boltzgen_chain_map` | How the binder chain was tied to its `source_context` entry. **Only `ok` and `ok_renamed` mean the link holds.** `ok` — `output_chain`, digits stripped, equals the written label. `ok_renamed` — `output_chain` is present and differs, i.e. the mmCIF writer re-pooled the label; the link rests on chain **position**, which upstream enforces with `strict=True` (`source_context.py`), so it still holds. `no_source_context` — the whole `source_context` was JSON `null` or unreadable: **provenance is LOST, not benignly absent**; treat the row as untrusted and `boltzgen_scaffold` is NA. A per-chain missing `output_chain` cannot occur upstream (`update_designed` writes it for every entry, `source_context.py:200-205`) and is a hard **error** if it ever does. | `ok` or `ok_renamed` |
| `boltzgen_scaffold_src` | **How far `boltzgen_scaffold` itself was validated.** `in_run_scaffold_list` — the name matched a scaffold **this group** was submitted with (read from `runs/<group>/.submit.json`, so a later submit with a different `--scaffolds` cannot re-point it). `distinct_from_target` — no allow-list existed (a `--spec` run), but the source differs from every target chain's. `not_in_run_scaffold_list` / `is_a_target_chain_source` / `no_source_context` — **`boltzgen_scaffold` is NA**; the positional lookup produced something that is not a scaffold of this run, or is the target itself. | `in_run_scaffold_list` |
| `boltzgen_n_design_res` | Number of designed (CDR) residues. Cross-checked **three ways** — the design-mask B-factor of the design mmCIF, the npz `design_mask`, and the CSV's own `num_design`. **Any disagreement is an error row, never a trimmed fit — and so is an absent `num_design`**, because the three-witness guarantee would otherwise silently become a two-witness one. | never 0; observed range unknown until a real run |
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
> front of them. **Nothing here has been run on real data, so no observed range for any column
> is known.** A threshold written into a skill is copy-pasted long after it has gone stale —
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
> zero designs here — but treat a conventional RMSD bar borrowed from `protein-anything` with
> suspicion on nanobodies.

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

**Gate 5 — is it developable? `boltzgen_liability_score`, `boltzgen_liability_num_violations`,
`boltzgen_free_cysteines`, `boltzgen_design_hydrophobicity`.** The pinned SHA sets
`liability_modality=antibody` for this protocol — **PyPI 0.3.2 does not** (issue #228, a PR
was made 2026-05-19 and never merged), which is one of the reasons this tool pins a git SHA.

---

## Invocation and flags

```bash
# root run: one BoltzGen run against a target not in any table yet
sapia run boltzgen <run_dir> -e modal \
    --target-structure inputs/trimer.cif --target-chains A,B,C \
    --binding A91,A128,B91,B128 --num-designs 50

# child run: one BoltzGen run per target already in a table
sapia run boltzgen <run_dir> -t table0 -i chainsel_path -e modal \
    --binding 'A{epitope_start}-{epitope_end}' --num-designs 200

sapia collect boltzgen <run_dir> -t <the table the run reserved>
```

| Flag | Default | Notes |
| --- | --- | --- |
| `--target-structure` | none | Root run only. With `-t` the target comes from `--input-column`. Group name is `<stem>_bg`. |
| `--target-chains` | `all` | prosapia chain mini-language (`,` separates, `:` is an inclusive letter range): `A`, `A,B,C`, `A:C`. **AUTHOR chain letters**, matched **case-sensitively** — the ones a viewer shows; the tool translates them (see the chain-ID trap). **An explicit list writes `include:` with the POLYMER subchains only, so any cofactor, metal or glycan on those chains is DROPPED — the default `all` keeps them.** The submit prints what it drops; the two modes design against materially different targets. |
| `--binding` | none | Epitope, chain-prefixed: **AUTHOR chain letter** + **1-based position in the DECLARED polymer sequence** (`label_seq_id`), e.g. `A91,A128,B91-95`. `{expr}` islands resolved per design up the lineage. Soft preference, not a constraint (#215). |
| `--not-binding` | none | Same syntax and frame. See #217 before putting both on one chain. |
| `--scaffolds` | `7eow,7xl0,gontivimab,isecarosmab,sonelokimab` | The VHH frameworks BoltzGen draws from, **at random, one per design**. Bare names resolve against `--scaffold-dir`; a name with `/` is a path. `8coh`, `8z8v` and `gefurulimab` also ship. One name pins the framework and kills the diversity. |
| `--scaffold-dir` | `$BOLTZGEN_SCAFFOLD_DIR` or `/opt/boltzgen/example/nanobody_scaffolds` | Inside the task container (Modal only). **Not checked at submit time** — it does not exist on the workstation. |
| `--num-designs` | **50** | **Upstream's default is 10000 and its advice for real work is 10,000–60,000.** 50 is upstream's own smoke-test size and is the default here *because per-design runtime on this image is unmeasured*. Raise it once a verification run has timed it. Issue #139: values above 10000 are rejected. |
| `--diffusion-batch-size` | upstream auto (1 below 100 designs, else 10) | **Every design in one batch shares the same sampled CDR length** — a large batch relative to `--num-designs` under-samples lengths. Same trap as rfd3's `--length min-max`. |
| `--num-workers` | `1` | → `--num_workers`. |
| `--spec` | none | A BoltzGen design YAML used verbatim (still `{expr}`-resolved). Mutually exclusive with the spec-building flags — passing both is refused, not silently ignored. **Use absolute paths inside it**: BoltzGen resolves relative paths against the *original* YAML's directory, not the resolved copy's. **No chain translation happens on a `--spec` YAML** — you write the ids BoltzGen resolves (`label_asym_id` for an mmCIF target), and `boltzgen_chain_id_map` is NA, so nothing in the table attests them. |
| `--extra-spec` | none | YAML/JSON mapping merged into every written spec. A key this tool also writes (`entities`) is an error. **`{expr}` is resolved but chain ids are NOT translated** — any chain id inside it reaches the YAML verbatim and is not covered by `boltzgen_chain_id_map`, so write the ids BoltzGen resolves (`label_asym_id` for an mmCIF target). |
| `--config-set` | none | Repeatable, `'<step> <key>=<value>'` → `--config <step> <key>=<value>`. Appended after the forced `analysis esmfold2_metrics=false`; **setting that key to anything but false is refused at submit time**. |
| `--extra-flags` | `""` | Verbatim `boltzgen run` tokens. **`--protocol`, `--steps`, `--output`, `--num_designs`, `--config` and `--cache` are refused at submit time** — they are the tool's contract and argparse is last-wins. |
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
| **A declared metric is absent from the CSV.** | One printed line per collect: `not in this CSV (left NA): …`, and those columns are `NA`. **NA is the absence of a measurement, never 0.** |
| **A design that cannot be mapped.** | `boltzgen_status` starts with `error:`, `boltzgen_path` is **empty**, `boltzgen_cif_path` still points at the refolded mmCIF so you can look, and every chain/sequence/scaffold/residue column is NA. The metrics are still there. |
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

| Column | Meaning | How to read it |
| --- | --- | --- |
| `boltzgen_status` | Per-design outcome | **`OK` is the only proof of success.** `error: …` carries its own diagnosis. |
| `boltzgen_mapping_seq_identity` | CSV sequence vs the binder chain of the structure `boltzgen_path` points at | `1.0` on every `OK` row by construction; the **measured** value is written on the identity-failure error row so you can see what condemned it, and is NA on error kinds that fail earlier. About the chain MAPPING — **not** sequence recovery against the parent framework. |
| `boltzgen_chain_map` | `ok` / `ok_renamed` / `no_source_context` | **Only the two `ok*` hold.** `no_source_context` means the whole `source_context` was null or unreadable — provenance LOST; `boltzgen_scaffold` is NA and the row should not be trusted. A per-chain missing `output_chain` cannot occur upstream and errors the row if it does. |
| `boltzgen_scaffold_src` | `in_run_scaffold_list` / `distinct_from_target` / `not_in_run_scaffold_list` / `is_a_target_chain_source` / `no_source_context` | **Read this before `boltzgen_scaffold`.** Only the first two carry a name; the last three mean `boltzgen_scaffold` is NA. `distinct_from_target` means there was no allow-list to check against (a `--spec` run). |
| `boltzgen_n_design_res` | Designed (CDR) residue count | Agreed by three independent witnesses or the row errors — and an **absent** third witness errors too. |
| `boltzgen_binding_resnames` | `A91=LYS,A128=GLU,…` at submit time | Eyeball against your reference. NA when no `--binding`; if one was given and the record is missing, the designs are **errors**, so NA always means "none asked for". |
| `boltzgen_binding_frame` | `declared_sequence` / `synthesized_from_model` (per chain, `+`-joined if mixed) | Which RESIDUE frame the epitope was read in. `synthesized_from_model` means the file declared no polymer sequence, so positions are the modelled ordinals — fine for a predictor's output, a smell for a deposited structure. NA when no `--binding`. |
| `boltzgen_binding_chain_frame` | `label_asym_id` (mmCIF) / `author_letter` (PDB) | Which CHAIN identifier BoltzGen resolved against. NA when no `--binding`. |
| `boltzgen_chain_id_map` | `H>A,I>B` — author letter > the id written into the YAML | **Read it.** On an mmCIF the two routinely differ, and a wrong translation moves the whole epitope. Identity on a PDB. NA when no `--binding`. |
| `boltzgen_not_binding_resnames` | The ANTI-epitope, same `A91=LYS` notation | A separate column from `binding_resnames` on purpose: opposite claims must not share one. NA when no `--not-binding`. |
| `boltzgen_path` | The **refolded complex**, as PDB | Target **+** binder. There is **no binder-alone structure** — see below. Empty on an error row. |
| `boltzgen_cif_path` | The same complex as the original mmCIF | Present on an error row **when the file exists**, so a failure is inspectable; empty when it does not, so a downstream `-i boltzgen_cif_path` cannot pick up a path that is not there. |
| `boltzgen_binder_chain` | PDB chain label of the nanobody in `boltzgen_path` | What you pass to `chainsel --chains` / `cms --binder-chains` / `epitope --design-chains`. |
| `boltzgen_sequence` | The designed chain's **full** sequence (framework + CDRs) | What a structure predictor should fold. |
| `boltzgen_cdr_sequence` | Only the designed residues, concatenated | Diversity bookkeeping; not a foldable sequence. |
| `boltzgen_scaffold` | Basename of the VHH scaffold `.cif` drawn for this design | The diversity axis. **Read `boltzgen_scaffold_src` first** — this column is NA whenever the name could not be validated. Group by it before concluding anything about frameworks. |
| `boltzgen_design_resnums` | Designed positions, prosapia syntax (`26:33,51:58,97:115`), **1-indexed within the binder chain** | The CDRs. |
| `boltzgen_fixed_resnums` | The complement, same syntax | **What `proteinmpnn --fixed-positions` consumes** for a CDR-only redesign. |
| `boltzgen_binder_len` | Binder chain length | CDR length is sampled, so this varies design to design. |
| `boltzgen_id` | BoltzGen's own `id` | The join key into the raw output on disk. |
| `boltzgen_design_ptm` | pTM of the designed chain | Fold confidence of the nanobody alone. |
| `boltzgen_design_iptm` | iPTM **averaged over every chain pair** | **Diluted by chain count — do not gate on it.** |
| `boltzgen_design_to_target_iptm` | Design↔target iPTM, per pair | **The honest interface confidence.** |
| `boltzgen_design_iiptm` | Intra-interface iPTM | Secondary. |
| `boltzgen_design_residue_iptm`, `boltzgen_target_ptm`, `boltzgen_ptm`, `boltzgen_iptm`, `boltzgen_protein_iptm` | Boltz-2's other confidence terms | Diagnostic; same dilution caveat on the complex-wide ones. |
| `boltzgen_min_design_to_target_pae`, `boltzgen_min_interaction_pae` | Minimum interface PAE, Å | Lower is better. No bar is given here: the observed range on this tool is unknown until it has been run. |
| `boltzgen_design_ipsae_min`, `boltzgen_design_to_target_ipsae`, `boltzgen_target_to_design_ipsae` | ipSAE variants | Present, but upstream's own ranking column (`esmfold2_ipsae_min`) is not — that needs ESMFold2. |
| `boltzgen_complex_plddt`, `boltzgen_complex_iplddt`, `boltzgen_complex_pde`, `boltzgen_complex_ipde` | Complex-level confidence | Note: upstream issue #202 claims pLDDT is missing from BoltzGen metrics; on the pinned SHA `complex_plddt` **is** in `eval_keys_confidence` and is collected. The refolded mmCIF also carries per-atom pLDDT in its B-factor column. |
| `boltzgen_bb_rmsd_target` | Backbone RMSD of the **target** between design and refold, Å | **Gate 1.** The column that catches issue #226. |
| `boltzgen_bb_rmsd_design` | Backbone RMSD of the **design** | Self-consistency of the binder's own fold. |
| `boltzgen_bb_target_aligned_rmsd_design` | Design RMSD after aligning on the target, Å | The **pose** test, stricter than `bb_rmsd_design`. |
| `boltzgen_bb_rmsd`, `boltzgen_bb_rmsd_design_target` | Whole-complex and design+target backbone RMSDs | Context for the two above. |
| `boltzgen_delta_sasa_refolded` | The **target's** buried area on the refolded complex, Å², 1.4 Å probe | **ONE-SIDED BSA, not summed over both partners** (`get_delta_sasa`: `target_area - target_bound`). Roughly half a both-sides number like `polyv_dsasa`; not comparable with it or with pyrosetta's `if_dSASA`. Meaningless until gate 1 clears. |
| `boltzgen_design_sasa_bound_refolded`, `boltzgen_design_sasa_unbound_refolded` | The two terms above, **of the TARGET** | Upstream's `design_` prefix is a misnomer — these are the target's bound and unbound SASA. Collected under upstream's spelling so they match the raw CSV. |
| `boltzgen_plip_hbonds_refolded`, `boltzgen_plip_saltbridge_refolded`, `boltzgen_plip_hydrophobic_refolded` | Interface contact counts | Counts, not energies. |
| `boltzgen_design_hydrophobicity`, `boltzgen_design_chain_hydrophobicity` | Hydrophobicity of the designed residues / the whole chain | Developability. |
| `boltzgen_liability_score`, `boltzgen_liability_num_violations`, `boltzgen_liability_details` | Antibody developability liabilities | **Antibody** rules on this SHA; PyPI 0.3.2 would silently apply *peptide* rules (#228). |
| `boltzgen_free_cysteines` | Unpaired cysteines in the design | The protocol already sets `filter_cysteine=true` for inverse folding. |
| `boltzgen_num_design`, `boltzgen_num_tokens`, `boltzgen_num_prot_tokens` | Upstream's own counts | `num_design` is the third witness behind `boltzgen_n_design_res`. |
| `boltzgen_helix`, `boltzgen_sheet`, `boltzgen_loop` | Secondary-structure fractions | A VHH is β-rich; a helix-heavy design is suspicious. |
| `boltzgen_target_id` | Upstream's regex-derived target id | Provenance only; derived from `id`, not independent. |

**Spec validation.** Every task runs `boltzgen check` on its written spec before `boltzgen
run`, and a non-zero exit kills the task. The colour-coded mmCIF it produces lands in
`<out_dir>/runs/<group>/spec_check/` — open it to see the binding site BoltzGen understood.

**On disk**, per target group:
`<out_dir>/runs/<group>/intermediate_designs_inverse_folded/` holding `<id>.cif` (the design,
with the design mask in its B-factors and **zeroed side chains on designed residues**),
`<id>.npz` (`design_mask`, `source_context`), `refold_cif/<id>.cif` (the refolded complex) and
`aggregate_metrics_analyze.csv`. The written spec YAMLs are under `<out_dir>/specs/`.
The join is by name — `file_name == id + ".cif"` — never by row order.

---

## Filters a caller can write

**No numeric cutoff appears below, on purpose.** The only thing filtered here is trust — which
is not a scientific choice but a statement about whether a row's chain labels and residue
numbers mean anything. The scientific gating is left to you, with your own run's distribution
in front of you; see the gate section above for the order and the direction.

```python
# -f module: keep only rows whose mapping was verified. Gates nothing scientific.
# This is the filter to apply BEFORE looking at any number, every time.
def apply_filter(df):
    # .isin(), not .str.startswith(): an error row leaves these NA, .str.* propagates
    # NA, and `&` then raises "Cannot mask with non-boolean array containing NA / NaN
    # values" -- so a startswith filter breaks on exactly the run you most need it for.
    return df[
        (df["boltzgen_status"] == "OK")
        & (df["boltzgen_chain_map"].isin(["ok", "ok_renamed"]))
        & (df["boltzgen_scaffold_src"].isin(
            ["in_run_scaffold_list", "distinct_from_target"]))
    ]
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
# Diversity bookkeeping: did the random scaffold draw actually spread?
# scaffold_src first -- an unvalidated scaffold name is NA, and a value_counts()
# over a column full of NA would report a spread that is not there.
def apply_filter(df):
    ok = df[(df["boltzgen_status"] == "OK")
            & (df["boltzgen_scaffold_src"] == "in_run_scaffold_list")]
    print(ok["boltzgen_scaffold"].value_counts())
    return ok
```

### Handing CDRs to ProteinMPNN

`boltzgen_fixed_resnums` is written in exactly the syntax `proteinmpnn --fixed-positions`
parses (`,` fragments, `start:end` inclusive, 1-indexed within the chain), so a CDR-only
redesign is `--chains-to-design <binder chain> --fixed-positions <that string>`.

> **Trap: `{expr}` cannot inject a whole string.** prosapia's expression language resolves
> only **integers** (`expr.py` raises on a non-integral column), so
> `--fixed-positions '{boltzgen_fixed_resnums}'` will not work. CDR boundaries vary per design
> because CDR *lengths* are sampled. Either group designs that share a `fixed_resnums` string
> and run proteinmpnn once per group, or fix the whole target chain instead and let the binder
> be fully redesigned.

---

## What it is blind to

* **Whether the complex is real.** Every number here is Boltz-2's opinion of its own output.
  Fold the designed sequence with an *independent* predictor (`boltz`, `alphafold3`) and
  compare back — that is self-consistency, and this tool cannot provide it.
* **There is no binder-alone structure.** `design_folding` is disabled for nanobody binders
  and `refold_design_cif/` is **never written** (issue #148). A binder-alone self-consistency
  check needs **`chainsel`** to cut `boltzgen_binder_chain` out of `boltzgen_path` first, then
  `usalign`.
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
  assay.
* **Anything about a non-nanobody binder.** The premise is one designed chain on a VHH
  framework; nothing here generalises to Fabs, peptides or de-novo backbones.
