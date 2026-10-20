---
name: alphafold3
description: >-
  How do I predict a binder-target complex with AlphaFold3 (v3.0.4), and which MSA, template and chain policy does each entity get? Per-chain MSA policy, templates, chain-letter rules and the columns it writes. Load before composing an alphafold3 run or interpreting its numbers. Traps: chain letters are positional in the `/` sequence and a wrong-letter MSA is not an error; changing MSA, template or chain policy needs a fresh `-l`; the tool validates nothing it stages.
---

# alphafold3

AlphaFold3 structure prediction over the sequences in a table. **`action: update`** — annotates the table it ran against, never mints a child. `-t/--table` is required. Default input column is `proteinmpnn_sequence` (a `/`-chainbreak sequence); override with `-i`.

Runs on both executors: Modal (own image, built from source) and vib (`SAPIA_ACTIVATE_ALPHAFOLD3` points at the lab's shared `alphafold3.sh`; `$AF3_DATABASE` is bound into the container).

## Premise, and when its numbers mean nothing

AF3 predicts the entities you hand it. The premise is that each entity's **sequence and its MSA describe the same molecule**, and that the chain letters in your flags are the chains you think.

- **Chain letters are positional in the `/`-chainbreak sequence**, after `--chains` selection: first segment is A, second B, and so on. The sequence carries no chain-letter information. **Which column you pass decides the letters** — `proteinmpnn_sequence` is the binder alone (one chain), `mkcomplex_sequence` is the whole complex. Count the `/` before writing any chain letter.
- A letter the complex does not have is an **error**. A letter that exists but is the wrong molecule is **not** — it will happily give the binder the target's MSA.
- **Chains sharing a sequence collapse into ONE entity** with a multi-letter id, so a per-chain policy is really per-**entity**. For a homo-oligomeric target that is correct: one MSA covers every copy. Name them together (`--msa-file A:C=target.a3m`); naming only some raises at submit. **Do not reach for `--chains` to silence that** — dropping a protomer predicts a *different complex*, still collects `OK`, and gives an `iptm` not comparable with the one you are ranking against.
- **It always returns numbers.** No input makes AF3 report that the question was meaningless. A confident prediction of the wrong complex looks exactly like a confident prediction of the right one.
- **Nothing checks that an A3M belongs to the chain you gave it.** `--msa-file` copies the file in unread. AF3's own `data/msa.py` compares the first record to the chain's sequence and raises at **run** time, so a mismatch costs a dead container per shard instead of a submit-time error. Read the first record yourself before a fan-out.
- **`--positions` slicing invalidates an MSA** for that chain: the sliced sequence is no longer the A3M's query. Nothing stops you at submit; every task dies in featurisation.
- **`--positions` with `--templates-json` is the silent version of that, and the worst thing here.** `queryIndices` index the *unsliced* chain. After slicing they are still in range, still equal in length to `templateIndices`, still schema-valid — so AF3 cannot catch it and nothing dies. You get a confident pose with the template register shifted by however many residues you removed. **Do not combine the two flags.** There is no honest automatic re-mapping, because `--positions` can concatenate non-contiguous fragments.
- **Every row must have the same chain layout, and nothing checks that.** `select_chains` *trims*: letters beyond a row's segment count are silently dropped, and with no `--chains` each row simply becomes however many chains its own sequence has. A pooled or merged table whose rows came through different parents can therefore be predicted as a 4-chain complex on some rows and a 2-chain one on others, every row `OK`, in one `iptm` column. **Because `iptm` is chain-count-diluted, that does not add noise — it inverts**: the rows that lost their target chains score *higher*. Check the `/` count across the whole input column before a fan-out, not on the first row.

It does **not** build MSAs. `--msa-file` takes an A3M you already have.

## Read these before any score

No numeric trust column exists — it is a predictor, not a mapper. The tool validates nothing it stages, so the evidence is **only** in what the run *wrote*:

| Where | Read it for |
| --- | --- |
| `<run_dir>/<table>/<leaf>/af3_inputs/<design>.json` | **The audit trail of what each chain was given.** Check `sequences[i].protein.id` against the chain you meant, and whether each entity has `unpairedMsaPath`, `unpairedMsa: ""`, or nothing. The single most useful check before a big run — and the only defence against a wrong `--msa-empty-chains` letter on a chain that was going to be empty anyway, which nothing else can catch. |
| `.../af3_msas/`, `.../af3_paired_msas/`, `.../af3_templates/` | Staged copies, named `<chains>_<original basename>` — **a plain copy, not content-addressed**. Empty directory = no chain got one. Two *different* files sharing a basename stage to the same name and are indistinguishable afterwards, so a rerun against a corrected A3M needs a fresh `-l` or a distinct filename. |
| `.../.meta.json` | The flags as given: `no_msa`, `msa_empty_chains`, `msa_file`, `paired_msa_file`, `templates_json`. **Paths, not contents** — it records what you asked for, and cannot tell you whether the file behind a path has since changed. |
| The manifest's 2nd field | `--norun_data_pipeline` or `--run_data_pipeline` — the whole database-free contract. **It is per row, so read every row, not the first**: shards differ, and one `--run_data_pipeline` among them is the signature of a design whose chains you did not fully name. |
| `.../results_shard_<i>/<design>/<design>_data.json` | **What AF3 actually folded**, after the pipeline: the entity list, each entity's `id` and `sequence`, and each entity's *resolved* `unpairedMsa`. This is the only place that distinguishes "the chain was left empty as intended" from "the chain was searched because nobody named it". **None of it is collected**, so the check is per design and by hand. |
| `<leaf>_status` | `OK` is the only proof. Anything else leaves every metric `NA`. |

**Changing MSA, template or chain policy: use a fresh `-l/--dir-label`. This is mandatory, not hygiene.** Shard membership is recomputed each run, so a design can move between shards and leave its old prediction on disk under the old shard. The collector indexes every `results_shard_*` into one dict keyed by design-directory name, so **a stale prediction wins whenever its shard sorts after the new one** — and the sort is lexicographic, i.e. `shard_10` before `shard_2`. Nothing marks the result as stale and `<leaf>_status` still says `OK`.

The staleness is **one-sided, and it is the output half**. The *inputs* are clean (`af3_shards/` is `rmtree`d before each submit, so AF3 never folds a JSON left over from an earlier run, and the manifest's data-pipeline field always describes what executed); the `results_shard_*` directories are not cleared, so an old prediction sits under its old shard index and wins if that index sorts later, with nothing in the table saying so.

**Smoke-test one design first: it is the whole quality gate.** Run it, open its `af3_inputs/<design>.json`, check each entity's `id` and whether it has `unpairedMsaPath` / `unpairedMsa: ""` / nothing, confirm `Skipping data pipeline...` in the `.out`, then scale. Every failure is a container failure, so one design costs one container and a fan-out costs one per shard.

## Invocation

```bash
# de-novo binder against a HOMOTRIMERIC target, database-free — the main use case.
# A,B,C are one entity taking one MSA; D is the binder, left empty.
# Note A:C= — the separator is '=', and the chains are the same mini-language as
# --chains, so 'ABC=' is an ERROR and 'A:C=' is how you name three chains.
# -i mkcomplex_sequence is REQUIRED: a four-segment A/B/C/D sequence only exists after
# mkcomplex. With the default -i this fails with "names chain(s) B,C,D, which this
# complex does not have".
sapia run alphafold3 <run_dir> -t table1 -i mkcomplex_sequence \
    --msa-file A:C=/runs/assets/af3_msa/target_unpaired.a3m \
    --msa-empty-chains D --shard-size 10

# single-chain target after mkcomplex (A = target, B = binder)
sapia run alphafold3 <run_dir> -t table1 -i mkcomplex_sequence \
    --msa-file A=/runs/assets/target.a3m --msa-empty-chains B

# the binder alone, straight off proteinmpnn (one chain, default -i is fine)
sapia run alphafold3 <run_dir> -t table1 --no-msa
```

Resources default to A100-80GB (Modal), 24 CPU, 64 G, 4 h; a 500-residue complex predicts in about a minute.

| Flag | Default | Meaning |
| --- | --- | --- |
| `--shard-size N` | `10` | JSON inputs per task. One shard = one container. |
| `--model-seeds N [N ...]` | `[42]` | AF3 model seeds. |
| `--chains A:D` | all | Narrow which chains to predict. **Permissive where the MSA flags are strict**: a letter the sequence does not have is silently dropped here, an error there. Selecting *nothing* is not an error either — it writes a JSON with zero entities. |
| `--positions 10:50/1:30` | whole chain | Residues to KEEP per chain. **Breaks any `--msa-file` on a sliced chain** (loudly, at run time) and **silently shifts any `--templates-json` indices** (never caught). No submit-time guard for either. |
| `--no-msa` | off | Every chain gets `unpairedMsa: ""`, `pairedMsa: ""`, `templates: []`. Cannot combine with the per-chain flags. |
| `--msa-empty-chains B` | unset | Only the named chains get the empty policy. Chain mini-language (`B`, `A,C`, `B:D`). Expanding to no chains is an error. |
| `--msa-file CHAINS=PATH` | unset | Repeatable. **`=`, not `:`.** `CHAINS` is the ordinary chain mini-language, so `A:C=t.a3m` is chains A, B, C and **`ABC=t.a3m` is an error**. Gives `unpairedMsaPath` to a plain copy under `af3_msas/`, plus `pairedMsa: ""`, `templates: []`. A `:` in the path is now safe. |
| `--paired-msa-file CHAINS=PATH` | unset | Repeatable. Requires `--msa-file` for the same chains; swaps their `pairedMsa` for `pairedMsaPath`. |
| `--templates-json FILE` | unset | Per-chain template entries. Keys use the mini-language (`"A:C"`). Every chain named must already have an MSA policy. |

**The no-flag default is wrong for a database-free setup, deliberately.** With no MSA flag AF3 searches everything and needs the genetic databases on disk. The default is preserved so existing runs do not change behaviour; it is not what you want on Modal unless `SAPIA_MODAL_VOLUME_AF3_DB` is populated.

`--templates-json` takes AF3's own schema:

```json
{"A": [{"mmcifPath": "ref.cif", "queryIndices": [0, 1, 2], "templateIndices": [0, 1, 2]}]}
```

Index lists are **0-based**, equal in length, and must account for residues unresolved in the mmCIF. Max 20 per chain; a multi-chain mmCIF is an AF3 error. A relative `mmcifPath` resolves against the JSON file and is staged with the path rewritten absolute; an entry giving an inline `mmcif` instead is passed through untouched.

**The tool validates none of this.** It copies each `mmcifPath` in and rewrites the key — the schema, the index lengths, the 20-template cap and the chain count are all checked by AF3 when it loads the input, i.e. inside the container. Two consequences worth knowing before you build a template:

- **AF3 requires a release date in the mmCIF**: `_pdbx_audit_revision_history.revision_date`, read by `templates.py::get_polymer_features`. Without it the task dies in featurisation, before inference, with `Template structure <id> must have a release date`. **A CIF written by `utils/boltz_template.py` (gemmi `make_mmcif_document()`) has no such loop and is therefore not AF3-ready** — the revision history has to be carried over from the source RCSB file.
- AF3 takes **one single-chain template per entity**. It does not arrange protomers the way Boltz's forced trimer template does, and **no log line says whether the model actually used the template** — the only way to know is a template-free baseline on the same designs.

## The data-pipeline contract

`run_alphafold.py` passes every `*_database_path` through `replace_db_dir`, which **checks the files exist and raises `FileNotFoundError`** — inside `if _RUN_DATA_PIPELINE.value:`. So with the pipeline on (its default), AF3 aborts at startup when the databases are absent **even if every chain carries its own MSA**.

The manifest records, **per shard**, whether every protein entity in every JSON *in that shard* is fully specified, and the task script passes `--norun_data_pipeline` exactly then, `--run_data_pipeline` otherwise. Consequences:

- **The decision is per shard, and all-or-nothing inside one.** One under-specified design poisons its whole shard (`--shard-size`, default 10) and leaves the rest of the run alone. A mixed table therefore gives a *mixed* outcome, which is the confusing case: on Modal the affected shards die and their designs collect as `missing`, which reads like "those designs failed" rather than "those designs had a chain nobody named".
- **Being specified is a property of the design, not of the flags.** It comes from that design's own chain map, so rows whose input column has a different number of `/` segments fall on different sides of it. That is how a table ends up mixed in the first place — check the `/` count across the column, not just on the first row.
- **A partly-specified run is NOT refused at submit.** Name three chains of a four-chain complex and the run goes out with the pipeline ON. **Count the chains in your flags against the `/` count in the input column before you submit** — nothing else will.
- **On Modal the affected shards fail loudly**: with no `SAPIA_MODAL_VOLUME_AF3_DB`, each dies at startup on `FileNotFoundError`, before any GPU work. Cheap, but it is dead containers rather than a submit-time error, and the surviving shards still collect `OK`, so the run looks partially successful.
- **On vib/SLURM it does not fail at all.** `$AF3_DATABASE` is there, so AF3 happily searches an alignment for the chain you meant to leave empty and the design collects as `OK`. That is the expensive case, and it is now undetectable from the table — only `af3_inputs/<design>.json` shows it.
- The task script passes the manifest's field through verbatim (no whitelist).

## Gotchas, each with its signature

| What you see | What it is |
| --- | --- |
| `chains A,B,C share one sequence and therefore form a SINGLE AlphaFold3 entity` | You named only part of a homo-oligomeric group. Name them together: `--msa-file A:C=<a3m>`. **Do not work around it with `--chains`.** |
| `names chain(s) B,C,D, which this complex does not have` | A multi-chain flag against a single-chain column. Add `-i mkcomplex_sequence`. |
| `--msa-file 'ABC=/p' … chain range 'ABC' must use single letters` | The chains half is the ordinary mini-language now. Write `A:C=/p` or `A,B,C=/p`. |
| `--msa-file 'A:C:/p' must be CHAINS=PATH, e.g. --msa-file A:C=t.a3m` | **The separator is `=`, not `:`**; a recipe written with the old `:` form hits this. |
| `--msa-file: chain(s) X already have an MSA policy` | Two flags claim the same chain — typically `--msa-empty-chains D` plus a range that ran past the target. |
| `FileNotFoundError: ${DB_DIR}/uniref90_...fa not found` at the start of `.err`, on **some** tasks only | Those shards ran **with** the data pipeline and there are no databases. At least one design in each had an entity with no policy; the other shards were fine. Read the manifest row per task, and `af3_inputs/<design>.json` for the designs in it — see the contract above. |
| A record of a reused ColabFold A3M measures one residue longer than the query | **A trailing NUL byte.** A3M files written by a ColabFold MSA server can end in a single `\x00` inside the **last record's sequence** (one record 150 long against a 149-residue query); Boltz consumes such files without complaint. **Strip it before handing the A3M to AF3**, and check the length of every record with lowercase insertion characters removed: the bad record is the last one, exactly what a spot check of the head misses. Whether AF3 would count it as a residue was not tested (a precaution, not an observed failure). |
| `First MSA sequence ... is not the query_sequence` at run time | AF3's own check (`data/msa.py`) and **the only one**: the tool does not read the A3M at submit, so there is no submit-time message for this. Usually means **the chain letters are in a different order than you assumed** (the binder is chain A, not the target): count the `/` in the input column. Fires even with `--norun_data_pipeline`, at the cost of a container per shard. *Unverified against the image's own source; the "a crash, not a wrong number" reading of an unchecked A3M rests on this one claim.* |
| `Template structure <id> must have a release date` in `.err` | The template mmCIF has no `_pdbx_audit_revision_history.revision_date`. Dies in featurisation, before inference. |
| Two designs report identical numbers | Their names sanitise to the same string, so AF3 wrote both into one directory and the collector maps both rows to it. There is no submit-time collision check; names outside `[A-Za-z0-9_-.]` are the risk. |
| A submit that printed nothing and queued nothing | Every error above is raised while building the manifest. **Under `sapia modal-shell --cmd` that still exits 0** — check `<script>_modal.json` exists, never the exit code. |
| `Cannot mask with non-boolean array containing NA / NaN values` from your `-f` | A NA `iptm` on an `OK` row (single-entity input). Use the NA-safe filter below. |
| `Submitting 0 designs` / `No designs to submit.` | Every row is already `OK` (use `--force`), or the input column is empty on every row. |
| First Modal run hangs a long time | **(Modal only)** Building AF3 from source: a C++20 CMake build fetching five git repos plus `components.cif`. Later runs are seconds, because the dependency sync is its own cached layer. |
| `ERROR: expected alphafold3 3.0.4 but the build produced 3.0.2` | **(Modal only)** The clone lost its tag and `setuptools_scm` fell back. A deliberate hard failure — a silent version swap would otherwise be invisible. |
## What it collects

Leaf-prefixed (`alphafold3_*`, or `alphafold3_<dir-label>_*`), from `<design>_summary_confidences.json`:

| Column | Meaning | How to read it |
| --- | --- | --- |
| `<leaf>_status` | **Trust column.** `OK`, `missing: ...`, `error: ...` | **The only proof of success.** Gate every filter on it. Non-`OK` rows are `NA` in every metric, never `0`. Note collect takes **no filter**: it stamps every row with a non-empty input column, so `missing` means "no output dir", not "not submitted". Rows with an empty input column get no status at all. |
| `<leaf>_path` | The best-ranked model, `<name>_model.cif` | AF3 picks the top model itself. Feed it to usalign / epitope / dockq. |
| `<leaf>_ranking_score` | AF3's overall ranking score | Higher is better. An **uncalibrated composite** — nobody here has characterised it. **Rank with it; do not threshold it** unless the bar came from this run's own distribution. Not comparable across targets. |
| `<leaf>_ptm` | Predicted TM-score, whole complex | Fold confidence. Says nothing about the interface. |
| `<leaf>_iptm` | Interface pTM, whole complex | **Diluted by chain count exactly as boltz's is** — on a multi-chain target the target–target pairs drag it toward the deposited assembly. *An* interface number, not *the binder's*. **NA on an `OK` row means a single-entity input** (no pair to score), not a bad prediction. |
| `<leaf>_fraction_disordered` | Fraction AF3 calls disordered | High on a de-novo binder means the binder chain did not fold. |
| `<leaf>_has_clash` | Clash flag — **a float, not a bool** | AF3 writes `0.0`/`1.0`, and it lands in the table as a float. Compare numerically (`> 0`), **never** against the string `"false"`. Non-zero is a hard reject; zero is AF3's coarse flag, not `epitope`'s steric bracket. |

Values are AF3's own, rounded to 2 dp in that JSON, and come from the top-ranked sample.

On disk, `<run_dir>/<table>/<leaf>/results_shard_<i>/<design>/` holds `<design>_model.cif`, `_summary_confidences.json`, `_confidences.json`, `_data.json`, `_ranking_scores.csv`, plus a `seed-<n>_sample-<k>/` directory per sample. **`chain_iptm`, `chain_ptm`, `chain_pair_iptm`, `chain_pair_pae_min`, per-residue pLDDT and the PAE matrix live in those raw files and are NOT collected.**

**Output chain labels:** AF3 splits a multi-letter entity into separate chains with separate `label_entity_id`s, so an `ABC` entity comes back as `label_asym_id` A, B, C (`auth_asym_id` identical). Only the staged `_data.json` keeps the grouped id. Downstream tools naming chains get plain letters.

AF3 names the design directory with `Input.sanitised_name()`: spaces become `_` and characters outside `[A-Za-z0-9_-.]` are dropped. **The tool's copy of this does not lowercase, and real AF3 v3 probably does** (`folding_input.py` builds from `string.ascii_lowercase` on a `.lower()`-ed name). If that is right, any design name containing an uppercase letter collects `missing: no output dir` despite a perfectly good prediction. **Treat as unverified: if a whole run comes back `missing` and your names carry uppercase, this is why.** The collector joins on both the sanitised and raw name. **The submit does not refuse when two names would collide**: one design's output then overwrites the other's and both rows collect the same numbers, both `OK`.

## Filters

**Write the mask NA-safe** — an `OK` row can have a NA `iptm`, and a naive comparison raises.

```python
import pandas as pd


def apply_filter(df):
    ok = df["alphafold3_status"] == "OK"
    iptm = pd.to_numeric(df["alphafold3_iptm"], errors="coerce")
    # has_clash is a FLOAT (0.0/1.0), not a bool and not "True"/"False".
    clash = pd.to_numeric(df["alphafold3_has_clash"], errors="coerce")
    return df[(ok & (iptm > 0.60) & (clash.fillna(1) <= 0)).fillna(False)]
```

**Why an `OK` row can have a NA `iptm`:** AF3 computes `iptm` only across *pairs* of entities, so a single-entity input (a monomer, or a homo-oligomer alone) gives `null`, which the collector turns into `pd.NA` — not `0`, not `None`. The prediction succeeded; there is no interface to score, and `ptm` is still there. Seeing this on a run you expected to be a complex means the chain letters or `-i` are wrong.

**The numeric bars above are placeholders.** `has_clash` is the only column here with a defensible absolute threshold. `iptm` is diluted by chain count, so a bar that is strict on a two-chain complex is loose on a four-chain one. Take any bar from this run's own distribution: it is a property of the run.

The real use is a **second opinion against boltz on the same table** — keep what both like, look hard at what they disagree about:

```python
import pandas as pd


def apply_filter(df):
    af3 = pd.to_numeric(df["alphafold3_iptm"], errors="coerce")
    blz = pd.to_numeric(df["boltz_iptm"], errors="coerce")
    ok = (df["alphafold3_status"] == "OK") & (df["boltz_status"] == "OK")
    return df[(ok & (af3 > 0.6) & (blz > 0.6)).fillna(False)]
```

Both are diluted the same way, so they are comparable in kind — but only within one chain layout.

## What it is blind to

- **Whether the binder binds.** Every column is the predictor's belief; two predictors agreeing is two beliefs, not evidence. Pose belongs to `epitope` (did it land on the epitope, does it clash) and `dockq` (did it keep the reference binding mode); self-consistency to `usalign` against the parent.
- **The binder's own confidence.** `iptm` and `ptm` are whole-complex. Nothing here is restricted to the binder chain or to the residue pairs that touch — `boltzconf` does both for boltz. For AF3, **`af3conf`** provides the binder restriction (binder pLDDT, the binder's `chain_ptm`, binder–target `chain_pair_iptm`) and ipSAE, which is PAE-restricted, **not** contact-restricted; it reads the raw files beside `alphafold3_path`.
- **Energetics and packing.** `pyrosetta` (`if_dG`, `if_dSASA`, `packstat`) and `cms` (CMS, SC).
- **Whether the MSA was the right one.** The tool does not read the A3M at all — it copies the file and points the chain at it. Depth, quality, organism, and even whether the first record is this chain's sequence are yours to check; AF3 catches only the last of those, in the container.
