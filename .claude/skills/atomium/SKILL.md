---
name: atomium
description: >-
  AtomiUM, a private ProteinMPNN-like sequence designer, on Modal: designs sequences onto a backbone table and mints one child row per draw. Load before composing an atomium run, setting any amino-acid bias, or feeding an atomium sequence to mkcomplex, epitope or a predictor. Traps: `atomium_sequence` IS the full designed chain (pinned residues included -- an older version of this skill said otherwise and was wrong); read `atomium_len_matches_parsed` and `atomium_fixed_held` before any sequence; `--bias-aa` is global and far stronger than exp(bias); and tables collected before 2026-10-14 carry one bogus `_a0` row per parent holding the NATIVE input sequence at status OK.
---

# atomium

**Custom tool** (lives in `tools/atomium/`, not bundled with prosapia).
**`action: create`** — mints a child table, one row per designed sequence. `-t` required
in practice (it designs onto an existing backbone table).

AtomiUM is a private PyG re-implementation of ProteinMPNN-style sequence design, with
**noise-conditioned weights**: one checkpoint per training noise level. It ships
ProteinMPNN's own helper scripts and takes the same jsonl inputs, so this tool is
deliberately the `proteinmpnn` tool with the same two mini-languages, the same
signature-grouping/bin-packing, and the same lineage contract. **Everything you know
about composing a proteinmpnn run transfers** — read that skill for the mini-languages.

Repo: `AndreiSokolovskii/develop_atomium`, branch `pure_wo_jit`, **private**.

## Premise, and when its numbers mean nothing

**One sequence per draw, for the chains you designed.** The FASTA AtomiUM writes
contains the **designed chains only** — a chain that is not designed never appears in
it. Every column here is a property of that string and of the input chain it was
designed onto.

**Above one designed chain, the `/` join and the order of the chains within it are
ASSUMED, not measured.** The probe below had a single designed chain; the convention is
read off ProteinMPNN's printing code. The collector is therefore written not to depend
on it: `len_matches_parsed` and `init_matches_parsed` compare **multisets**, and
`init_order_matches` reports the element-wise answer **without ever failing a row**. The
first real multi-chain run settles the question from the table — if it comes back `1`,
the order is `assigned_pdbs.jsonl`'s; if `0` while `init_matches_parsed` is `1`, it is
not.

**`atomium_sequence` is the FULL designed chain, pinned residues included.** Measured
2026-10-14 on real output (run_dir `20261009_130629_5M2J_vhh_redmotif_lowtemp_scaleup`,
child table `table1_atomium_probe`, Volume `sapia-runs-toon`, `--chains-to-design D`
with 17 pinned positions), from the raw files:

| design | chain-D length | generated length | `_init_sequence` |
| --- | --- | --- | --- |
| `scaleup0_diff_1` | 144 | **144** (both draws) | byte-identical to native chain D |
| `scaleup1_diff_0` | 120 | **120** (both draws) | byte-identical to native chain D |

The 17 pinned characters read `TSNYWYENTNLIKRSPF` in the native chain **and in all four
generated sequences**. `_init_sequence` was byte-**un**equal to native-minus-the-pinned
(which would have been 127 and 103). No `/` or `:` appeared in any entry.

> **An earlier version of this skill claimed the opposite** — that under `--fixed-positions`
> the FASTA omits the fixed residues, and that a column `atomium_sequence_full` carried
> the complete chain. **That was false, and that column never existed.** Do not reinsert
> anything, and do not look for `sequence_full`, `init_matches_parsed` in its old sense,
> `fixed_residues` as a reinsertion audit, or `len_full_matches_parsed`.

Where the numbers stop meaning anything:

- **More than one designed chain.** `atomium_fixed_held` is **NA**, not 0: the
  position-to-chain mapping of the `/`-joined string is not established here, and a guess
  is exactly the error that column exists to catch. `len_matches_parsed` still holds (it
  compares the multiset of per-chain lengths, so it is order-independent). Read
  `atomium_init_order_matches` on such a run and report what it says.
- **A design with no `parsed_pdbs.jsonl` in its `grp_` dir** (a hand-edited or
  half-deleted out_dir): every trust column is NA and the row is an error. Nothing is
  inferred from the FASTA alone.
- **`atomium_seq_rec` across runs.** It is recovery against *the input sequence*, over
  the designable positions. Two runs that pinned different numbers of positions are not
  comparable on it.

## Trust columns — read these before any sequence

| Column | Expected | How to read it |
| --- | --- | --- |
| `atomium_status` | `OK` | **The only proof of success.** Anything starting `error:` means a check below failed; the row still carries the raw FASTA fields but must not be used as a sequence. |
| `atomium_init_matches_parsed` | `1` | AtomiUM's own `_init_sequence` echo equals the parsed input chain(s) this row claims. `1` is the proof that the row's frame is the chain you think it is. `0` only ever appears on an error row; **NA** means the check could not be made (no echo, no parsed entry). |
| `atomium_len_matches_parsed` | `1` | The draw's per-chain lengths equal the parsed designed chains' lengths. `0` ⇒ truncation or a chain-composition surprise, and the row is an error. **NA** on a row whose input could not be read. |
| `atomium_fixed_held` | `1` when positions were pinned | Every pinned position carries the input residue in the generated sequence. **NA** when `n_fixed == 0` or when more than one chain was designed — never a vacuous `1`. |
| `atomium_n_fixed` | what you meant to pin | How many positions this design **actually** pinned, read from the task's own `fixed_pdbs.jsonl`. Compare it against your intent on the first row of every run. **NA**, never `0`, when that file exists but has no entry for this design or for a designed chain — a gap there is an error row, not a confident zero. |
| `atomium_input_pdb` | the parent's structure path | **Which structure this sequence was designed on**, read from the symlink the submitter staged. Cross-check it against the parent row's `rfdiffusion3_path`. Caveat: the trust columns are computed against *the same submit's* `parsed_pdbs.jsonl`, so a FASTA left behind by an earlier failed submit is self-consistent and reads as trusted — this column does not close that race, it makes it visible. |
| `atomium_init_order_matches` | `1` (single chain) | **Informational, never fails a row.** Element-wise version of `init_matches_parsed`. On a multi-chain run, `0` here with `init_matches_parsed == 1` means the FASTA's chain order is not `assigned_pdbs.jsonl`'s — the open question above. |

A row where all of these read `1` (or NA for `fixed_held` with nothing pinned) is a row
whose sequence is the designed chain, at the input's length, with the pinned residues
intact. That is what used to require a GPU probe.

## Verified invocation

```bash
sapia run atomium <run_dir> -t table0 --num-seq-per-target 2
sapia collect atomium <run_dir> -t <the table the run reserved>
```

`default_input_column` is **`rfdiffusion3_path`** — unlike the bundled proteinmpnn,
whose default is the older `rfdiffusion_path` and silently submits nothing after an
rfd3 run. Coming from a different parent, pass `-i` explicitly.

## Flags

| Flag | Default | Note |
| --- | --- | --- |
| `--model-noise` | `n05` | **The one flag that is not in proteinmpnn.** Which weight file to design with: `n00`…`n07` (training sigma 0.0 Å…0.7 Å), plus the `n05_v2` retrain. Validated at submit time, so a typo fails fast instead of per container. |
| `--num-seq-per-target` | `2` | Sequences per backbone **per temperature**. |
| `--sampling-temp` | `0.1` | One or more space-separated temperatures: `'0.1 0.2'`. See the count rule below. |
| `--seed` | `37` | **`0` means pick a random seed**, it does not mean seed zero. |
| `--designs-per-task` | `10` | Designs bin-packed per task; identical-param designs share one batched call. |
| `--chains-to-design` | `''` (all chains) | Chain mini-language (`A:C,E`). **A letter this design's own structure does not have now refuses the submit** — it used to design nothing and still collect `OK`. |
| `--chains-to-design-column` | `''` | Per-design chain list from a table column. **Mutually exclusive with `--chains-to-design`** (it used to win silently). A missing/empty value refuses the submit. |
| `--fixed-positions` | `''` | Position mini-language, 1-indexed **within each parsed chain**. Groups (`/`) map one-to-one onto the designed chains. Requires a chain list. |
| `--fixed-positions-column` | `''` | Per-design spec from a table column. **Mutually exclusive with `--fixed-positions`.** Missing/empty refuses the submit. |
| `--tied-positions` | `''` | Same grammar; mutually exclusive with `--symmetry`. |
| `--symmetry` | `None` | `auto\|N`. Ties all chains and lets one position group describe a single asymmetric unit (broadcast across the designed chains). |
| `--bias-aa` | `''` | Global `AA:bias` pairs. Read the calibration section — it is not `exp(bias)`. |
| `--set` | `[]` | Raw flags forwarded to `atomium.py` (repeatable). |

Default Modal resources **(Modal only)**: L4, 8 CPU, 16 GiB, 1 h timeout. Weights ship
inside the image — no cache Volume, no download step. The builder does **not** force
`gpus_per_task = 0`; this is a GPU tool.

**Submit-time refusals** (all raise before anything is queued or staged, and
per-design failures are **collected into one raise naming every offender**, so a wide
table is repaired in a single edit): a position outside `1..chain_length`; a number of
`/`-groups that does not equal the number of designed chains; both flag and column of a
pair; a column absent from the table; a per-row value that is empty/NA; a chain letter
absent from that design's PDB; and **a chain spec that is non-empty but expands to no
chains** (`[]`, `,`, `[ ]`). That last one is the nastiest of the set: an empty
`--chain_list` means *design every chain* to the helper scripts, so it would have
designed the whole structure, pinned **nothing**, and collected `status OK` with
`n_fixed == 0` and `fixed_held` NA — a motif-graft design whose motif was never held,
indistinguishable in the table from a legitimate no-pin run.

### Picking `--model-noise`

The suffix is the coordinate noise the checkpoint was **trained to denoise**, so it is a
statement about how much you trust the input geometry:

- **Low (`n00`–`n02`)** — trusts the backbone as given. For crystal structures or
  already-relaxed models.
- **Mid (`n03`–`n05`)** — the useful default range for **de novo backbones**. `n03` is
  AtomiUM's own default; this tool defaults to `n05`.
- **High (`n06`–`n07`)** — tolerates rough or coarsely sampled backbones.

Noise level is a *cheap* axis to scan: run the same parent table twice with different
`--table-label`s and compare downstream self-consistency, not the sequences themselves.

```bash
sapia run atomium <run_dir> -t table0 --table-label n03 --model-noise n03
sapia run atomium <run_dir> -t table0 --table-label n05 --model-noise n05
```

## The sampling count is a product

`BATCH_COPIES = num_seq_per_target * len(temperatures)` — each temperature gets a **full**
`--num-seq-per-target` set. So `--num-seq-per-target 4 --sampling-temp '0.1 0.2 0.3'` is
**12 sequences per backbone**, not 4. Easy way to accidentally triple a run.

`atomium.py`'s own `--batch_size` affects **chunking only** and does not change the
sampling count; this tool exposes no `--batch-size`. (A previous note here said it was
"live at the pinned commit `dee5a4d…`" — **`dee5a4d` is an atomiumph commit, not this
tool's pinned checkout**, which is `ATOMIUM_COMMIT = def1d6a…`. Its status at the pinned
commit is unverified. It was dead at `292e29c`.)

## `--bias-aa` is far stronger than `exp(bias)` — calibrate it, never reason about it

`--bias-aa` takes space-separated `AA:bias` pairs and writes a single run-wide `bias_AA.jsonl`
(`make_bias_AA.py` is purely generative, no structure input). Two things follow, and both have
bitten:

**It is GLOBAL, not positional.** One logit offset applied to every position of every design in
the run. There is no per-position bias path in the manifest builder. So it cannot place a residue
*somewhere specific* — it raises that residue everywhere, buried core included. `--fixed-positions`
does not help either: it *freezes* whatever the input structure already has, it does not choose.
**To get a residue at a particular place, oversample and select on a measured column.**

**Its magnitude does not follow `exp(bias)` at low sampling temperature.** Measured on EGFR dIII
binders, 2310 sequences per arm, `--model-noise n05`, default `--sampling-temp 0.1`:

| `--bias-aa` | global His frequency | mean His per ~80-residue binder |
| --- | --- | --- |
| none | 0.615% | 0.49 |
| `H:1.39` | **14.238%** | ~11.4 |

`exp(1.39) = 4.0`, so a plain logit offset would predict ~2.5%. The observed shift is **23×, not
4×** — sequences with up to **52** histidines. The mechanism is almost certainly the temperature:
at T = 0.1 the distribution is nearly deterministic, so a modest shift flips the biased residue to
the argmax at every position where it was anywhere near competitive.

**Practical rule: fit the multiplier from two runs before trusting a third.** Those two points
imply `exp(2.26 × bias)` here, so the bias landing the mean at 1.0 is ≈ `H:0.31`, not 1.39. A
bias tuned at one `--sampling-temp` does not transfer to another.

**Nor does it transfer across amino acids, and the SIGN matters more than the magnitude.**
Measured 2026-10-01 on 1107 EGFR dIII binder sequences, `--model-noise n05`,
`--sampling-temp 0.1`, `--bias-aa C:-1.0`, against a 10,944-sequence unbiased baseline:

| | unbiased | `C:-1.0` | ratio |
| --- | --- | --- | --- |
| Cys frequency | 1.41% | **0.774%** | **1.8×** |
| mean Cys per ~75-residue binder | 1.11 | **0.584** | **1.9×** |

`exp(2.26 × bias)` predicted **10×**. Observed **1.9×** — and even plain `exp(bias)` would predict
2.7×, so **suppression came in weaker than a naive logit offset**, while the His enrichment above
came in 5.7× *stronger* than naive. Fitting this arm alone gives `exp(0.64 × bias)`.

**The asymmetry is the rule: at low temperature, enrichment is amplified and suppression is
attenuated.** So budget a negative bias of roughly `exp(0.6 × bias)` and a positive one of roughly
`exp(2.3 × bias)` as starting points, and expect to need a filter as well when suppressing.

**A bias cannot break a COUPLED choice, only reduce how often it is made.** In the same measured
arm, cysteines still came in pairs: **93.7% of sequences carried an EVEN number** — essentially
unchanged from the unbiased 94%. Where a residue's placement is structurally coupled rather than
sampled per position, **gate on the measured count instead**.

**And a strong bias can make an "exactly one of X" criterion WORSE, not better.** The count of
sequences carrying exactly one copy peaks when the mean is 1.0; pushing past that overshoots
(210 designs / 25 backbones unbiased → 28 / 8 at `H:1.39`).

**Always verify the bias actually applied**: `ls <out_dir>/bias_AA_*.jsonl`, then
`cat` the one whose path appears after `--bias_AA_jsonl` in the staged task file
(`<out_dir>/atomium_tasks/task_*.tsv`). An unbiased arm has no such file.

The bias jsonl is **content-addressed** (`bias_AA_<digest>.jsonl`) as of 2026-10-14, so
two arms in one leaf no longer overwrite each other. **On a leaf written before that**
the file is the fixed name `bias_AA.jsonl` and the *last* submit's bias is the only one
on disk — so for older rows the audit above can report a bias those rows were not
designed with. Check the date of the leaf before trusting it.

## The image is built from a private repo with Modal Secrets *(Modal only)*

`modal_image.py` clones the repo at **build** time using the `github-token` and
`github-username` Modal Secrets, then deletes `.git` so the credentialed remote is not
baked into the image. It has to happen at build time: prosapia's Modal executor hardcodes
the task secret to the run's `.env` and exposes no per-tool `secrets()` hook.

- **You need access to those two Secrets** in the Modal workspace, or the image build fails.
- The checkout is **pinned to a commit** (`ATOMIUM_COMMIT`); bump it deliberately.
- First run pays a **long image build inside the submit call** (torch 2.11 + cu128 and
  ~470 MB of weights). Allow 15+ minutes before assuming a first submit is stuck.
- `torch_cluster` is a **required** dependency (`model_lib.py` calls
  `torch_geometric.nn.radius_graph`). Without it every task dies at model construction.

## Gotchas, each with its signature

| Gotcha | What you will see |
| --- | --- |
| **The `_a0` remedy must be applied to CHILD tables too.** Any `update` tool run on an affected table (boltz, usalign, epitope, cms …) annotated the `_a0` row in place, so those numbers are measurements **of the native input sequence** at `status OK`; any `create` tool minted `<parent>_a0_f*` / `_a0_a*` children from it. | Grep every descendant table for row names containing `_a0`; drop them wherever they appear, not just in the atomium table. A design-count that is 1.5× what you asked for (3 rows where `--num-seq-per-target 2`) is the signature at every generation. |
| **Tables collected before 2026-10-14 contain a bogus row per parent.** The collector matched the echo's header with `==` against `_init_sequence`, but the real header is `>_init_sequence, chains=D`, so the test never fired and the echo was minted as a design. | `--num-seq-per-target 2` produced **3** rows per parent: `<parent>_a0`, `_a1`, `_a2`. The `_a0` row carries the **native input sequence** at `atomium_status == OK`, with `atomium_sample` and `atomium_temperature` **NaN** and `atomium_iteration == 0`. Drop it with `df[df["atomium_iteration"] != 0]`, or re-collect with `--force` and then require `atomium_len_matches_parsed == 1` (the stale `_a0` rows are never re-emitted, so they stay NA there). |
| A design already staged under this leaf in a different subgroup | `ValueError: N design(s) are already staged under this output dir ... in a DIFFERENT subgroup`. It means a second submit to the same table+leaf changed a design's chains/fixed/tied signature or its batch. Use `--force` (replaces only subgroups wholly restaged by this submit), a fresh `-l/--dir-label`, or delete the named `grp_` dir. |
| An extra `/`-group in a position spec | `ValueError: ... has 2 '/'-separated chain group(s) but 1 chain(s) are designed`. Before 2026-10-14 this pinned only the first group and said nothing. |
| A chain letter the structure lacks | `ValueError: ... names chain(s) Z, which are not polymer chains of <pdb> (present: A, D)`. Before 2026-10-14 the run designed nothing and the row still collected `OK` carrying the input sequence back. |
| Both a flag and its column | `ValueError: --fixed-positions ... and --fixed-positions-column ... are mutually exclusive`. Before 2026-10-14 the column won silently. |
| `seq_rec` exactly `0.0` | AtomiUM's `sequence_recovery` returns `0.0`, not NA, when nothing is designable. Read it beside `atomium_n_fixed` and `atomium_parsed_len`. |
| Row names | `<parent>_a1`, `_a2`, … — the index is over **all** FASTA entries, echo included, so draws start at `_a1`. Deliberate: row names are lineage keys and renumbering would re-key every existing table. |

## What it collects

Child rows named `<parent>_a1`, `<parent>_a2`, … each linked to its parent backbone.
(The `_a` suffix distinguishes them from proteinmpnn's `_f` rows.) All columns are
leaf-prefixed `atomium_` (or `atomium_<dir_label>_`).

**Trust columns first** — `status`, `init_matches_parsed`, `len_matches_parsed`,
`fixed_held`, `n_fixed`: see the table near the top.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `status` | **TRUST.** `OK`, or `error: …` | The only proof of success. |
| `init_matches_parsed` | **TRUST.** `1`/`0`/NA | `_init_sequence` equals the parsed designed chain(s). |
| `len_matches_parsed` | **TRUST.** `1`/`0`/NA | The draw's chain lengths equal the parsed chains'. |
| `fixed_held` | **TRUST.** `1`/`0`/NA | Pinned positions carry the input residue. NA when nothing was pinned or >1 chain designed. |
| `n_fixed` | **TRUST.** int / NA | Positions actually pinned, from the task's `fixed_pdbs.jsonl`. `0` is a measurement (the file is absent ⇒ nothing was pinned); **NA** means the file exists but could not account for this design or chain — an error row. |
| `input_pdb` | **TRUST.** The staged structure's path / NA | Which backbone this sequence was designed on. |
| `init_order_matches` | `1`/`0`/NA | Element-wise `init_matches_parsed`. Informational only; never changes a status. |
| `sequence` | The designed chain(s), `/`-joined | **This is the binder.** Feed it to mkcomplex / epitope / a predictor. Not Boltz's default input column — pass `-i atomium_sequence`, or `-i mkcomplex_sequence` on a binder table. |
| `seq_len` | Residues in `sequence`, separators excluded | Equals `parsed_len` on a healthy row. |
| `parsed_len` | Residues in the parsed designed chain(s) | The input length this draw was designed onto. |
| `n_chains_designed` | int / NA | How many chains the FASTA covers. `1` for an ordinary binder redesign. |
| `fixed_positions` | Audit string, `D:3,11` (`;` between chains) | Chain and 1-based positions **in the order the helper listed them**, not sorted. Compare against the spec you launched with. **NA** (not `""`) when nothing was pinned. |
| `fixed_residues` | **Chain-prefixed**: `D:HH` | The residues at the pinned positions in the input chain, same order. **Should be what you meant to pin** (all `H` for a His-pinning run). NA when nothing was pinned or >1 chain designed. The prefix is not decoration: a bare two-residue pin of Asn+Ala spells `NA`, which is in `read_csv`'s default `na_values` and would come back from the TSV as a null. |
| `sample` | AtomiUM's sample id | **Not unique per row**: `n % num_seq_per_target`, so with several temperatures the same id recurs. Identify a draw by (`sample`, `temperature`) or by the row name. |
| `temperature` | The sampling temperature of this draw | |
| `seq_rec` | Recovery against the INPUT sequence | **Not a quality score.** Measured 0.21–0.54 on rfd3 backbones, varying systematically **by backbone** — it appears to track backbone designability. Do not rank by it; do not dismiss it. |
| `iteration` | Index of this entry in the FASTA | `1`-based for draws (the `0` entry is the skipped echo). A value of `0` on a row means it was collected by the pre-2026-10-14 collector. |
| `path` | The FASTA | One file per design, shared by all its draws. |

On disk: `<out_dir>/grp_<digest>/seqs/<parent>.fa`, beside the task's own
`parsed_pdbs.jsonl`, `assigned_pdbs.jsonl` and `fixed_pdbs.jsonl` — which is where every
trust column comes from. `<out_dir>` is `run_dir/<child table>/<leaf>`, so a different
`--table-label` or `-l/--dir-label` is a different directory tree entirely.

Subgroup dirs and task files are **content-addressed over the designs and their
chains/fixed/tied signature**, so two submits collide only when they stage *the same
designs with the same signature* — and then sharing the directory is correct, not a
collision (`test_identical_resubmit_is_idempotent` asserts exactly that). **The cost of
that design:** `--model-noise`, `--sampling-temp`, `--seed` and `--bias-aa` are **not**
in the digest. Re-submitting the same designs into the same child table and leaf at a
different sampling arm therefore reuses the directory and **overwrites their FASTAs**,
and the next collect rewrites `<parent>_a1…` in place while any child table still points
at those row names. **Give every arm its own `--table-label` _and_ `-l/--dir-label`** —
the noise-scan example above is safe only because its `--table-label`s differ.

### There is no score column — this changes how you rank

AtomiUM's FASTA reports only sequence recovery. **There is no `score` or `global_score`**,
so the model-likelihood ranking that `proteinmpnn_score` supports **does not exist here**.
Rank atomium designs by what comes **after** them — predictor confidence and, above all,
self-consistency (`usalign` back to the parent backbone).

## A filter you can actually write

```python
# keep only draws whose sequence is provably the designed chain, with its pins intact
def apply_filter(df):
    ok = (df["atomium_status"] == "OK") & (df["atomium_len_matches_parsed"] == 1)
    held = df["atomium_fixed_held"].fillna(1) == 1          # NA = nothing was pinned
    return df[ok & held & (df["atomium_n_fixed"] == 17)]    # 17 = this campaign's pin count
```

The pin count is a property of the campaign, not of the table — record it in the campaign
log. On an older table, add `& (df["atomium_iteration"] != 0)` to drop the bogus `_a0`
echo rows.

## What it is blind to

- **Whether the sequence folds to the backbone.** Nothing here is a structure. That is
  `boltz`/`alphafold3` followed by `usalign` self-consistency against the parent.
- **Whether the binder still binds.** Use `epitope` (`--sequence-column atomium_sequence`),
  `cms`, `pyrosetta` interface metrics.
- **Sequence quality.** `seq_rec` is similarity to the input, not a likelihood; there is
  no score column to rank on.
- **Whether the positions you pinned were the right ones.** `fixed_held` proves they were
  *held*; whether they were the correct residues in the correct numbering is the upstream
  tool's guarantee (`motifrmsd` makes it explicitly).
- **Redundancy across draws.** Many draws from one backbone are near-identical; use
  `seqcluster`.
