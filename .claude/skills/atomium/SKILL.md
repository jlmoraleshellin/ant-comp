---
name: atomium
description: >-
  AtomiUM, a private ProteinMPNN-like sequence designer, on Modal. Load before composing an atomium run, setting any amino-acid bias, or feeding an atomium sequence to mkcomplex or a predictor. Traps: `--bias-aa` is global and far stronger than exp(bias), and with `--fixed-positions` only `atomium_sequence_full` carries the pinned residues.
---

# atomium

**Custom tool** (lives in `tools/atomium/`, not bundled with prosapia).
**`action: create`** — mints a child table, one row per designed sequence.

AtomiUM is a private PyG re-implementation of ProteinMPNN-style sequence design, with
**noise-conditioned weights**: one checkpoint per training noise level. It ships
ProteinMPNN's own helper scripts and takes the same jsonl inputs, so this tool is
deliberately the `proteinmpnn` tool with the same two mini-languages, the same
signature-grouping/bin-packing, and the same lineage contract. **Everything you know
about composing a proteinmpnn run transfers** — read that skill for the mini-languages.

Repo: `AndreiSokolovskii/develop_atomium`, branch `pure_wo_jit`, **private**.

## Verified invocation

```bash
sapia run atomium <run_dir> -t table0 --num-seq-per-target 2
sapia collect atomium <run_dir> -t <the table the run reserved>
```

`default_input_column` is **`rfdiffusion3_path`** — unlike the bundled proteinmpnn,
whose default is the older `rfdiffusion_path` and silently submits nothing after an
rfd3 run. Coming from a different parent, pass `-i` explicitly.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--model-noise` | `n05` | **The one flag that is not in proteinmpnn.** Which weight file to design with: `n00`…`n07` (training sigma 0.0 Å…0.7 Å), plus the `n05_v2` retrain. Validated at submit time, so a typo fails fast instead of per container. |
| `--num-seq-per-target` | 2 | Sequences per backbone **per temperature**. |
| `--sampling-temp` | `0.1` | One or more space-separated temperatures: `'0.1 0.2'`. See the count rule below. |
| `--seed` | 37 | **`0` means pick a random seed**, it does not mean seed zero. |
| `--designs-per-task` | 10 | Designs bin-packed per task; identical-param designs share one batched call. |
| `--chains-to-design`, `--fixed-positions`, `--tied-positions`, `--symmetry`, `--bias-aa`, `--set` | | **The flag GRAMMAR is proteinmpnn's**, including the guardrail that positions require a chain list. **The OUTPUT is not.** With `--fixed-positions` — and for any chain that is not designed — AtomiUM's FASTA **does not contain the fixed residues at all**, so `atomium_sequence` is not the binder. Take **`atomium_sequence_full`**; see the trap section below. |

Default Modal resources: **L4**, 8 CPU, 16 GiB, 1 h timeout. Weights ship inside the
image — no cache Volume, no download step.

### Picking `--model-noise`

The suffix is the coordinate noise the checkpoint was **trained to denoise**, so it is a
statement about how much you trust the input geometry:

- **Low (`n00`–`n02`)** — trusts the backbone as given. For crystal structures or
  already-relaxed models.
- **Mid (`n03`–`n05`)** — the useful default range for **de novo backbones**, which carry
  generator-specific geometry quirks. `n03` is AtomiUM's own default; this tool defaults
  to `n05`.
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

There is **no `--batch-size`**: `atomium.py` accepts the flag but never reads it. Do not
pass it through `--set` expecting an effect.

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
attenuated.** Near-deterministic sampling means a positive bias flips the residue to argmax
everywhere it was merely competitive, while a negative bias only bites where its margin was
already thin. Suppression saturates against the model's conviction. **So budget a negative bias
of roughly `exp(0.6 × bias)` and a positive one of roughly `exp(2.3 × bias)` as starting points,
and expect to need a filter as well when suppressing.**

**A bias cannot break a COUPLED choice, only reduce how often it is made.** In the same measured
arm, cysteines still came in pairs: **93.7% of sequences carried an EVEN number** (818 with zero,
only 36 with one, 193 with two, 28 with three, 25 with four) — essentially unchanged from the
unbiased 94%. Two cysteines stayed five times commoner than one. The bias talked some sequences
out of building a disulfide at all; it never broke a pair apart. Where a residue's placement is
structurally coupled rather than sampled per position, **a per-position logit offset is the wrong
instrument — gate on the measured count instead, and treat the bias as a way to improve the yield
of that gate, not to replace it.**

**And a strong bias can make a "exactly one of X" criterion WORSE, not better.** The count of
sequences carrying exactly one copy peaks when the mean is 1.0. Pushing the mean past that
overshoots: the same selection went from 210 designs / 25 backbones (unbiased, mean 0.49) to
28 / 8 at `H:1.39`. Loose criteria (`>= 1`) improve; strict ones collapse.

**Always verify the bias actually applied**: `cat <out_dir>/bias_AA.jsonl` and confirm
`--bias_AA_jsonl` appears in the staged task argv. An unbiased arm has no such file. Without that
check a no-op bias produces a run that looks entirely normal and silently duplicates its control.

## The image is built from a private repo with Modal Secrets

`modal_image.py` clones the repo at **build** time using the `github-token` and
`github-username` Modal Secrets, then deletes `.git` so the credentialed remote is not
baked into the image.

**This has to happen at build time.** prosapia's Modal executor hardcodes the task
secret to the run's `.env` (`executors/modal.py`) and exposes no per-tool `secrets()`
hook, so a named Secret is simply not reachable from a running task. Cloning in the
task (as an ad-hoc Modal script would) also re-clones a private repo on every container.

Consequences:

- **You need access to those two Secrets** in the Modal workspace, or the image build
  fails. They are workspace-level, not personal.
- The checkout is **pinned to a commit** in `modal_image.py` (`ATOMIUM_COMMIT`). The
  branch is under active development, so bump that constant deliberately — a rebuild
  will not drift on its own.
- First run pays a **long image build inside the submit call** (torch 2.11 + cu128 and
  ~470 MB of weights). Allow 15+ minutes before assuming a first submit is stuck; later
  runs are seconds.

`torch_cluster` is a **required** dependency and easy to drop: nothing imports it
directly, but `model_lib.py` calls `torch_geometric.nn.radius_graph`, which is a
torch-cluster binding. Without it every task dies at model construction.

## What it collects

Child rows named `<parent>_a1`, `<parent>_a2`, … each linked to its parent backbone.
(The `_a` suffix distinguishes them from proteinmpnn's `_f` rows.)

Columns (leaf-prefixed `atomium_`): `sequence`, `sample`, `temperature`, `seq_rec`,
`iteration`, `sequence_full`, `n_fixed`, `fixed_positions`, `fixed_residues`,
`init_matches_parsed`, `len_full_matches_parsed`, plus `_status` and `_path` (the FASTA).

`atomium_sequence` is what the structure predictors consume **when no position was
fixed**. It is **not** Boltz's default input column (that is `proteinmpnn_sequence`), so
a Boltz run after atomium needs `-i atomium_sequence` — or `-i mkcomplex_sequence` on a
binder table, where mkcomplex still belongs in between. **With `--fixed-positions`, read
the next section before using `atomium_sequence` for anything.**

### With `--fixed-positions`, `atomium_sequence` is NOT the binder

**MEASURED** on run `table2_merged_r2cap12_his_fix` (6020 children of 602 parents, label
`his_fix`, `--chains-to-design D --fixed-positions '{hisbond_af3_centre1_idx1},{hisbond_af3_centre2_idx1}'`): atomium.py writes **only the designable positions** into the FASTA. So every child's `sequence` is shorter than its input chain **by the number of fixed positions** (−2 on all 6020 rows) and **does not contain the fixed residues at all** — the very histidines the run was pinning. `_init_sequence` is the input chain with the same positions deleted, so it does not contradict this. For 602 of 602 parents the parent carried H at both fixed positions and the echo equalled the parent chain with exactly those two deleted.

Two consequences:

- **`atomium_seq_rec` is computed over the designable positions only**, so it is not comparable with a `seq_rec` from a run without fixed positions.
- **Downstream steps must take `atomium_sequence_full`, never `atomium_sequence`.** Handing `sequence` to mkcomplex and AF3 folds a binder two residues short with no pinned histidine — a plausible-looking prediction of a molecule nobody designed.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `init_matches_parsed` | **TRUST, read first.** `1` when AtomiUM's own `_init_sequence` echo equals the parsed input chain with exactly the fixed positions removed | `1` is the proof that the positions reinserted are the ones AtomiUM actually held fixed. It is never `0` on an `OK` row — a mismatch makes the row an **error row** — so what you are checking is that the row is `OK`. Vacuously `1` when nothing was fixed. |
| `len_full_matches_parsed` | **TRUST.** `1` when `len(sequence_full)` equals the parsed input chain's length | Same contract: `0` only ever appears on an error row. Vacuously `1` when nothing was fixed. |
| `sequence_full` | The designed sequence **with the fixed residues reinserted** — the whole binder chain | **This is the column every downstream step wants.** Equals `sequence` exactly when `n_fixed == 0`. **Empty on an error row**, never a shortened fallback. |
| `n_fixed` | How many positions were held fixed in the designed chain | `2` on the `his_fix` run, `0` on every older table. NA on an error row whose fixed list could not be read. |
| `fixed_positions` | Audit string, `D:49,12` — chain and 1-based positions **in the order given**, not sorted | Compare against the `{hisbond_af3_centre<k>_idx1}` values the run was launched with. The order is the islands' order; the reinsertion does not depend on it. |
| `fixed_residues` | The residues put back, in that same order (e.g. `HH`) | **Should be all `H` for a His-pinning run.** Anything else means the position pinned was not a histidine in the input chain — the run pinned the wrong residue, and the error is upstream in the position columns. |

**The reinsertion, and what makes it safe.** The collector reads `grp_<g>/parsed_pdbs.jsonl` (one line per design, `seq_chain_<X>` = the full input chain), `grp_<g>/fixed_pdbs.jsonl` (one line per group, `{design: {chain: [1-based positions]}}`, absent when nothing was fixed) and `grp_<g>/assigned_pdbs.jsonl` (one line per group, `{design: [[designed chains], [fixed chains]]}`). It takes each fixed position's residue from the parsed chain and fills **every other position, in order**, from the designed string — so a descending fixed list like `[49, 12]` cannot scramble the result. Files are read once per group and cached (163 groups, parsed lines up to ~54 KB).

**Error behaviour — a doubtful row never carries a sequence.** The row collects with a status starting `error:` and an **empty** `sequence_full` when: the echo does not match (`error: _init_sequence does not equal chain D with positions 49,12 removed …`); the designed string is not `len(chain) − n_fixed` (`error: designed sequence is 55 residues but the parsed chain has 58 with 2 fixed …`); a fixed position is outside the chain or listed twice; the parsed entry or its `seq_chain_<X>` is missing; or **more than one designed chain** (`error: 2 designed chains (D,A); the FASTA joins them with '/' …` — multi-chain reinsertion is deliberately **not implemented** rather than guessed). The FASTA's own columns (`sequence`, `sample`, `temperature`, `seq_rec`) are still filled on such a row — they are exactly what AtomiUM wrote; only `sequence_full` is withheld.

**Older tables are unaffected.** A design with no fixed positions gets `sequence_full == sequence`, `n_fixed == 0`, both trust columns `1`, and **every pre-existing column keeps its name, value and meaning** — re-collecting an old leaf changes nothing but adds the six columns.

**Known upstream, not fixed here.** Three residual risks live in prosapia's position parser and in AtomiUM itself; this collector does not paper over any of them:

- **`_parse_positions` does not check the group count against the designed chains.** `--chains-to-design D --fixed-positions '49/12'` gives two groups for one chain and **silently pins only 49** — no error anywhere. **Compare `atomium_n_fixed` against what you meant to pin** on the first row of every run; it is the only place the discrepancy shows.
- **`sequence_recovery` returns `0.0`, not NA, when nothing is designable.** A `seq_rec` of exactly `0.0` can therefore mean "no designable position" rather than "recovered nothing" — read it beside `n_fixed` and the chain length.
- **The multi-chain `/` split mis-assigns residues when fixed positions precede the last designed chain.** This collector **refuses** that case (the error row quoted above) instead of producing a plausible wrong `sequence_full`; the underlying parser does not.

**What has NOT been checked** (the reinsertion was exercised only on a synthetic output dir built from the real values of one sampled design, plus four deliberately broken variants):

- **A fixed-positions run on a MULTI-CHAIN design.** No such run exists here; the collector refuses it with an error rather than guessing how the `/`-joined FASTA lines up with the fixed list. If a campaign needs it, that is an edit, not a reinterpretation.
- **A PDB with insertion codes, or with chains that `parse_multiple_chains.py` renumbers**, so that the 1-based positions in `fixed_pdbs.jsonl` do not index `seq_chain_<X>` the way assumed. `init_matches_parsed` is what would catch it — the row would error rather than reinsert at the wrong place — but no such input has been run through.
- **Fixed positions on a chain that is not the designed chain**, and the case where `assigned_pdbs.jsonl` is absent while fixed positions exist (`--fixed-positions` requires `--chains-to-design` at submit, so it should be impossible).
- **Whether `atomium_seq_rec` is weighted the way you would expect over the designable subset** — only that its denominator is that subset, inferred from the file format, not read in atomium.py.

### There is no score column — this changes how you rank

AtomiUM's FASTA reports only sequence recovery. **There is no `score` or `global_score`**,
so the model-likelihood ranking that `proteinmpnn_score` supports **does not exist here**.

- **`seq_rec` is not a quality metric, but it is NOT meaningless either** (corrected
  2026-09-30 — this entry previously claimed rfd3 backbones are "effectively poly-glycine"
  and the column close to meaningless, and **that was wrong**). rfd3 emits real side chains
  and real residue identities, so `seq_rec` is a genuine comparison. Measured on 32 designs
  from 8 rfd3 backbones: values **0.21–0.54**, nowhere near zero, and varying **systematically
  by backbone** (one parent 0.209–0.254 across its draws, another 0.433–0.537). It therefore
  appears to track **backbone designability** — how strongly the geometry dictates a sequence —
  rather than sequence quality. Do not rank designs by it, and do not dismiss it; whether it
  predicts anything downstream is still open.
- **`sample` is not unique per row.** AtomiUM numbers samples `n % num_seq_per_target` and
  walks temperatures with `n // num_seq_per_target`, so with several temperatures the same
  `sample` id recurs once per temperature. Use the row name, or the
  (`sample`, `temperature`) pair, to identify a draw.

Rank atomium designs by what comes **after** them — predictor confidence and, above all,
self-consistency (`usalign` back to the parent backbone). That is the honest test for any
sequence designer, and here it is the only one available.
