---
name: atomium
description: >-
  Which sequences does AtomiUM, a private noise-conditioned ProteinMPNN-like designer, give for a backbone? Runs on Modal and vib. Load before composing an atomium run, setting any amino-acid bias, pinning positions with `--fixed-positions`, or feeding an atomium sequence to mkcomplex or a predictor. Traps: `--bias-aa` is global and far stronger than exp(bias) for enrichment; whether `atomium_sequence` carries pinned residues depends on the AtomiUM build, so probe two designs before any fan-out.
---

# atomium

**Custom tool** (`tools/atomium/`), **`action: create`**: mints a child table, one row per designed sequence (`<parent>_a1`, `<parent>_a2`, …; the `_a` suffix distinguishes them from proteinmpnn's `_f` rows). AtomiUM is a private PyG re-implementation of ProteinMPNN-style design with **noise-conditioned weights** (one checkpoint per training noise level). It ships ProteinMPNN's helper scripts and takes the same jsonl inputs, so this tool is deliberately the `proteinmpnn` tool with the same two mini-languages, the same signature-grouping/bin-packing and the same lineage contract: **everything you know about composing a proteinmpnn run transfers** (read that skill for the mini-languages). Repo: private, branch `pure_wo_jit` (`ATOMIUM_REPO` / `ATOMIUM_COMMIT` in `modal_image.py`). On vib, `SAPIA_ACTIVATE_ATOMIUM` points at the lab's shared `atomium.sh` (one GPU per task).

```bash
sapia run atomium <run_dir> -t table0 --num-seq-per-target 2
sapia collect atomium <run_dir> -t <the table the run reserved>
```

`default_input_column` is **`rfdiffusion3_path`** (the bundled proteinmpnn defaults to the older `rfdiffusion_path` and silently submits nothing after an rfd3 run). From another parent pass `-i`.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--model-noise` | `n05` | **The one flag not in proteinmpnn.** Weight file: `n00`…`n07` (training sigma 0.0–0.7 Å) plus the `n05_v2` retrain; validated at submit |
| `--num-seq-per-target` | 2 | Sequences per backbone **per temperature** |
| `--sampling-temp` | `0.1` | One or more space-separated temperatures (`'0.1 0.2'`); see the count rule below |
| `--seed` | 37 | **`0` means a random seed**, not seed zero. The fixed default makes a **rerun deterministic**: re-running the same designs with the same flags returns byte-identical draws, so it is a *reproduction*, never a second sample. To enlarge a set raise `--num-seq-per-target` or pass a different `--seed`; do not count an earlier run's rows plus a rerun of the same rows as independent evidence |
| `--set '--omit_AAs CH'` | none | **A hard omit** (the `--set` help example in `run_atomium.py`; `atomium.py` carries the `omit_AA` path): the listed residues are excluded at every position, unlike the global, attenuated `--bias-aa`. Not yet run here: count the omitted residues in the first output before relying on it |
| `--designs-per-task` | 10 | Designs bin-packed per task; identical-param designs share one batched call |
| `--chains-to-design`, `--fixed-positions`, `--tied-positions`, `--symmetry`, `--bias-aa`, `--set` | | **The flag GRAMMAR is proteinmpnn's**, including the guardrail that positions require a chain list; positions are 1-based within each chain. **Per-design positions:** prosapia's `{...}` islands resolve integers only (an integer, a bare column name, or `+ - * //` of those, looked up the lineage), so a STRING column such as a `fixed_positions` list cannot be substituted; have the upstream tool write one integer column per position and use `--fixed-positions '{fix1},{fix2},{fix3}'` (a NA island fails that design at submit, loudly), or pass a literal. **Whether the OUTPUT holds the pinned residues depends on the build: see "With `--fixed-positions`" below** |

Default Modal resources: **L4**, 8 CPU, 16 GiB, 1 h; weights ship inside the image (no cache Volume).

**Picking `--model-noise`:** the suffix is the coordinate noise the checkpoint was trained to denoise, i.e. how much you trust the input geometry. Low (`n00`–`n02`): trusts the backbone as given (crystal structures, relaxed models). Mid (`n03`–`n05`): the useful range for **de novo backbones** (`n03` is AtomiUM's default; this tool defaults to `n05`). High (`n06`–`n07`): rough or coarsely sampled backbones. Noise level is a cheap axis to scan: run the same parent twice with different `--table-label`s and compare downstream self-consistency, not the sequences.

**The sampling count is a product:** `BATCH_COPIES = num_seq_per_target * len(temperatures)`, each temperature getting a **full** `--num-seq-per-target` set, so `--num-seq-per-target 4 --sampling-temp '0.1 0.2 0.3'` is **12 sequences per backbone**, not 4. There is **no `--batch-size`** (`atomium.py` accepts it but never reads it).

## `--bias-aa` is far stronger than `exp(bias)` for enrichment: calibrate it, never reason about it

`--bias-aa` takes space-separated `AA:bias` pairs and writes one run-wide `bias_AA.jsonl` (`make_bias_AA.py` is purely generative, no structure input).

- **It is GLOBAL, not positional**: one logit offset on every position of every design, buried core included; there is no per-position bias path. `--fixed-positions` does not help (it freezes what the input already has, it does not choose). **To get a residue at a particular place, oversample and select on a measured column.**
- **Enrichment is amplified at low sampling temperature.** Typical (`--model-noise n05`, `--sampling-temp 0.1`, ~2300 sequences per arm, ~80-residue binders): global His frequency 0.615% unbiased (0.49 His per binder) vs **14.2%** at `H:1.39` (~11 His per binder, up to 52). `exp(1.39) = 4` would predict ~2.5%: the shift is **23×, not 4×**, because at T = 0.1 the distribution is nearly deterministic and a modest shift flips the biased residue to the argmax wherever it was competitive. Those two points fit `exp(2.26 × bias)`, so the bias landing the mean at 1.0 is ≈ `H:0.31`, not 1.39. **Fit the multiplier from two runs before trusting a third**; a bias tuned at one `--sampling-temp` does not transfer to another, nor across amino acids.
- **Suppression is attenuated, so the SIGN matters more than the magnitude.** Typical (`C:-1.0`, same settings, ~1100 sequences against an ~11,000-sequence unbiased baseline): Cys frequency 1.41% → 0.774% (**1.8×**; mean Cys per binder 1.11 → 0.584), where `exp(2.26 × bias)` predicted 10× and even plain `exp(bias)` 2.7×; this arm alone fits `exp(0.64 × bias)`. A negative bias only bites where its margin was already thin. **Budget about `exp(0.6 × bias)` for a negative bias and `exp(2.3 × bias)` for a positive one as starting points, and expect to need a filter as well when suppressing.**
- **A bias cannot break a COUPLED choice, only reduce how often it is made.** Cysteines still came in pairs (93.7% of sequences carried an EVEN number, vs 94% unbiased; two Cys stayed five times commoner than one): the bias talked some sequences out of building a disulfide at all and never broke a pair apart. For a structurally coupled residue **gate on the measured count; use the bias only to improve that gate's yield.**
- **A strong bias can make an "exactly one of X" criterion WORSE.** The count of sequences with exactly one copy peaks when the mean is 1.0; pushing past it overshoots (210 designs over 25 backbones unbiased, mean 0.49, fell to 28 over 8 at `H:1.39`). Loose criteria (`>= 1`) improve; strict ones collapse.
- **Verify the bias applied**: `cat <out_dir>/bias_AA.jsonl` and confirm `--bias_AA_jsonl` is in the staged task argv. An unbiased arm has no such file; a no-op bias otherwise gives a normal-looking run that silently duplicates its control.

## Image and build (Modal)

`modal_image.py` clones the private repo at **build** time using the `github-token` and `github-username` Modal Secrets (workspace-level; the build fails without access), then deletes `.git` so the credentialed remote is not baked in. It has to be build time: prosapia's Modal executor hardcodes the task secret to the run's `.env` and exposes no per-tool `secrets()` hook, so a named Secret is unreachable from a running task, and cloning in the task would re-clone on every container. The checkout is **pinned to a commit** (`ATOMIUM_COMMIT`): bump it deliberately, a rebuild will not drift. The first run pays a **long image build inside the submit call** (torch 2.11 + cu128 and ~470 MB of weights): allow 15+ minutes. **`torch_cluster` is a required dependency and easy to drop**: nothing imports it directly, but `model_lib.py` calls `torch_geometric.nn.radius_graph`, a torch-cluster binding; without it every task dies at model construction.

## What it collects

Leaf-prefixed `atomium_`: `sequence`, `sample`, `temperature`, `seq_rec`, `iteration`, plus `_status` and `_path` (the FASTA). **`sequence_full`, `n_fixed`, `fixed_positions`, `fixed_residues`, `init_matches_parsed` and `len_full_matches_parsed` are NOT produced** by this collector (the reinsertion collector is only needed if a build's FASTA omits pinned residues, see below). `atomium_sequence` is what structure predictors consume; it is **not** Boltz's default input column (`proteinmpnn_sequence`), so a Boltz run after atomium needs `-i atomium_sequence`, or `-i mkcomplex_sequence` on a binder table where mkcomplex still belongs in between.

- **The echo row.** AtomiUM's echo of the native input is skipped in either header form, matched case-insensitively on the header's **first comma-separated field**. The real header is **always** `>_init_sequence, chains=<letters>` (`atomium.py` joins the designed chain letters with `/` unconditionally; `chains=` lists the chains AtomiUM designed, i.e. every parsed chain when `--chains-to-design` was not given). An older collector matched the whole header and collected the echo as an `_a0` row per backbone (the parent's input sequence, with `sample`/`temperature`/`seq_rec` all NA). Indices are not renumbered, so designed rows start at `_a1`. **The fix is forward-only**: a collect never deletes rows already in a table and `--force` will not remove them, so a table collected before the fix keeps its `_a0` rows. After `sapia merge` they read `<parent>_a0_<label>`: match `_a0` **as a segment** (not a suffix, never a substring: `_a01` would match) or, more robustly, test `atomium_sample` being empty (true of echo rows and no designed row). Cleaning an existing table is a text-level row drop on the `.tsv`, not a re-collect.
- **`--table-label` alone does not rename a single column.** The leaf prefix comes from `-l/--dir-label`, so `--table-label g8` alone writes `atomium_status`, `atomium_sequence`, … into `table<n>_..._g8`, **not** `atomium_g8_status`. Passing both prefixes every column with the dir label, almost never wanted for a fork.
- **`sample` is not unique per row.** AtomiUM numbers samples `n % num_seq_per_target` and walks temperatures with `n // num_seq_per_target`, so with several temperatures the same `sample` recurs once per temperature: identify a draw by the row name or the (`sample`, `temperature`) pair.

### With `--fixed-positions`: usable, but probe what `atomium_sequence` holds

On the current pin (`ATOMIUM_COMMIT` in `modal_image.py`) **`atomium_sequence` is the full designed chain with the pinned residues in place and no reinsertion step is needed**: the generated length equals the input chain's, `_init_sequence` is byte-identical to the input chain (not the input minus the pinned residues), and the pinned residues sit at the pinned positions (checked with 21 pins on a 4-chain complex with `--chains-to-design D`, and on a ~780-row pinned-His redesign). Applying the reinsertion spec below to such a build would corrupt every row.

**Probe two designs after any pin bump**: the behaviour is a property of the AtomiUM build and **nothing flags it in the table** (rows collect `OK` either way). On the raw FASTA and on `atomium_sequence`, (1) the length equals the input chain's length, (2) the residues at the pinned positions equal the input chain's. **A length short by the number of pinned positions means the old behaviour**: atomium.py then writes only the designable positions, so `sequence` lacks the pinned residues (the very histidines being pinned), `_init_sequence` is the input chain with those positions deleted, `atomium_seq_rec` is computed over the designable positions only (not comparable with a run without fixed positions), and handing `sequence` to mkcomplex and a predictor folds a binder short by N with no pinned residue, a plausible prediction of a molecule nobody designed. That build needs a reinsertion step in the collector (not restored here): read `grp_<g>/parsed_pdbs.jsonl` (`seq_chain_<X>` = the full input chain), `fixed_pdbs.jsonl` (`{design: {chain: [1-based positions]}}`) and `assigned_pdbs.jsonl` (`{design: [[designed chains], [fixed chains]]}`), take each fixed position's residue from the parsed chain and fill every other position **in order** from the designed string (so a descending fixed list like `[49, 12]` cannot scramble it), and **error rather than guess** when the echo is not the input chain minus the fixed positions, the designed string is not `len(chain) − n_fixed`, a fixed position is outside the chain or repeated, a parsed entry is missing, or **more than one chain was designed** (multi-chain reinsertion is not implemented).

Upstream risks no collector here papers over:

- **`_parse_positions` does not check the group count against the designed chains.** `--chains-to-design D --fixed-positions '49/12'` gives two groups for one chain and **silently pins only 49**; check the pinned count in the output against what you meant on the first design of every run.
- **`sequence_recovery` returns `0.0`, not NA, when nothing is designable**: a `seq_rec` of exactly `0.0` can mean "no designable position"; read it beside the chain length and the pin count.
- **The multi-chain `/` split mis-assigns residues when fixed positions precede the last designed chain.**

Not checked: a fixed-positions run on a MULTI-CHAIN design; a PDB with insertion codes or chains `parse_multiple_chains.py` renumbers (1-based positions may not index `seq_chain_<X>` as assumed); fixed positions on a chain that is not the designed chain; and whether `atomium_seq_rec` is weighted as expected over the designable subset.

### There is no score column: this changes how you rank

AtomiUM's FASTA reports only sequence recovery: **no `score` or `global_score`**, so the likelihood ranking `proteinmpnn_score` supports does not exist here.

- **`seq_rec` is not a quality metric, but not meaningless either.** rfd3 emits real side chains and residue identities, so `seq_rec` is a genuine comparison; typical values 0.21–0.54, varying **systematically by backbone** (one parent 0.21–0.25 across its draws, another 0.43–0.54). It appears to track **backbone designability** (how strongly the geometry dictates a sequence), not sequence quality: do not rank by it, do not dismiss it.
- **On a poly-glycine input, `seq_rec` is the glycine fraction and nothing else** (e.g. a graft stub whose free positions really are all glycine: `atomium_seq_rec` equals the designed sequence's glycine fraction outside the pinned positions to five decimals). It is still a true recovery figure, against glycine, and says how often the designer *chose* glycine; it carries no information about design quality. rfd3 output has real side chains, a graft stub does not.

Rank atomium designs by what comes **after** them: predictor confidence and, above all, self-consistency (`usalign` back to the parent backbone), the only honest test for any sequence designer here.
