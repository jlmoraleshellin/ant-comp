# The revised pipeline, and what blocks a Modal test

**Date:** 2026-10-07
**Change:** no tied-positions multi-state design. Design against human only, then
screen for cross-reactivity at the Boltz step.

## The pipeline

| # | Step | Tool | Lands in |
| --- | --- | --- | --- |
| 1 | binder backbones vs human TNF | `laproteina --task-name 90_TNFa_HUMAN_xreact` | `table0` (create, root) |
| 2 | sequences for those backbones | `atomium -i laproteina_path --chains-to-design <B>` | `table1` (create) |
| 3 | rebuild the **human** complex | `mkcomplex -l human --prepend-seqs <human> --repeat 3` | `table1`, col `mkcomplex_human_sequence` |
| 4 | rebuild the **mouse** complex | `mkcomplex -l mouse --prepend-seqs <mouse> --repeat 3` | `table1`, col `mkcomplex_mouse_sequence` |
| 5 | predict human complex | `boltz -l human -i mkcomplex_human_sequence` | `table1`, cols `boltz_human_*` |
| 6 | predict mouse complex | `boltz -l mouse -i mkcomplex_mouse_sequence` | `table1`, cols `boltz_mouse_*` |
| 7 | cross-species gate | a `-f` filter, or just reading the table | — |

Steps 3-6 run the same tool twice on one table. That is supported: `-l/--dir-label`
makes the column leaf `<tool>_<label>`, so the two runs cannot collide. **The label
must be repeated on `collect`.**

### Target sequences for mkcomplex

Taken from the exact structures LPC and Boltz will see, because a sequence that
does not match the structure is the classic off-by-N (every later residue pairing
shifts, and the numbers look plausible and are wrong).

Human, 1TNF chain A, residues 6-157, **152 aa**:

```
RTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLLFAESGQVYFGIIAL
```

Mouse, 2TNF chain A, residues 9-157, **148 aa**:

```
SDKPVAHVVANHQVEEQLEWLSQRANALLANGMDLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVLLTHTVSRFAISYQEKVNLLSAVKSPCPKDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQVYFGVIAL
```

`--repeat 3` emits three consecutive copies — the homotrimer. Prepending puts the
target first, so the **binder is the last chain**, which is also the last chain
index in Boltz's confidence JSON.

---

## The gap: the "simple filter" is not simple

This is the one thing the plan change did not account for, and it decides whether
the campaign can measure what it is trying to measure.

**`boltz_iptm` will not tell you whether the binder binds.** ipTM averages every
interchain pair. The complex here is 4 chains (3 target protomers + binder) = 6
pairs, of which **3 are target-target** — TNF trimerising with itself, which it
will do beautifully whatever the binder does. Half the number is noise for this
purpose, and the half that matters is not separable from it.

The boltz skill is explicit, from a measured campaign:

> use `confidence_score` as a gate for "did the target assemble", then stop using
> it. It carries no binder signal once past that gate.

A filter on `boltz_human_confidence_score >= X and boltz_mouse_confidence_score >= X`
would therefore select **designs where TNF folded nicely in both species** — which
is every design — and say nothing about cross-reactive binding.

### What actually carries the signal

`pair_chains_iptm[binder][target]` — and **prosapia does not collect it.** It is in
the per-prediction JSON on disk:

```
<out_dir>/boltz_results_shard_*/predictions/<name>/confidence_<name>_model_0.json
```

- `chains_ptm` — keyed by positional chain index; the binder is the last one.
- `pair_chains_iptm` — full N x N, **not symmetric** (observed asymmetry up to
  0.079), so symmetrise by averaging both directions.

Reference points from the earlier SlyB campaign: a believed interface is roughly
**0.5-0.6** binder-target ipTM and **< 10 A** interface PAE. That campaign's best
were 0.357 and 19.5 A, and it concluded no validated binder.

### Three ways to close it

- **(a) A small `update` tool** that reads the confidence JSON and writes
  `pair_iptm_binder_target` (plus `binder_ptm`) as columns. ~1 file of real work,
  mirrors `chainsel`'s shape. Then the cross-species filter is genuinely simple:
  both labels above threshold.
- **(b) The binder-pLDDT delta.** Compare the binder's mean pLDDT alone vs. in the
  complex; a real interface raises it. The SlyB campaign found this *more*
  informative than any ipTM flavour (all 27 designs lost 7.5-23.3 points). But it
  needs an extra Boltz run on the binder alone, so it is not cheaper.
- **(c) Ship without it** and read the JSONs by hand for the handful of survivors.
  Fine for a smoke test, useless at scale.

**Recommendation: (c) for the smoke test, (a) before any real batch.** The smoke
test is about whether the machinery runs at all; the metric only has to be right
before it decides something.

---

## Blockers

### Needs you

1. **Modal is not configured.** No `~/.modal.toml`, no tokens in `.env`. Nothing
   runs until `modal setup` (browser login).
2. **atomium's Modal Secrets.** Its image clones a **private** repo at build time
   using workspace-level Secrets `github-token` and `github-username`. If you are
   not in that Modal workspace the image build fails and **step 2 is dead** — and
   there is no workaround from here, because prosapia's Modal executor exposes no
   per-tool secrets hook. Confirm with the repo owner. `proteinmpnn` is the
   fallback if the answer is no.

### Needs Modal first, then I can do it

3. **`lpc-checkpoints` Volume** — `complexa.ckpt` (2.7 GB) + `complexa_ae.ckpt`
   (3.8 GB). Public on NGC, no account. Pull them *inside* a Modal container: this
   network intercepts TLS to Modal blob storage, so `modal volume put` from here is
   unreliable.
4. **`lpc-assets` Volume** — `configs/targets/targets_dict.yaml` plus
   `data/target_data/tnf_xreact/`. The image clones LPC fresh, so without this the
   target does not exist in any container.

### Unknown until the first run

5. **The binder's chain id** in LPC's output. Needed for `atomium
   --chains-to-design`, and to know which index is the binder in Boltz's JSON. One
   collected PDB answers it.

---

## Two decisions before spending

### MSA policy

`--use-msa-server` is off by default. For the **de-novo binder** that is correct —
no MSA exists. For the **TNF target it is wrong**: TNF has a deep MSA, and denying
it one weakens both target assembly and interface confidence, so good designs read
as bad ones.

The previous campaign left this unresolved and flagged it as the confound that
undermined its own negative result. Do not repeat that. Note it calls an external
server, so it needs your explicit go-ahead.

### Cost

Boltz is the expensive step and this plan runs it **twice per design**, on a
~600-residue 4-chain complex. Against a $50 budget that is the line item to watch.
Start with a handful of designs and measure before scaling.

---

## Smoke test, when unblocked

```bash
sapia new_run --label tnf_smoke

sapia run laproteina <run_dir> --task-name 90_TNFa_HUMAN_xreact \
    --num-jobs 1 --nsamples 2 --design-prefix tnfh
sapia collect laproteina <run_dir> -t table0
# -> inspect one PDB: confirm the binder's chain id before going further
```

Everything after that waits on what the backbones look like.
