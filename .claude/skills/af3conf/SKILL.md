---
name: af3conf
description: >-
  Does AlphaFold3 believe the binder's fold and its interface, per binder rather than per whole complex? Reads the JSONs an alphafold3 leaf already wrote and writes binder pLDDT, `chain_ptm`, symmetrised `chain_pair_iptm` and ipSAE under boltzconf's column names. Load before composing an af3conf run, comparing `af3conf_*` with `boltzconf_*`, or trusting a per-chain AF3 number. Trap: the summary's chain arrays carry no chain letters, so read `chain_pae_min_maxdiff` and `chain_pae_min_margin` first.
---

# af3conf

Per design, the **binder-specific** confidence of one AlphaFold3 prediction: binder pLDDT, the binder's own `chain_ptm`, the symmetrised binder–target `chain_pair_iptm` (overall and per target chain) and **ipSAE** (Dunbrack's interface-restricted pTM, per target chain and in both directions). **`action: update`**: annotates the table it reads, in place; `-t` is required. CPU only (the manifest builder forces `gpus_per_task = 0`, so no `-g 0`). It is the AF3 counterpart of `boltzconf` and runs on finished `alphafold3` leaves, with no re-prediction.

## Why it exists

`alphafold3` collects whole-complex numbers only. Its `iptm` is diluted by chain count exactly as `boltz_iptm` is: on three target copies plus one binder, most averaged pairs are target–target pairs the binder has nothing to do with. AF3 writes the binder-specific numbers (`chain_ptm`, `chain_pair_iptm`, the PAE matrix, per-atom pLDDT) in raw files beside the structure that nothing else collects.

## Premise, and when its numbers are meaningless

Every assumption below is **checked and errors rather than guessing**, except the two marked unchecked.

| Assumption | Signature when it fails |
| --- | --- |
| The input column is AF3's own `<prefix>_model.cif` | `error: … is not named <prefix>_model.cif` (a non-mmCIF is refused; a backbone or `chainsel` output lands here) |
| `<prefix>_confidences.json` and `<prefix>_summary_confidences.json` sit beside it (AF3 drops the `_model`) | `error: <file> is missing beside the structure`; never a partial row or another model's number |
| One token per residue | `error: N tokens but M residues …`. A **multi-atom ligand or modified residue anywhere** (named chain or not) is tokenised per atom and refuses the row, N − M = ligand atoms − 1. A one-atom ion is one token and passes here |
| Every token is the CIF residue of the same index | `error: token i is chain X … but residue i … is chain Y` or `… numbered differently`, element by element on chain and `label_seq_id`. Every residue of a **named** chain needs an integer `label_seq_id`, which is what refuses a named ion (`… NAMED chain C has no integer label_seq_id`; relies on AF3 writing `.` for a non-polymer, standard mmCIF but not checked on real ligand output) |
| The named chains are proteins | `error: chain(s) B (has DA) … look like a nucleic acid`. A nucleotide is one token, so the token checks pass on DNA/RNA, and without this check ipSAE would use the protein `d0` floor of 1.0 instead of 2.0 and come back **inflated with every trust column clean**. An unnamed nucleic-acid chain is fine |
| Exactly one binder chain | refused at submit (a two-chain binder could only be an average that reads like one chain) |
| The summary's positional chain index *k* is the *k*-th chain of `token_chain_ids` | `error: chain_pair_pae_min in the summary disagrees by up to X A …`; can be blind, see the trust columns |
| `chain_ptm` and `chain_pair_iptm` follow the order of `chain_pair_pae_min` | **Unchecked**: only `chain_pair_pae_min` can be recomputed |
| The summary describes the same complex | `error: 'chain_pair_iptm' is K x K but the prediction has N chain(s)` (likewise `chain_ptm`, `chain_pair_pae_min`); `error: the summary JSON has no 'ranking_score'` |
| Every value is finite | `error: … holds N non-finite or null value(s)` (a NaN would silently defeat the mapping checks, since `nan > tol` is False) |
| The rows of `-t` are the molecules that were predicted | the submit **raises** if the structure column (or `--sequence-column`) is absent from `-t`, unless `--allow-lineage`: a child table would otherwise score its parent's prediction under its own name |

It does **not** know whether the chain you named is the binder, whether the pose is right, or whether the binder binds: every column is the predictor's belief about its own output, and a binder confidently docked on the wrong face gets a high ipSAE. A wrong-but-present `--binder-chains` letter gives a clean row about the wrong pair, and **only `binder_n_res` against the DESIGNED binder length shows it**. `binder_seq_match` does not: it compares the named chain with the same letter's segment, so a target letter passed as the binder matches itself. Keep the `binder_n_res` length gate in every filter, with or without `--sequence-column`.

## Trust columns: read these before the science

| Column | Expected | Read it as |
| --- | --- | --- |
| `chain_pae_min_maxdiff` | **≤ 0.05 Å** (rounding alone allows 0.055); > 0.1 is an error row | The check on the chain index → letter map: `chain_pair_pae_min` recomputed from the full PAE over the `token_chain_ids` blocks, compared under the identity order. Blind to any permutation the `chain_pair_pae_min` matrix is invariant to within 0.1 Å (identical copies in a symmetric pose, but also chains whose minimum PAEs happen to agree) |
| `chain_pae_min_margin` | **> 0.1 Å** for the check above to mean anything | Smallest maxdiff over every transposition of two chains. **≤ 0.1 means the order check could NOT have caught a swap on this row**; not an error. It is not only a symmetric-target issue: on any two-chain complex it is the binder↔target swap, `max(\|R00−R11\|, \|R01−R10\|)` over the block minima R, and on a confident interface the two off-diagonal minima can sit within PAE's 0.1 Å resolution. What stays unchecked then is `chains_ptm_binder` (could be the target's `chain_ptm`) and, with several target chains, which `iptm_<x>` is which copy (the symmetrised two-chain `iptm` is unchanged by the swap). Bounded to transpositions |
| `binder_seq_match` | `1` (with `--sequence-column`) | `1` when the CIF chain named by `--binder-chains` equals the **same letter's** segment of the sequence column; `0` is a trust value, not an error. Catches a prediction of a different design at this path; does not say which chain is the binder. **NA = unchecked**. Kept on error rows when measured before the error |
| `ranking_score` | equals `alphafold3_ranking_score` | The summary's own score, echoed. Weak: it catches a prediction regenerated at the same path only if the score moved by ≥ 0.01, and never a stale shard (alphafold3's collector reads the same stale summary) |
| `plddt_bfac_maxdiff` | **0.0**; > 0.05 is an error row | Max \|`atom_plddts` − CIF B-factor\| over every atom; with the element-wise `atom_chain_ids` check it validates the atom order `binder_plddt` reads through. Says nothing about the chain-index map |
| `n_tokens`, `n_res_struct` | equal on every OK row | A mismatch is an error; present so the equality is visible |
| `binder_n_res` | the designed binder length | **The only catch for a wrong `--binder-chains` letter**; two chains of the same length defeat it, so know the binder letter |
| `n_target_res` | the exact expected target total | AF3 writes every residue it was given, so any other value is a wrong input or `--target-chains` |
| `stem_matches_name` | `1` | `0` = AF3 sanitised the name or `-i` points at a per-sample model; not an error. **`1` does not prove identity**: two names that sanitise to one string, or a stale prediction in a later-sorting `results_shard_*`, leave it clean. `binder_seq_match` catches only the collision; **no af3conf column catches a stale shard** of the same design, and the defence is a fresh alphafold3 `-l` per MSA, template or seed policy change |
| `input_source` | `own table` | `lineage` only with `--allow-lineage`: an ancestor's prediction was scored |
| `chains` | `binder\|targets` (`D\|A,B,C`) | Echo of what was asked; kept on error rows |

Typical, on a real AF3 prediction with no confident interface: the identity order agrees to 0.05 Å, rotating the chain order by one puts it 6–13 Å off, and the margin is a couple of Å. **On a confident complex the margin can be far smaller**, which is why it is a column. One `chain_pair_pae_min` error in a batch means distrust the batch's OK rows: a systematic order change errors loudly on rows with a large margin and passes rows at or under the tolerance, which are often the confident ones.

### How the mapping differs from boltzconf

- AF3 names chains: `confidences.json` gives the letter of every token and atom, checked element by element against the CIF. Boltz's arrays carry no letters, so boltzconf can only check the token→residue order via pLDDT against B-factors.
- Both summaries are positional (`chain_pair_iptm` is a plain list of lists). boltzconf checks only the chain count; af3conf checks `chain_pair_pae_min`'s order on every row with margin > 0.1 Å and assumes the order of `chain_ptm`/`chain_pair_iptm`. **Also run the boltzconf batch-level ranking check**: `iptm_<x>` and `ipsae_max_<x>` come from different files and should rank the target chains the same way.
- `binder_plddt` is the mean **CA-atom** pLDDT (the ipsae.py AF3 convention); boltzconf's is Boltz's per-token value. Same kind of number, not the same definition.
- AF3 splits a multi-letter entity into plain chains (`ABC` → A, B, C), so `--target-chains A:C` names them.
- `chain_pair_iptm` is symmetrised, a no-op on real AF3 output (exactly symmetric), kept so `iptm_<x>` has boltzconf's definition. pLDDT is 0–100; a 0–1 array is an error, not a rescale.

## ipSAE: a verbatim port of boltzconf's

`ipsae_d0`, `ipsae_asym`, `reject_nucleic_acid` and their constants are source-identical copies of `tools/boltzconf/boltzconf_worker.py` (reference: `DunbrackLab/IPSAE/ipsae.py`, md5 `a48df7adc64afa36b4b475429b0a8aff`; provenance in boltzconf's worker docstring). `test_af3conf.py` asserts the copies are source-identical and give identical output on one PAE matrix; change the arithmetic in both tools together. The arithmetic: `pae[i][j] < cutoff`, **strictly, on the raw asymmetric matrix**; `d0` from `n0res` with L **raised to 26 first** and floored at 1.0; the chain-pair value is a **max over residues of a per-residue mean**, asymmetric, with `ipsae_max_<x>` the larger direction. Differences from the reference (same as boltzconf): **NA instead of 0.0** for a pair with no sub-cutoff residue pair, and counts by token index. Agreement with the reference pins file layout and arithmetic, not behaviour on a binding pose. `ipSAE_d0chn`, `ipSAE_d0dom`, `ipTM_d0chn`, pDockQ, pDockQ2 and LIS are not collected.

## Questions, and the bar for each

| Question | Column | Bar |
| --- | --- | --- |
| Is the binder a confidently folded domain? | `binder_plddt`, `binder_plddt_min`, `chains_ptm_binder` | Mean pLDDT 80+ is the usual reference point for a believable fold. Independent of placement |
| Does AF3 believe the interface? | `iptm_binder_max` (best target chain), `iptm_binder_top2` (two best), `iptm_binder_mean` (all) | **`_max` or `_top2` for a binder meant to span protomers**; `_mean` is diluted by every target chain the binder never touches. With one target chain `_top2` equals `_max`. No calibrated bar exists for AF3 here |
| Same, normalised for interface size? | `ipsae_binder_max`, `ipsae_binder_top2`, `ipsae_max_<x>`, always with `ipsae_n0res_*` | **No published pass mark**: the paper's 0.3 is the bottom of the 0.3–0.7 overlap band, not a cut, and non-interaction is 0.0–0.1. Take the bar from your own batch's distribution |
| Which side is AF3 confident about? | `ipsae_binder_to_<x>` vs `ipsae_<x>_to_binder` | ipSAE is asymmetric; identical values on every row would mean the raw matrix was not used |

### Comparability keys

- `ipsae_pae_cutoff`: ipSAEs at different cutoffs are different quantities; never pool them, and give a rerun at another cutoff its own `-l`.
- `--target-chains` (read it off `chains`): `ipsae_binder_mean`'s denominator and `ipsae_binder_top2`'s width are the number of named target chains, with an unreached chain contributing a hard 0.0, so `A:C` and `A` runs do not share a distribution of the aggregates. Per-chain columns are unaffected.
- af3conf vs boltzconf: same names and arithmetic on a different predictor's PAE and pLDDT. Compare **as ranks within one batch**, never as one number scale, and only with the same cutoff and target-chain set. AF3 stores PAE to 0.1 Å (boltz's npz is full precision), which can move a pair within 0.05 Å of the cutoff across it.

## Invocation and flags

```bash
# Smoke-test ONE design first. -w filters on a COLUMN (-w name=... does not work): use that design's exact alphafold3_path.
sapia run af3conf <run_dir> -t <table> -i alphafold3_path --binder-chains D --target-chains A:C \
    --sequence-column mkcomplex_sequence -l verification -w alphafold3_path=<exact path of one design>
sapia collect af3conf <run_dir> -t <table> -l verification

# The whole table, cutoff stated rather than inherited:
sapia run af3conf <run_dir> -t <table> -i alphafold3_path --binder-chains D --target-chains A:C \
    --sequence-column mkcomplex_sequence --ipsae-pae-cutoff 15.0
sapia collect af3conf <run_dir> -t <table>
```

| Flag | Default | Note |
| --- | --- | --- |
| `-t/--table` | none | **Required** (`update`) |
| `-i/--input-column` | **`not applicable` sentinel, so effectively required** | AF3's **predicted structure** column (`alphafold3_path` or `alphafold3_<label>_path`), not a sequence column |
| `--binder-chains` | none | **Required, exactly one chain.** AF3 letters are positional in the sequence it was given; after `mkcomplex` the binder is the **last** chain |
| `--target-chains` | none | **Required.** Chain mini-language (`A,B,C`, `A:C`), disjoint from the binder. **The order is the order of the per-chain columns** |
| `--sequence-column` | **unset, so `binder_seq_match` is NA** | The `/`-chainbreak column alphafold3 ran on, normally `mkcomplex_sequence`. Segment *i* is chain `chr(ord('A') + i)` (`prosapia.utils.select_chains`, the rule `alphafold3` letters its chains by). Resolved own-table-first, lineage only with `--allow-lineage` (unlike epitope, which always walks it). An unresolvable sequence or missing binder segment is an **error**; a mismatch is `binder_seq_match = 0`. If alphafold3 ran with `--positions` slicing the binder, the CIF binder is shorter than the segment and the **correct** chain gives `0` |
| `--ipsae-pae-cutoff` | `15.0` | Å; a pair counts when `pae[i][j] <` this, strictly. Upper end of the paper's "10 or 15 Å" band, same as boltzconf; echoed into `ipsae_pae_cutoff` |
| `--allow-lineage` | **off** | Lets structure and sequence columns come from an ancestor table; then read `input_source` and `stem_matches_name` |
| `--designs-per-task` | `20` | |
| `-l/--dir-label` | empty | Own output dir and column prefix, for a second alphafold3 leaf or a second cutoff. **A label with any `_`-separated segment equal to or starting with `iptm`, `ipsae` or `ipae` (case-insensitive) is refused at submit** (the collector resets per-chain columns by prefix, so `af3conf_x` would read and blank a sibling `-l x_iptm` leaf's columns) |

No `--contact-cutoff`: the contact-restricted `ipae_*`, `n_if_*` and `pae_binder_self` are boltzconf-only; for contacts use `epitope`. Resources: 2 CPU, 4 GiB, 30 min per task; Modal image is `debian_slim` + `gemmi==0.7.5` + `numpy==2.5.3` (boltzconf's pins). On vib, `SAPIA_ACTIVATE_AF3CONF` points at the shared `prosapia_python.sh`; an unset entry shows as `sapia: set SAPIA_ACTIVATE_AF3CONF in your .env` in `.err`.

## Gotchas, with their signatures

- **Siblings are found from the file, not the design name.** Every row `error: … is not named <prefix>_model.cif` or `…_confidences.json is missing beside the structure` means the wrong column was passed.
- **Only the top-ranked model is read** (files at the design dir root, copied from the best `seed-*_sample-*`); pointing `-i` at a per-sample CIF works but gives `stem_matches_name = 0`.
- **Stale OK after the alphafold3 leaf changes underneath.** af3conf rows stay `OK` when their prediction is replaced. After any in-place alphafold3 rerun, run AND collect af3conf with `--force`, or use a fresh `-l`: `--force` on the run alone is not enough because `sapia collect` without `--force` drops rows already `OK`. Filters must also require `alphafold3_status == "OK"` and the path equality below.
- **`af3conf_path == alphafold3_path` selects ZERO rows on Modal.** `af3conf_path` is absolute (`/runs/outputs/…`), `alphafold3_path` relative (`outputs/…`); on vib, where run_dirs are absolute, they are equal. Compare after stripping the `/runs/` prefix (helper in the filter below). **Do not use `os.path.abspath`**: it resolves against the filter process's cwd, which on Modal is the volume's real mount, not `/runs`, so every row fails while the filter looks correct.
- **A blank own input cell keeps whatever af3conf status the row already had, including a stale `OK`** (`filter_ready` drops the row from collect). Signature: `alphafold3_path` blank or `alphafold3_status` not `OK` on a row whose `af3conf_status` is `OK`; the `alphafold3_status` gate removes it.
- **`ipsae_binder_top2` and `_mean` carry the reference's 0.0 for an unreached chain, while `ipsae_max_<x>` for that chain is NA.** Signature: `_top2 ≈ _max / 2` on a multi-copy target means the binder reached one protomer, not two. All three aggregates are NA only when no target chain has a qualifying pair.
- **`d0res` saturates at 1.0 for every `n0res ≤ 26`**: `ipsae_d0res_max_<x> == 1.0` means "small interface"; read `ipsae_n0res_*` for the size.
- **Memory on large complexes.** `json.load` materialises `pae` and `contact_probs` as Python lists and the peak is per design: about 150 MB at 1100 tokens, n²-dominated (extrapolated: ~2 GB at 5k tokens, the 4 GiB limit near 7k). Signature of an OOM: every design in that task collects `missing`. Lower `--designs-per-task` does not help; raise `--mem`.
- **Rerun semantics.** `sapia run … --force` deletes the per-design TSVs first, so a killed task collects as `missing`, never the old `OK`. `sapia collect` without `--force` skips rows already `OK`, so a rerun at a new cutoff keeps the old numbers unless the collect also gets `--force`; prefer a fresh `-l` per alphafold3 leaf or cutoff.
- **Per-chain columns are reset on every collect** (every `af3conf[_<label>]_iptm_*` / `_ipsae_*` column of this leaf is blanked first, including on `missing`/`error` rows and on a rerun with fewer target chains), which is why those label stems are refused.

## What it collects

Prefix `af3conf_` (or `af3conf_<label>_`); trust columns are in the table above. `_status` is `OK` / `error: …` / `missing`: **`OK` is the only proof of success**, and only together with `alphafold3_status == "OK"` and the path equality. Every other row is NA in every numeric column (`binder_seq_match` aside). `_path` is the scored structure (absolute; empty on non-OK rows). Per-chain columns are dynamic: the collector reads their names from each design's TSV header, which sits at `<run_dir>/<table>/af3conf[_<label>]/<name>.tsv`.

| Column | Meaning | How to read it |
| --- | --- | --- |
| `binder_plddt` | mean CA pLDDT over the binder, 0–100 | Fold confidence; blind to placement |
| `binder_plddt_min` | worst binder residue's CA pLDDT | A frayed terminus shows here while the mean stays high |
| `chains_ptm_binder` | AF3's `chain_ptm` for the binder, 0–1 | Named after boltzconf's column; its chain order is assumed, not checked |
| `iptm_binder_mean` / `_max` / `_top2` | mean / highest / mean of two highest per-chain symmetrised binder–target `chain_pair_iptm` | See the questions table; `_mean` diluted, `_top2 = _max` with one target chain |
| `iptm_<x>` | symmetrised `chain_pair_iptm`, binder vs chain *x* | One per `--target-chains` entry, lower-cased; never NA on an OK row; unchecked among chains the margin cannot separate |
| `ipsae_pae_cutoff` | the cutoff used, Å | Comparability key |
| `ipsae_binder_max` | highest `ipsae_max_<x>`, 0–1 | **The headline filter column**; NA only when no target chain has a qualifying pair |
| `ipsae_binder_top2` / `ipsae_binder_mean` | mean of the two highest / of all `ipsae_max_<x>` | Contain the reference's 0.0 for an unreached chain; `_mean` is more diluted than `iptm_binder_mean`: **never gate on it** |
| `ipsae_binder_to_<x>` / `ipsae_<x>_to_binder` | ipSAE with the frame on the binder / on chain *x* | NA when no pair is under the cutoff, never 0.0 |
| `ipsae_max_<x>` | the larger direction | NA when both are NA |
| `ipsae_n0res_binder_to_<x>` / `ipsae_n0res_<x>_to_binder` | SUPPORT: `n0res` of the best residue | A high score on `n0res = 2` is one residue, not an interface; float; NA with the score |
| `ipsae_d0res_max_<x>` | SUPPORT: d0 of the winning direction, Å | 1.0 for every `n0res ≤ 26` |
| `ipsae_nres_binder_<x>` / `ipsae_nres_<x>` | SUPPORT: binder residues with ≥ 1 sub-cutoff PAE to *x* / chain-*x* residues reached | Real counts; 0 is a real zero and is exactly when `ipsae_binder_to_<x>` is NA |

**NA is never 0.** Values that can legitimately be a real 0 are `binder_seq_match`, `ipsae_nres_*` and the binder-level aggregates. Reading an ipSAE NA as 0 is legitimate only on `OK` rows: mask on `af3conf_status` first, as the boltzconf skill shows.

## Filters

Take cuts from your own batch's distribution (`table-analyst` gives quantiles). Both filters gate identity first, and the `binder_n_res` gate is the only catch for a target letter passed as the binder.

```python
# filters/af3_ipsae_gate.py
def _abs(col):
    # af3conf_path is '/runs/outputs/...', alphafold3_path 'outputs/...' on Modal: strip the prefix.
    def one(p):
        if not isinstance(p, str) or not p:
            return None
        return p[len("/runs/"):] if p.startswith("/runs/") else p.lstrip("/")

    return col.map(one)


def apply_filter(df):
    ok = (
        (df["alphafold3_status"] == "OK")
        & (df["af3conf_status"] == "OK")
        & (_abs(df["af3conf_path"]) == _abs(df["alphafold3_path"]))
        & ((df["af3conf_ranking_score"] - df["alphafold3_ranking_score"]).abs() < 1e-6)
        & (df["af3conf_binder_seq_match"] == 1)            # NA (unchecked) fails here
        & (df["af3conf_input_source"] == "own table")
        & df["af3conf_binder_n_res"].between(40, 60)       # the DESIGNED binder length
        & (df["af3conf_ipsae_pae_cutoff"] == 15.0)         # one cutoff per distribution
    )
    good = df[ok]
    bar = good["af3conf_ipsae_binder_max"].quantile(0.75)
    return good[
        (good["af3conf_ipsae_binder_max"] >= bar)
        & (good["af3conf_ipsae_n0res_binder_to_a"] >= 8)   # support behind the score
        & (good["af3conf_binder_plddt"] >= 70.0)
    ]
```

If `chains_ptm_binder` or a per-chain value (`iptm_a` vs `iptm_b`) decides the selection, also require `af3conf_chain_pae_min_margin > 0.1`, or accept that the chains it names may be swapped. To keep what both predictors like, rank within the batch and combine: require both leaves `OK`, the same `ipsae_pae_cutoff` and `chains`, `binder_n_res` in the designed band for both, then `rank(pct=True)` of `af3conf_ipsae_binder_max` and `boltzconf_ipsae_binder_max` both ≥ 0.75. Substitute labelled column names (`alphafold3_<label>_path`, …) for a labelled leaf.

## Invariants a worker can check

Against numbers af3conf did not compute, on every OK row of a mkcomplex complex:

- `len(segment[binder letter])` of `mkcomplex_sequence` == the CA count of that chain in the CIF == `binder_n_res`, and that is the designed length; the sum of all segment lengths == `len(token_chain_ids)` in `<prefix>_confidences.json` == `n_res_struct` (no unnamed chains).
- `ranking_score` == `ranking_score` in the raw `<prefix>_summary_confidences.json` beside `alphafold3_path`; `af3conf_path` and `alphafold3_path` differ only by the `/runs/` prefix.
- On a design with a confident interface, recompute from the raw JSONs max \|summary `chain_pair_pae_min` − block minimum of `pae`\| under the identity order (≤ ~0.055, equals `chain_pae_min_maxdiff`) and the smallest over every two-chain transposition (equals `chain_pae_min_margin`; > 0.1 for the order check to have meant anything).
- **Break it on purpose:** a **target** letter as `--binder-chains` must give `binder_seq_match = 1`, a clean `OK`, and `binder_n_res` equal to a target chain's length (only the length gate catches it); a non-AF3 structure column must give `error:` with NA on every row, never a number.

## What it is blind to

- **Whether the binder binds.** Agreement between two predictors is two beliefs, not evidence.
- **Where it bound.** A confident pose on the wrong epitope scores well: `epitope`. Whether it kept the reference mode: `dockq`; self-consistency against the parent backbone: `usalign`.
- **Contacts, interface area, shape, energy.** No contact-restricted PAE or count here (`epitope`, boltzconf); `cms` (SC, CMS) and `pyrosetta` (`if_dG`, `if_dSASA`); neighbouring-protomer occupancy: `ringfit`.
- **Identity.** Without `--sequence-column`, nothing ties the structure to this design's sequence; with it, only for the named chain, and no column catches a stale shard of the same design.
- **Chain swaps the PAE minima cannot separate** (`chain_pae_min_margin ≤ 0.1`), and the order of `chain_ptm`/`chain_pair_iptm` everywhere.
- **Other AF3 samples and seeds** (only the top-ranked model is read), the reference's other scores, and AF3's whole-complex numbers and per-chain `chain_iptm` (`iptm`, `ptm`, `has_clash` come from `alphafold3`).
