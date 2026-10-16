---
name: posrmsd
description: >-
  Did a prediction keep specific residues where its backbone put them? posrmsd Kabsch-fits
  structure A (the prediction) onto structure B (the reference/backbone) over an EXPLICIT list of
  1-indexed FILE-ORDER chain positions, such as motifrmsd_design_fixed_positions. It writes the
  CA RMSD, a per-segment RMSD from that one fit, and a residue-identity trust check. Load before
  composing a posrmsd run or reading any posrmsd_ column. The trap: it assumes A and B share
  residue ordering and never aligns anything, so a wrong chain or shifted numbering still gives a
  plausible RMSD. Read seq_match_frac first. Covers lineage resolution of --col-b and
  --positions-column, the A-onto-B order, the prefix contract shared with usalign,
  missing_input vs missing vs error, and why positions are not PDB residue numbers.
---

# posrmsd

**Custom tool** (lives in `tools/posrmsd/`, not bundled with prosapia).
**`action: update`**: it annotates the `-t` table in place, and `-t` is required.
**CPU-only**: the manifest builder forces `gpus_per_task = 0`, so you do not pass `-g 0`.

For each row it takes structure A (`--col-a`) and structure B (`--col-b`, resolved up the lineage).
It picks **one chain in each** and takes the CA atom at every position in the spec. A position
means the Nth polymer residue of that chain **in file order**. It then **superposes A's CAs onto B's
CAs** (SVD/Kabsch, reflections excluded) and reports the CA RMSD after that single fit.
**B is the reference frame; A is moved.** Pass the prediction as A and the backbone as B.

## Premise, and when its numbers are meaningless

**The premise: A's chain and B's chain have the same residue ordering, so position N in A is
the same residue as position N in B.** The correspondence is never searched for, aligned or
checked against sequence. That holds by construction when A is a structure prediction of a
sequence that was designed on B's backbone, for example Boltz of a ProteinMPNN sequence against
the rfd3 binder chain it came from: same length, same order, no insertions.

| Situation | What happens |
| --- | --- |
| A predicted from a sequence designed on B (same length, same order) | the measurement means what it says |
| A and B differ in length on the chosen chains | `warn: chain lengths differ …` with the numbers kept. **`== "OK"` excludes the row.** The premise is not guaranteed |
| A and B have equal length but different ordering (wrong chain, domain swap, unrelated construct) | **it returns a plausible RMSD for the wrong residue pairs.** Only `seq_match_frac` and the two `*_motif_seq` columns catch this |
| the position spec was written in PDB residue numbers on a chain not numbered 1..N | **it measures the wrong residues.** `*_resnum_offset` ≠ 0 or `*_gapless` = 0 flags it |
| a position past the end of either chain | `error: … has N polymer residues but the spec asks for position(s) …`, everything NA |
| a structure with several polymer chains and no `--chain-a`/`--chain-b` | `error: … has 4 polymer chains (A, B, C, D) and --chain-b was not given …`. **A chain is never guessed** |

**Use usalign instead** to compare two structures with no known residue correspondence; it finds
its own alignment. **Use motifrmsd instead** to ask whether rfd3 rebuilt a grafted motif against
the original reference. That correspondence comes from rfd3's `diffused_index_map`, not from shared
ordering. posrmsd fills the gap between those two: a known, shared ordering plus an explicit
position list.

**Positions are FILE-ORDER positions, not PDB residue numbers.** This is the same position
language as ProteinMPNN's `--fixed-positions` and `motifrmsd_design_fixed_positions`: 1-indexed,
inclusive `start:end` ranges or single positions, comma-separated, **one chain group** (a `/` is
refused). The spec is parsed at submit time by prosapia's own `parse_chain_positions`. The
concatenated result is then checked against `parse_positions` on the same string, so the semantics
cannot drift from ProteinMPNN's. `{expr}` islands resolve up the lineage. Waters, ions and ligands
in the chain are not counted. Polymer residues are counted, including non-standard ones such as MSE.

## Trust columns: read these before the RMSD

Every column is `posrmsd[_<dir-label>]_<prefix>_<name>` (see "Prefix contract" below). The
examples below use the full name for the default prefix `boltz_vs_chainsel_binder`.

| Column (`…_<prefix>_`) | Expected in the designed-on-B case | How to read it |
| --- | --- | --- |
| `status` | `OK` | **The only proof of success.** `warn:` = numbers computed but the chain lengths differ. `error:` = nothing computed; the message names the structure (A/B), the chain and the index. `missing_input:` = not submitted, because `--col-b` or `--positions-column` resolved to nothing on this row and all its ancestors. `missing` = no result file (not selected by `--col-a`, still running, or killed) |
| `seq_match_frac` | **1.0** when ProteinMPNN held these positions fixed and the indexing is right | fraction of positions where A's one-letter residue equals B's. **`X` (unknown residue) never counts as a match**, so a poly-UNK file reads 0 rather than a false 1. A value below 1.0 has two causes: ProteinMPNN changed a residue it was told to keep, or the indexing/chain choice is wrong. The motif-seq pair tells them apart. **`seq_match_frac == 1.0` does NOT certify the spec points at the graft** — it only certifies that A and B are indexed the same way and that B's residue at each spec position survived into A unchanged. If the position spec itself is wrong (off by one, or pointing at a non-graft window), ProteinMPNN still "holds" whatever B has there, A copies it, and `seq_match_frac` reads a clean 1.0 while the real graft residues elsewhere were redesigned. **To confirm the spec itself is the graft, cross-check `b_motif_seq` (commas stripped) against the parent row's `motifrmsd_ref_motif_seq` (commas stripped)** — that is the one column tying the position spec back to the actual reference motif identity, independent of what ProteinMPNN did |
| `a_motif_seq` / `b_motif_seq` | identical | the residues actually read, segments comma-joined in spec order (`YWVNPTLT,NTNGLITKYP,…`). **One or two scattered differences** suggest a held-position violation. **Every segment different, or the same letters shifted by one** means wrong chain or an off-by-one in the spec |
| `a_chain` / `b_chain` | the binder letters you meant | the chain actually indexed in each file (gemmi's chain name: `auth_asym_id` for mmCIF, the chain column for PDB) |
| `a_len` / `b_len` | equal | polymer residues on each chain. If they differ, the status is `warn:` |
| `n_positions` | the spec's size | always the full spec on OK/warn rows; a shortfall is an error, never a smaller fit |
| `a_resnum_offset` / `b_resnum_offset` | 0 for rfd3/Boltz output | residue number minus file position, if constant over the selected positions. **NA = not constant** (a numbering gap or insertion code inside the span); that is not "not computed", so check `*_gapless` |
| `a_gapless` / `b_gapless` | 1 | 1 = the whole chain runs first, first+1, … with no gaps or insertion codes. Only with offset 0 **and** gapless 1 do file position, residue number and ProteinMPNN's own position all coincide |

**ProteinMPNN's position ≠ file position on a gapped chain.** ProteinMPNN fills residue-number gaps
when it parses a PDB, while posrmsd counts residues in the file. The two agree only for a gapless
chain, which is what `*_gapless` = 1 certifies. rfd3 renumbers every output chain from 1 without
gaps. A chain extracted with `chainsel --renumber-from` does **not** necessarily satisfy this.

## The question it answers, and the threshold

**"Over these explicitly given positions, how far is A's CA trace from B's after a superposition
on those positions alone?"** For the campaign that commissioned it: did the Boltz prediction of
each ProteinMPNN sequence keep the held-fixed graft residues in the geometry of the backbone they
were designed on?

**No house threshold.** A sub-Å bar is conventional for motif fidelity in benchmarks, but nothing
here has been measured on this repo's data. Pick a bar from the first real distribution and record
it in the campaign log with the run it came from. **Gate on `worst_seg_rmsd`, not `rmsd`**, when the
spec has more than one segment: one loop at 2 Å averaged with three at 0.5 Å reads ~1 Å globally.

## Invocation and flags

```bash
sapia run posrmsd outputs/RUN -t table1 \
    --col-a boltz_path \
    --col-b chainsel_binder_path \
    --positions-column motifrmsd_design_fixed_positions \
    --chain-a B

sapia collect posrmsd outputs/RUN -t table1 \
    --col-a boltz_path --col-b chainsel_binder_path
```

**Not yet run on Modal or vib.** The command above was run end to end only through prosapia's
real run and collect drivers with a stub executor, on synthetic files.

| Flag | Default | Meaning |
| --- | --- | --- |
| `--col-a` | **required** | structure A, the prediction: the structure that is **moved**. It **selects the rows**: a row with no value here is not submitted and collects `missing`. Read from the `-t` table only, not up the lineage |
| `--col-b` | **required** | structure B, the reference/backbone: the **fixed frame**. **Resolved up the lineage** with the same `ctx.lookup` walk usalign uses for `--col-b` (the row's own value, else the nearest ancestor's). A child's own value shadows its parent's |
| `--positions` | `""` | one literal spec for every row, e.g. `'11:18,45:54,87:92'`. **Mutually exclusive with `--positions-column`; exactly one is required** |
| `--positions-column` | `""` | a column holding a per-design spec, e.g. `motifrmsd_design_fixed_positions`. **Resolved up the lineage** |
| `--chain-a` | `""` | the chain of A to index. **Empty: A must hold exactly one polymer chain**, otherwise that row is an error naming the chains. A Boltz complex therefore always needs it |
| `--chain-b` | `""` | the same for B. **An rfd3 binder-scaffolding backbone has 4 chains (A, B, C target + D binder), so it needs `--chain-b D`**, or chainsel the binder out first |
| `--output-prefix` | `<col_a minus _path>_vs_<col_b minus _path>` | names the results subdir and every column. **collect must be given the same value**, or the same `--col-a`/`--col-b`. **Deviation from usalign**: the default strips a trailing `_path` suffix only (`str.removesuffix`), where usalign's `.replace('_path', '')` removes every occurrence anywhere in the string. The two differ only when `_path` appears mid-name (e.g. `boltz_path_relaxed` → usalign's default prefix is `boltz_relaxed`, posrmsd's is `boltz_path_relaxed`); porting a usalign-style command in that case gives a different default prefix, not a silently wrong one — collect fails loudly with `Results dir not found` if the two disagree |
| `-l/--dir-label` | `""` | standard; it adds `_<label>` to the leaf (`posrmsd_<label>_<prefix>_…`) and to the out dir |
| `-i/--input-column` | `not applicable` | **not used, and refused if given**: rows are selected by `--col-a` |

Refused **at submit, with nothing queued**. As an `update` tool it registers no table, so a
refusal leaves no phantom registry row, only the out dir with its `.meta.json` and an empty
`posrmsd_logs/`. These are refused:

- both or neither of `--positions` / `--positions-column`;
- `-i`;
- a chain flag that is not a single id;
- a `--col-a` that is not a column of `-t`;
- a `--col-b` or `--positions-column` that resolves on **none** of the selected rows (a typo or
  the wrong table);
- a prefix with characters other than letters, digits and `_`;
- a brace-free literal `--positions` that is malformed.

Per-row problems are written as status rows and the rest of the batch is submitted. Those cover a
column spec that has a `/`, an open-ended range, a fragment such as `18:11` that expands to
nothing, a position below 1, a position listed twice, fewer than 3 positions, or an unresolvable
`{expr}`.

Resources: 2 CPU, 4 GB, 20 min per task; the fit takes milliseconds.
**(Modal only)** The image is `debian_slim + gemmi + numpy`; it has no prosapia, because the spec
reaches the task already parsed, as JSON.
**(vib)** There is **no `SAPIA_ACTIVATE_POSRMSD` in `.env` yet**. Every vib task will die at
activation until one is added. The sibling tools use
`SAPIA_ACTIVATE_POSRMSD="${SAPIA_VIB_WORKSPACE}/activation/vib/prosapia_python.sh"`, the same
script as motifrmsd.

## Prefix contract (same as usalign)

All columns are prefix-keyed and the status is **per comparison**:
`posrmsd_<prefix>_status`, not `posrmsd_status`. One out dir can therefore hold several comparisons
on one table. Files are `<run_dir>/<table>/posrmsd[_<label>]/<prefix>/<name>.tsv` and
`<name>_sup.cif`. Collect re-reads every row of the table on every call, so it is idempotent.

**The prefix is derived only from `--col-a`/`--col-b` (or an explicit `--output-prefix`) — never
from the spec or the chain flags.** A rerun under the *same* prefix after changing
`--positions`/`--positions-column`, `--chain-a` or `--chain-b`, without `--force`, does **not**
"overwrite the first" — it skips every row already `OK` under that prefix (prosapia's normal
rerun semantics) and only the not-yet-`OK` rows get the new spec/chains. The result is a **mixed
column**: some rows measured against the old spec, some against the new, all under one name with
no flag anywhere distinguishing them (only inspecting `seg_ranges`/`a_chain`/`b_chain` per row would
reveal it). **Always pass `--force`, or use a new `--output-prefix`, whenever the spec or the chain
flags change** — never rerun bare.

## Gotchas, each with its signature

- **Every row `error: … has 4 polymer chains (A, B, C, D) and --chain-b was not given`**: B is the
  full rfd3 complex. Pass `--chain-b D`, or compare against a chainsel-extracted binder.
- **`seq_match_frac` 1.0 on most rows and low on a few, with `a_motif_seq` scattered single
  differences**: real held-position violations by ProteinMPNN on those rows. **Low on every row with
  segments shifted by one letter**: the spec is off by one against the file order (residue numbers
  used as positions?). Check `b_resnum_offset`.
- **`seq_match_frac` 0.0 with `b_motif_seq` all `X`**: B carries UNK residues at those positions.
  Identity cannot be confirmed; that is not a mismatch.
- **`missing_input: --col-b 'x' has no value on this row or any ancestor`** on some rows only: those
  rows' lineage does not reach a table that has the column (a different parent, an orphan row). On
  **all** rows the submit is refused instead, with
  `ValueError: --col-b 'x' has no value on any of the N selected row(s)`.
- **`Submitting N designs` smaller than the table**: rows without `--col-a`, rows already `OK`
  under this prefix (rerun with `--force`), plus the `missing_input:` / `error: position spec` rows
  printed just above that line. Each of those has a status row on disk.
- **`FileNotFoundError: Results dir not found: …/posrmsd/<prefix>` at collect**: collect was not given
  the same `--col-a`/`--col-b` or `--output-prefix` (or `-l`) as the run.
- **`warn: chain lengths differ` on every row**: A and B are not the same construct, for example a
  full complex chain against a trimmed chain. The numbers exist, but the same-ordering premise is
  yours to vouch for.

## What it collects

All columns are `posrmsd[_<label>]_<prefix>_<name>`. Integer columns read back as floats when any row
is NA (pandas), so `== 1` and `== 1.0` both work.

| `<name>` | Kind | How to read it |
| --- | --- | --- |
| `status` | trust | see the trust table; **`OK` is the only success** |
| `a_chain`, `b_chain` | trust | the chains indexed |
| `a_len`, `b_len` | trust | polymer lengths; they must be equal for `OK` |
| `n_positions` | trust | positions compared |
| `seq_match_frac` | trust | **read first**; ~1.0 expected in the held-fixed use |
| `a_motif_seq`, `b_motif_seq` | trust | residues read, segments comma-joined |
| `a_resnum_offset`, `b_resnum_offset` | trust | constant resnum − position, or NA if not constant |
| `a_gapless`, `b_gapless` | trust | 1 = file position equals ProteinMPNN position |
| `rmsd` | science | CA RMSD (Å) after the one Kabsch fit of A onto B over all positions |
| `max_dev` | science | largest single-CA deviation (Å) after that fit; one bad residue hides in `rmsd` but not here |
| `n_segments` | science | number of comma-separated fragments in the spec. **The fragments ARE the segments**: never merged or re-split |
| `seg_ranges` | science | the fragments as resolved, `/`-joined: `11:18/45:54/87:92` |
| `seg_rmsds` | science | their RMSDs, `/`-joined, same order |
| `seg0_rmsd`, `seg1_rmsd`, … | science | one column per segment, **from the same global fit** (never refit per segment). Present only for indices some row has; NA'd on re-collect if a row no longer has them |
| `worst_seg_rmsd` | science | the max over segments; **the filter-safe scalar**, on every OK/warn row whatever the segment count |
| `sup_path` | output | structure A moved into B's frame (mmCIF), for viewing |

Every `error:`, `missing_input:` and `missing` row carries **NA in every metric column, never 0**.
An RMSD of 0 would read as a perfect match.

## Filter

```python
P = "posrmsd_boltz_vs_chainsel_binder"   # posrmsd[_<label>]_<prefix>

def apply_filter(df):
    return df[
        (df[f"{P}_status"] == "OK")
        & (df[f"{P}_seq_match_frac"] == 1.0)      # trust: positions really held, indexing right
        & (df[f"{P}_worst_seg_rmsd"] < 1.0)       # <- campaign bar, NOT a house threshold
    ]
```

Rows with NA in any of these fail every comparison, so they are dropped. That is correct, because NA
is the absence of a measurement.

## What it is blind to

- **Everything outside the listed positions.** It cannot tell whether the rest of the chain folded
  as designed. Use `usalign --mm 0` on chainsel-extracted chains for that.
- **The pose relative to a target.** The fit is on the positions alone, so a motif can be perfect
  while the binder sits in the wrong place. Use `usalign --mm 1`, `epitope`, or `dimerfit`/`ringfit`
  depending on the target.
- **Side chains.** It measures CA only.
- **Whether the positions are the ones you meant.** The spec is taken as given, and a correct RMSD
  over the wrong residues is still the wrong answer. Check the motif sequences against the motif you
  grafted.
- **Prediction confidence.** A confident wrong prediction and a confident right one are told apart
  here only through geometry. Read the predictor's pLDDT/PAE separately.

## Traps found while building it

- **Two parsers would drift.** A second implementation of the position grammar could disagree with
  ProteinMPNN's about whitespace, `[...]`, empty fragments or `{expr}`. The builder therefore
  parses each fragment with prosapia's `parse_chain_positions` and asserts the concatenation equals
  `parse_positions(spec)`. A disagreement raises
  `posrmsd's segment split … disagrees with prosapia.parse_positions`.
- **`18:11` is silently empty in prosapia.** `parse_chain_positions` returns no positions for a
  reversed range and raises nothing. Here that is refused per row:
  `fragment '18:11' … expands to no positions`. In ProteinMPNN the same string fixes nothing and
  exits 0.
- **Overlapping fragments** (`1:5,4:8`) are legal in prosapia (it does not de-duplicate) but would put
  one CA into the fit twice: `lists position(s) 4, 5 more than once`.
- **A chain letter that reappears after other chains in a PDB** (its waters listed at the end) is
  one chain to the worker. gemmi 0.7.5 already merges such chain parts on read (checked), and the
  worker calls `merge_chain_parts()` again as an idempotent guard. Either way the waters are
  non-polymer, so they never count as positions.

## Status

**Verified on real data, twice, on Modal.** Built 2026-10-10. Its 178 self-test checks pass (`uv run
python tools/posrmsd/test_posrmsd.py`). `tool-reviewer` (2026-10-10) found no code defect that
produces a believable wrong number, and flagged two things fixed above: (1) `seq_match_frac == 1.0`
does not by itself certify the spec is the graft (see the trust table); (2) a rerun under one prefix
after changing the spec/chains silently mixes old and new measurements under one column unless
`--force` or a new `--output-prefix` is used.

Two production runs on `outputs/20261009_130629_5M2J_vhh_redmotif_lowtemp_scaleup`, 2704 rows each,
both `status OK` on every row:

- **2026-10-10**, `table1` (ProteinMPNN, 24 pinned positions): `rmsd` median 2.078 Å.
- **2026-10-14**, `table1_red17` (atomium, 17 pinned positions): `rmsd` median 2.007 Å,
  `seq_match_frac == 1.0` on all 2704, `b_motif_seq` equal to the expected motif identity on
  2704/2704, and segment RMSDs reconstructing the global `rmsd` to within 0.0007 Å.

Both runs reproduce the headline caveat this tool exists to expose: **a motif-restricted RMSD can be
worse than the whole-chain self-consistency RMSD even when sequence identity at the pinned positions
is perfect on every row** (2.0-2.1 Å motif vs ~1.6 Å whole-chain in both arms). Sequence held and
geometry held are independent claims.

**Comparing `worst_seg_rmsd` across runs is only valid when the specs have the same segmentation.**
A spec of the same residues with interior positions removed fragments into more, shorter segments —
24 pins as 3 contiguous segments became 17 pins as 9 segments including 5 singletons — and after one
global fit a singleton's segment RMSD is just that residue's displacement. The worst of 9 short
segments is drawn from a noisier distribution than the worst of 3 long ones, so the metric changes
meaning between specs even when the geometry does not.
