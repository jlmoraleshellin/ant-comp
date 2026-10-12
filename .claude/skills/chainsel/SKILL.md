---
name: chainsel
description: >-
  Rewrites a structure column's chains into a new file (`chainsel_path`): extract the binder, merge split protomers, or cut a chain by residue position. Load before a binder self-consistency comparison, ringfit on a split prediction, or feeding an RFdiffusion symmetric design to epitope, atomium or usalign. Trap: a wrong chain or cut does not error, so read `chain_n_res` and `chain_mean_b`.
---

# chainsel

**Custom tool** (lives in `tools/chainsel/`, not bundled with prosapia).
**`action: update`** — a sub-structure is a *property* of a design that already exists,
so it annotates the table it reads, in place, and mints no child. `-t` is required.
**CPU-only**; the manifest builder forces `gpus_per_task = 0`, so you do not pass `-g 0`.

Extracts the chains you name from each design's structure into a new file — optionally
**merging** several source chains into one output chain, or **splitting** one source
chain into several by residue position — and collects its path as
`chainsel_path`, an ordinary structure column that `usalign`, `pyrosetta`, `ringfit` or
a viewer then consumes like any other.

## Why it exists

Every structure comparison in this workbench is **whole-file**, and a binder complex is
mostly target:

- **`usalign --mm 1`** (the default) aligns the whole multi-chain file. On a
  target+binder complex the score is dominated by the large, fixed target and reads
  ~0.95–1.0 **whatever the binder did** — including when the binder is a completely
  different fold.
- **`usalign --ter 2`** isolates only the **first** chain. In every complex we build
  (`mkcomplex` puts the fixed target chains first) the binder is the **last** chain, so
  `--ter 2` compares target to target.

So before chainsel there was **no way to compare a binder to a binder**. That matters
because a binder campaign has two separate questions, and one number cannot answer both:

| Question | What to compare | How |
| --- | --- | --- |
| **Did it fold?** | binder chain alone, prediction vs. design | `chainsel` both, then `usalign --mm 0` |
| **Did it stay?** | the pose in the target's frame | whole complex: `usalign --mm 1`, or `ringfit` |

A binder can fold perfectly and sit in the wrong place, or hold its pose in a prediction
that folded into a different topology. Measure them separately.

### The second reason: split protomers

A target trimmed of a membrane belt is **two chains per protomer**. SlyB here is
residues 19–59 (41 aa) and 107–155 (49 aa), with a real 13.4 Å gap between the
segments, so a predictor is handed two entities per protomer and returns a complex
where protomer 1 is chains `A` + `C` and protomer 2 is chains `B` + `D`.

`ringfit` takes **exactly two** `--target-chains` — one per protomer — and there is no
syntax for "t1 is the union of A and C". On a split-protomer file the per-protomer
metrics `bsa_t1`, `bsa_t2` and `bridge_ratio` (does the binder straddle the seam?)
cannot be computed *at all*. `--merge-groups "A+C:A,B+D:B"` writes the same coordinates
back as one chain per protomer, and `ringfit --target-chains A,B` works again.

### The third reason: fused subunits (split mode)

RFdiffusion's symmetric mode writes **one chain per subunit that fuses the designed
binder and the fixed target**: chains A, B, C of 212 residues each, positions 1–60 the
diffused binder (GLY, B-factor 0) and positions 61–212 the fixed target (B-factor 1),
with residue **numbers** running 1–60 and then 1–152 — so 60 (chain, number) pairs
repeat inside every chain. `epitope`, `atomium` and `usalign` need binder and target as
separate chains with unique numbering. `--split-groups` cuts them apart. See
[Splitting a fused chain](#splitting-a-fused-chain---split-groups).

## Invocation

```bash
sapia run chainsel <run_dir> \
    --table table1 \
    --input-column boltz_path \
    --chains E --rename-to A \
    --dir-label pred_binder
sapia collect chainsel <run_dir> --table table1 --dir-label pred_binder
```

| Flag | Default | Notes |
| --- | --- | --- |
| `-i/--input-column` | **effectively required** | `default_input_column` is the sentinel `"not applicable"` (the `usalign` convention). No column is a defensible default — the tool is equally at home on `rfdiffusion3_path`, `boltz_path`, `alphafold3_path` — so the builder **raises** at submit time unless you pass a column that the table actually has. |
| `--chains` | `""` | Chain IDs to KEEP, comma-joined — **one output chain each**. **The output order follows the order given**, so `--chains E,A` writes E first. **Required unless `--merge-groups` or `--split-groups` is given.** |
| `--rename-to` | `""` | New IDs, same count and order as `--chains`. Empty keeps the originals. Use it to normalise a binder to one fixed letter across tables whose lettering differs. Two kept chains may **not** be renamed onto the same ID — that guard stays a hard error, because there it would be an accidental merge. |
| `--merge-groups` | `""` | **Merge several source chains into one output chain.** `"<chain>[+<chain>...]:<out_id>"`, groups comma-joined: `"A+C:A,B+D:B"`. Replaces `--chains`/`--rename-to` (it states both). The only way to get two source chains into one output chain. |
| `--out-format` | `pdb` | `pdb` or `cif`. `pdb` is what USalign and Rosetta want. With `pdb`, an output chain ID longer than one character is refused at submit time (PDB has a single chain-ID column). |
| `--keep-het` | off | Default **strips** HETATM — ligands, ions, waters. `--keep-het` keeps them. Modified polymer residues (MSE and friends) are kept either way: the test is "tabulated amino/nucleic acid", not the HETATM flag. |
| `--split-groups` | `""` | **Cut source chains into several output chains by 1-based POSITION.** `"<chain>:<first>-<last>><out_id>[+<first>-<last>><out_id>...]"`, source chains comma-joined: `"A:1-60>D+61-212>A,B:1-60>E+61-212>B,C:1-60>F+61-212>C"`. Replaces `--chains`/`--rename-to`/`--merge-groups`/`--renumber-from` (any combination is a submit-time error). **Always renumbers** each output chain 1..N. Accepts `{expr}` placeholders. See the split section. |
| `--renumber` | off | Renumber each **output** chain 1..N in written order, insertion codes dropped — **continuous across a merge seam** (`A+C:A` with 41+49 residues comes out 1..90). Off preserves the source numbering. **A redundant no-op with `--split-groups`**, which always renumbers. |
| `--renumber-from` | `""` | Restore auth numbering **per source chain**: `"<chain>:<first>[-<last>]"`, comma-joined — `"A:19-59,C:107-155"`. Mutually exclusive with `--renumber`. **An error with `--split-groups`.** |
| `-l/--dir-label` | `""` | **In practice required** — see the traps. |

Reads PDB or mmCIF transparently (gemmi). Alternative conformations and hydrogens are
dropped on read; only the **first model** is used. The non-split modes **refuse** a kept
chain whose residue numbering repeats inside it (see the traps) — use `--split-groups`
for those.

## Merging chains (`--merge-groups`)

```
--merge-groups "A+C:A,B+D:B"
                └┬┘ └┬┘
   source chains ─┘   └─ the output chain ID
```

Merging is **never** implicit. `--rename-to` still refuses to map two kept chains onto
one ID; only `--merge-groups` merges, and it says so in the flag.

**Residue order inside a merged chain is the order the group writes its sources** —
all of A's residues, then all of C's — and is **never re-sorted by residue number**. A
predictor numbers each chain 1..N independently, so sorting by number would interleave
the two segments catastrophically.

### The numbering contract — read this before every merge

Exactly one of three regimes applies, and the output always carries one of them:

| Flags | What comes out |
| --- | --- |
| neither (default) | Every residue keeps its **source** number and insertion code. A merged chain whose segments both start at 1 therefore **collides** → that design is an **error**, not a file with repeated numbers. |
| `--renumber` | Each **output** chain is numbered **1..N in written order**, insertion codes dropped. `A+C:A` with a 41- and a 49-residue segment → **1..90 continuous across the seam** (not 1..41 then 1..49). |
| `--renumber-from "A:19-59,C:107-155"` | Each **source** chain's segment is numbered **consecutively from its own first number**, insertion codes dropped. Output chain A carries **19..59 then 107..155** — the target's real auth numbering, which is what `ringfit --resnum-match resnum` pairs against the reference. |

One sentence, for the record: **`--merge-groups` writes each group's sources in the
order given, and the residue numbers in the output are the source numbers unchanged
(default, an error if they collide inside a merged chain), or 1..N per output chain
with `--renumber`, or `<first>`, `<first>+1`, … per source chain with
`--renumber-from`.**

Two things `--renumber-from` does *not* do: it does not reproduce numbering **gaps
inside** a segment (it writes `first`, `first+1`, … in file order), and it does not
guess — every kept source chain must be listed. Give the `-<last>` whenever you know
it: the residue count actually written is checked against the range, so a segment of
the wrong length is an error rather than a silently shifted numbering.

### Guards, all at submit time (the run refuses; nothing is queued)

| Combination | Result |
| --- | --- |
| none of `--chains`, `--merge-groups`, `--split-groups` | error — chainsel needs a selection |
| `--chains` **and** `--merge-groups` | error — they are alternatives |
| `--rename-to` **and** `--merge-groups` | error — a group's output ID is the part after `:` |
| `--renumber` **and** `--renumber-from` | error — both set the numbering, pick one |
| a source chain in two groups, or twice in one | error — each source chain is written once |
| two groups with the same output ID | error — write them as one group (`A+C:A`) |
| `--rename-to` mapping two chains onto one ID | error (unchanged) — say `--merge-groups` if you mean it |
| a malformed group (`A+C`, `A+C:`, `A:B:C`) | error naming the entry |
| `--renumber-from` naming a chain this run does not keep | error naming it |
| `--renumber-from` missing a kept source chain | error — a half-renumbered file mixes two regimes |
| `--renumber-from` with `end < start` or a non-integer | error naming the entry |
| `--out-format pdb` with an output ID longer than 1 char | error (unchanged) |

Per design, recorded as data in `chainsel_status` (not a crash, no file written):

- **a chain named in a group but absent** → `error: chain C not in <file> (have: A,B)`;
- a segment length that contradicts its `--renumber-from` range →
  `error: chain C of <file> has 48 residues but --renumber-from says 107-155 (49)`;
- a merged chain whose residue numbers collide →
  `error: merged chain A of <file> (A+C) has duplicate residue numbers (1,2,...)`.

## Splitting a fused chain (`--split-groups`)

```
--split-groups "A:1-60>D+61-212>A,B:1-60>E+61-212>B,C:1-60>F+61-212>C"
                │ └┬─┘ │ └──┬─┘ │
                │  │   │    │   └── output ID of segment 2
                │  │   │    └────── positions 61-212
                │  │   └─────────── output ID of segment 1
                │  └─────────────── positions 1-60
                └────────────────── source chain
```

Comma separates **source chains**; `+` separates one chain's **segments**; each segment
is `<first>-<last>><out_id>`.

**Premise.** Each named source chain is a concatenation of pieces you want as separate
chains, and you know where each piece starts and ends **by position**. `<first>` and
`<last>` are **1-based positions in the chain's residue order in the file, NOT residue
numbers** — in a fused RFdiffusion chain the numbers repeat (binder 1–60, then target
1–152), so a number cannot address a residue. Positions count the residues chainsel
keeps: polymer residues only by default (waters/ligands are stripped first), every
residue with `--keep-het`.

A position is a residue **as delimited in the file text**: a new one starts at any
change of chain/number/insertion code/name, or when an atom name + altloc repeats.
Hydrogens and alternative conformations are then dropped **within** each residue; no
residue is ever removed by that. Two consequences:

- **Microheterogeneity counts as two positions.** Two alternative residues of different
  names at one site (ALA altloc A + SER altloc B, same number) are two positions, and
  **both are written** (each keeps its first conformer). The checksum will usually flag
  it as one residue too many. Design outputs do not contain this; deposited PDB
  entries can.
- **A residue made only of hydrogens is an error**, not a position:
  `error: chain A of <file> has residue(s) with no atoms left after dropping
  hydrogens/alternative conformations, at position(s) 30`. (gemmi 0.7.5's
  `remove_hydrogens` leaves such a residue in place with zero atoms, which would count
  in `chain_n_res` and write nothing.)

**The numbering contract, in one line: split mode ALWAYS renumbers — every output chain
is numbered 1..N in written order, insertion codes dropped.** The source numbering
repeats, so preserving it is not an option. `--renumber` is accepted as a redundant
no-op; `--renumber-from` is an error. On the RFdiffusion layout the target chains come
out 1..152 — the same numbers the target carried in the source — and the binder chains
1..60. That is RFdiffusion's numbering, not the original target PDB's: if the target
was renumbered by the generator, hotspots in the original numbering still need the
offset (see the epitope skill).

**Output order is ascending output ID, whatever the listing order**, residues in source
order within each chain. The spec above writes A, B, C (targets) then D, E, F (binders).
Listing `C:...` first changes nothing. With `--out-format cif` and multi-character IDs
the order is plain string sort.

**The coverage is a checksum.** For each source chain the segments must **tile it
exactly**: start at position 1, run on contiguously in the order listed, no gap, no
overlap, and **end at the chain's real residue count**. A chain not named is
**dropped** (as with `--merge-groups`); there is no syntax to drop *part* of a chain,
on purpose — a trimmed fit is exactly the silent failure this tool exists to prevent.

### Trust columns — read these before using the file

| Column | Expected on the RFdiffusion layout | What a different value means |
| --- | --- | --- |
| `chain_mean_b` | `A=1.00,B=1.00,C=1.00,D=0.00,E=0.00,F=0.00` | Mean B over each output chain's written atoms, 2 decimals. RFdiffusion writes B 0 on diffused (binder) residues and B 1 on fixed (target) ones. A value strictly between 0.00 and 1.00 **does** mean a mixed chain — but the converse does **not** hold: two decimals cannot show one leaked residue on a chain of about 200 residues or more (152/153 = 0.99 shows; 200/201 = 1.00 does not), and the leak shows only on the chain that *received* the residue. **Never read it without `chain_n_res`.** |
| `chain_n_res` | `A=152,B=152,C=152,D=60,E=60,F=60` | Residues written per output chain. A count you did not design means the cut is in the wrong place. |
| `n_chains` / `chains` | `6` / `A,B,C,D,E,F` | Fewer chains: a source chain was not named and was dropped. |
| `n_res` | `636` | Total; equals the sum of `chain_n_res`. |

Both `chain_*` columns are **NA in every other mode** (`--chains`, `--merge-groups`) —
never 0.

### Guards (submit time — the run refuses, nothing is queued)

| Combination | Result |
| --- | --- |
| `--split-groups` with `--chains`, `--rename-to`, `--merge-groups` or `--renumber-from` | error naming every clashing flag |
| `--split-groups` with `--renumber` | accepted; a no-op |
| a gap (`1-60>D+62-212>A`), a start other than 1, or segments out of order | error: `gap -- segment ... leaves position N uncovered` |
| an overlap (`1-60>D+60-212>A`) | error: `overlaps the previous one` |
| a segment with no `>` (`A:1-60>D+61-212`), `A:60-1>D`, `A:1-60>`, `A:0-60>D`, `A:60>D`, a non-integer, a trailing `+` | error naming the segment |
| a source chain named twice | error — give all its segments in one entry |
| the same output ID twice anywhere in the spec | error — split mode never merges |
| `--out-format pdb` with an output ID longer than 1 char | error (`--out-format cif` allows it) |

With `{expr}` placeholders (`A:1-{binder_len}>D+{binder_len + 1}-212>A`) the spec is
resolved and checked **per design** at submit time; a design that fails aborts the
whole submit with `design '<name>': ...`.

### Per-design errors (recorded as data in `chainsel_status`, no file, every column NA)

- the segments do not end at the chain's real residue count →
  `error: chain A of <file> has 212 residues but --split-groups 'A:1-60>D+61-211>A' covers 211`;
- a named source chain is absent → `error: chain Z not in <file> (have: A,B,C)`;
- a residue left with no atoms after dropping hydrogens (see the premise above);
- the reader returned a different set or order of residues than chainsel delimited →
  `error: chain A of <file>: the reader returned ... -- refusing to cut by position`
  (should never happen; it is the guard on the position bookkeeping below).

### Why split mode reads the file differently

Measured with gemmi 0.7.5, on exactly this layout: gemmi's PDB **and** mmCIF readers
merge a residue into an **earlier** residue of the same chain with the same number,
insertion code and name — adjacent or not — so a target GLY numbered 12 is folded into
binder GLY 12 as one residue with 8 atoms; and `remove_alternative_conformations` then
**deletes** every later residue whose number repeats an earlier one. On a fused chain
that silently removes target positions 61–120. Split mode therefore rewrites each
residue's number to its **position** in the file text before gemmi parses it (a new
residue starts at any change of chain/number/insertion code/name, or when an atom
name + altloc repeats inside the current residue — so two back-to-back `GLY 3` records
stay two residues), checks that gemmi returned exactly those positions in order, and
only then drops alt-confs and hydrogens. Atoms, coordinates and B-factors are cloned
unchanged. With an mmCIF output each output chain gets its own `label_asym_id` and
entity (a clone would otherwise keep its source chain's, so binder D and target A cut
from one chain would share `label_asym_id` A — verified, and fixed); identical chains
are **not** merged into one entity.

## Traps

**`--chains` / `--merge-groups` refuse a chain whose numbering repeats inside it**
(since 2026-10-07). These modes use gemmi's ordinary reader, which merges or deletes
such residues: before the guard, `--chains A` on the fused RFdiffusion layout returned
`chainsel_status = OK` with **`n_res = 152` instead of 212**. Now that design comes back
`error: numbering repeats inside chain A of <file>: residue GLY 1 holds the same atom
twice (two residues merged on read); this mode would write a wrong chain -- use
--split-groups` (or `... N residues (ALA/SER) are numbered 1; the reader would drop all
but one -- use --split-groups`). What is checked, exactly:

- only the chains you keep, and only the residues that would be written (polymer only
  unless `--keep-het`);
- "the same atom twice" looks at **heavy atoms of polymer residues** only — a kept
  ligand whose carbons are all named `C` (docking/OpenBabel output) and hydrogens with
  truncated, repeated names are legitimate and pass;
- "numbered N" fires when **two or more** residues at one number have blank-altloc
  atoms, or two of them share an altloc letter. **Microheterogeneity** — alternative
  residues at one site with their own altloc letters, at most one of them with
  blank-altloc atoms (e.g. ARG 12 with a blank backbone and an altloc-A side chain,
  plus LYS 12 with an altloc-B side chain) — is exempt and behaves as before. A
  fused-chain collision always fires, because binder and target residue both have
  blank altlocs — **unless the colliding residue carries altlocs on every atom** (see
  the known gaps below).
- **A stripped residue sharing a number is ignored only when it comes AFTER the
  polymer residue in the file** (the PDB/mmCIF convention: polymer first, waters
  last). If a water numbered 5 comes *before* polymer residue 5, gemmi's reader drops
  the **polymer** residue, the guard does not see it (the water is stripped), and the
  design comes back `OK` one residue short. Signature: `n_res` one below the chain
  length. Not guarded; check `n_res` on unusual inputs.

**Known gaps in the guard** (reviewer, round 3; documented, not fixed — the code under
verification is frozen). Each comes back `OK` and short:

1. **One blank-altloc residue beside a fully altloc'd one at the same number.** A
   residue whose atoms all have a blank altloc, next to another at the same number
   whose every atom has an altloc, is accepted as microheterogeneity; the reader then
   deletes the later residue. Status `OK`, `n_res` one short.
2. **A fully altloc'd GLY merged into a blank-altloc GLY.** No (atom name, altloc) pair
   repeats and there is no group of two residues, so neither check fires. The result is
   believed to carry two N atoms in one residue; **not verified**.
3. **A merge of non-polymer residues under `--keep-het` is no longer refused** (check
   (a) looks at polymer residues only): `n_res` undercounts by one.

All three need altlocs or kept ligands that share a number with something else. That
is implausible for RFdiffusion, rfd3, Boltz and AF3 output, none of which write
altlocs; it matters for deposited PDB entries. The proposed tightening is deferred:
make (b) also fire when any residue in the group has no altloc atom at all, and make (a)
flag a polymer residue whose heavy-atom name appears both with a blank altloc and with
a letter.

On every input the guard accepts, the output is byte-identical to before (standing test
`test_non_split_output_is_byte_identical_to_before`, 18 cases). To keep a whole fused
chain, use one segment: `--split-groups 'A:1-212>A'`. **Rows collected `OK` before
2026-10-07** from a source with repeating numbering may be short. Re-run them with
`--force` on both `run` and `collect`, or under a new `-l`: a plain re-run skips rows
that are already `OK` (see the next trap).

**A corrected spec does not replace an old cut by itself.** Reruns skip rows that are
already `chainsel[_<label>]_status == OK`, so re-running with a fixed `--split-groups`
under the same `-l` leaves those rows with the **old** file and the **old**
`chain_n_res` / `chain_mean_b` strings. Pass `--force` to `run` and `collect`, or use a
new `-l` (new columns, the old ones kept for the record).

**A wrong cut that still tiles is not caught.** `A:1-59>D+60-212>A` covers all 212
residues, so it succeeds: one binder residue lands in the target chain. Nothing errors.
Only the trust columns show it — `chain_n_res` reads `A=153,D=59` and `chain_mean_b`
reads `A=0.99` instead of `1.00`. `chain_n_res` is the one that always shows a wrong
cut whose lengths you know; `chain_mean_b` misses one leaked residue on a chain of
~200+ residues. Read both on every row before using the split file; the filter below
does it.

**An absent chain is an error, never a shorter file** — including a chain named inside
a merge group. This is the whole point of the
tool's error discipline. If `--chains E` and the structure has no chain E, that design
comes back `chainsel_status = "error: chain E not in <file> (have: A,B,C)"`, with an
empty path and NA counts, and **no file is written**. A silently partial extraction
would be a perfectly well-formed PDB that every downstream TM-score, RMSD, BSA and
energy would then be computed on — wrong, plausible, undetectable.

**But a chain that is present and WRONG is not caught.** Extracting chain `A` when the
binder is chain `E` succeeds: you get a clean file of the target. Nothing errors.
Guard it by reading `chainsel_n_res` against the binder length you designed, and
`chainsel_chains` against what you asked for. If `n_res` is the target's length, you
extracted the target.

**Chain lettering differs between tools.** rfdiffusion3 and a structure predictor need
not agree on which letter the binder gets: the predictor assigns letters in the order of
the `/`-separated entities `mkcomplex` built, so the binder is the **last** letter there,
while the designed backbone may call it something else. Look at an actual file before
choosing `--chains` — do not assume. `--rename-to` then normalises the two extractions
to the same letter.

**`-l/--dir-label` is effectively required.** The tool leaf (`chainsel`) names both the
output dir and every column, so a second chainsel run on the same table with a different
chain set would overwrite the first's `chainsel_path`. Label each extraction
(`-l design_binder`, `-l pred_binder`) and the columns become
`chainsel_design_binder_*` / `chainsel_pred_binder_*`. **The same `-l` must be passed to
`collect`** — it resolves the same output dir.

**Failures still exit 0.** Errors are recorded as data, so `n_tasks` `.exit` files all
`0` does not mean the step succeeded. The truth is the `chainsel_status` column after
collect.

**Rows can be dropped from the manifest silently.** A design whose input cell points at
a file that does not exist is printed as `MISSING ... (skipping)` and left out.
`Submitting N designs` under the table's row count means exactly that — compare them.

**`--renumber` changes what "residue 12" means.** USalign aligns by sequence, so it does
not care; anything that pairs residues **by number** across two files does
(`ringfit --resnum-match resnum`, a hotspot list, a mutation table). Renumber only when
the consumer needs matching numbering, and never when you will look residues up by their
original numbers. On a merged chain `--renumber` runs **across the seam**, so the second
segment's residue 1 becomes residue 42 — use `--renumber-from` instead whenever the
downstream tool matches by number against a reference.

**A merge lowers `n_chains` but not `n_res`.** `--merge-groups "A+C:A,B+D:B"` on a
4-chain prediction collects `chainsel_n_chains = 2` and `chainsel_chains = "A,B"`,
while `chainsel_n_res` still counts every residue written — 41 + 49 = **90 on chain
A**. If `n_res` dropped, a segment went missing; if `n_chains` did not drop, the merge
did not happen.

## Columns collected (`chainsel_` prefix, or `chainsel_<label>_` with `-l`)

| Column | Meaning |
| --- | --- |
| `_path` | **The extracted structure.** The point of the tool — feed it to the next tool as `--col-a` / `--col-b` / `-i`. |
| `_chains` | The chain IDs **actually written**, comma-joined, in output order (AFTER `--rename-to` or a group's `:<out_id>`). This is what you pass to `ringfit --target-chains`. |
| `_n_chains` | How many **output** chains were written = the number of `--chains` IDs, or the number of **groups**. A merge makes it smaller than the source chain count. |
| `_n_res` | Residues written, summed over the output chains. A merge does **not** change it — a merged chain contributes all its segments (41 + 49 = 90). Your sanity check against the expected length. |
| `_n_atoms` | Heavy atoms written, summed the same way (hydrogens and alt-confs are already gone). |
| `_chain_n_res` | **Trust column, split mode only (NA otherwise).** Residues written per output chain, in output-ID order: `A=152,B=152,C=152,D=60,E=60,F=60`. A string; compare it whole or parse it (below). |
| `_chain_mean_b` | **Trust column, split mode only (NA otherwise).** Mean B-factor over each output chain's written atoms, 2 decimals, same format: `A=1.00,...,D=0.00`. On RFdiffusion output it says which chain is binder (0.00) and which target (1.00); a value strictly between means a mixed chain, but a clean-looking value does not prove a clean cut (one residue on a ~200+ residue chain does not move the second decimal) — read with `_chain_n_res`. |
| `_status` | `OK`, `error: <reason>`, or `missing`. The only success signal. |

**NA (empty) means not applicable, not zero** — the extraction failed, or no TSV was
written. A blank `_n_chains` is not "zero chains", and a blank `_path` is not "nothing
to extract".

On disk: `<run_dir>/<table>/chainsel[_<label>]/<name>.{pdb,cif}` plus `<name>.tsv`.

**A chainsel mmCIF (`--out-format cif`) is not a Boltz forced template.** It has no
`_entity_poly_seq`, in every mode. From a **PDB** source `label_seq_id` is `.` on every
atom (verified, gemmi 0.7.5); from an mmCIF source it is the source's `label_seq_id`
(`--chains`/`--merge-groups`) or 1..N per chain (`--split-groups`). Boltz templates
need `_entity_poly_seq` and a `label_seq_id` running 1..N with no `.` (see the boltz
skill), so build them with `utils/boltz_template.py`, not from a chainsel file.
The per-design TSV carries `chain_n_res`/`chain_mean_b` only in split mode; in the other
modes it is exactly the seven columns it always had.

## Filter: gate on a clean split before epitope / atomium

```python
def apply_filter(df):
    # -l split  ->  chainsel_split_*. NA rows compare False, so errors drop out.
    return df[
        (df["chainsel_split_status"] == "OK")
        & (df["chainsel_split_chain_n_res"] == "A=152,B=152,C=152,D=60,E=60,F=60")
        & (df["chainsel_split_chain_mean_b"]
           == "A=1.00,B=1.00,C=1.00,D=0.00,E=0.00,F=0.00")
    ]
```

The expected strings are a property of the campaign (binder and target lengths), not
of the table: record them in the campaign log. When binder lengths vary, parse instead:
`dict(kv.split("=") for kv in s.split(","))` gives `{"D": "60", ...}`.

## Worked example: split RFdiffusion symmetric designs into binder and target chains

```bash
sapia run     chainsel <run_dir> -t table0 -i rfdiffusion_path \
                  --split-groups 'A:1-60>D+61-212>A,B:1-60>E+61-212>B,C:1-60>F+61-212>C' \
                  -l split
sapia collect chainsel <run_dir> -t table0 -l split

# then, on the split file (chains A,B,C = target, D,E,F = binder):
#   epitope --design-chains D,E,F --target-chains A,B,C -i chainsel_split_path
#   atomium --chains-to-design D,E,F --symmetry 3 -i chainsel_split_path
```

The input column is whatever holds the RFdiffusion PDBs on that table — **check it**;
`-i` has no default. Check one design's `chainsel_split_chain_n_res` and
`chainsel_split_chain_mean_b` before running anything on the split files.

## What it is blind to

- **Whether the cut is the one you meant**, beyond what the trust columns show: on
  a source without a B-factor convention (a prediction, a relaxed model) `chain_mean_b`
  says nothing, and only `chain_n_res` against the designed lengths guards the cut.
- **Geometry.** It never checks that the pieces are folded, in contact or in the right
  frame — that is `usalign` (fold), `epitope`/`ringfit` (pose).
- **Chain identity by content.** It trusts the letters and positions you give it; a
  present but wrong chain is extracted happily (see the traps).

## Worked example: did the binder fold?

Setup: `table0` holds rfd3 designed complexes (`rfdiffusion3_path`, target chains `A,B`
+ binder chain `C`); `table1` holds the ProteinMPNN sequences, their `mkcomplex`
complex sequence, and the Boltz prediction of that complex (`boltz_path`, chains `A,B`
target + `C` binder — **check yours**).

```bash
# 1. binder chain out of the DESIGNED complex — on table0, where rfdiffusion3_path lives
sapia run     chainsel <run_dir> -t table0 -i rfdiffusion3_path \
                  --chains C --rename-to A -l design_binder
sapia collect chainsel <run_dir> -t table0 -l design_binder

# 2. binder chain out of the PREDICTED complex — on table1, where boltz_path lives
sapia run     chainsel <run_dir> -t table1 -i boltz_path \
                  --chains C --rename-to A -l pred_binder
sapia collect chainsel <run_dir> -t table1 -l pred_binder

# 3. binder vs binder, monomer mode
sapia run     usalign <run_dir> -t table1 \
                  --col-a chainsel_pred_binder_path \
                  --col-b chainsel_design_binder_path \
                  --mm 0 --output-prefix binder_fold
sapia collect usalign <run_dir> -t table1 \
                  --col-a chainsel_pred_binder_path \
                  --col-b chainsel_design_binder_path \
                  --output-prefix binder_fold
```

Three things make that work:

- **chainsel runs on the table that holds the column.** It reads `ctx.df` — the table it
  is updating — and does **not** walk the lineage. `rfdiffusion3_path` is on `table0`, so
  the first extraction targets `table0`.
- **usalign's `--col-b` *does* walk the lineage.** From a `table1` row it finds
  `chainsel_design_binder_path` on the parent `table0` row. Nothing is copied forward.
- **`--mm 0` is monomer mode**, right for a one-chain file. `--mm 1` (usalign's default)
  would be wrong here in the opposite direction: it expects multiple chains.

Read the result as `usalign_binder_fold_TM1` / `_RMSD`: RMSD < 2 Å with TM > 0.9 means
the sequence really folds to the binder you designed. **Run the whole-complex comparison
too** (`usalign --mm 1` on `boltz_path` vs `rfdiffusion3_path`, or `ringfit`) and report
the pair — a good binder fold with a bad pose is a binder that folded and then went
somewhere else, and it is a failed design either way.

Verify before scaling: extract one design, check `chainsel_n_res` equals the binder
length and `chainsel_chains` is what you asked for, and open `chainsel_path` in a viewer.

## Worked example: a split protomer, so ringfit can score the seam

Setup: SlyB is trimmed of its membrane belt, so each target protomer is two segments —
**19–59 (41 aa)** and **107–155 (49 aa)** — with a real 13.4 Å gap. Boltz was given
four target entities plus the binder, so `boltz_path` holds **five** chains: `A`+`C` =
protomer 1, `B`+`D` = protomer 2, `E` = binder (**check yours** — the predictor letters
them in the order `mkcomplex` built them). `ringfit` needs exactly two target chains, so
as it stands `bsa_t1`, `bsa_t2` and `bridge_ratio` cannot be computed at all.

```bash
# 1. merge each protomer into one chain, restoring the auth numbering,
#    and carry the binder through untouched as its own (degenerate) group.
sapia run     chainsel <run_dir> -t table1 -i boltz_path \
                  --merge-groups 'A+C:A,B+D:B,E:E' \
                  --renumber-from 'A:19-59,C:107-155,B:19-59,D:107-155' \
                  -l protomers
sapia collect chainsel <run_dir> -t table1 -l protomers

# expected on every row:
#   chainsel_protomers_n_chains = 3          (A, B, E — not 5)
#   chainsel_protomers_chains   = "A,B,E"
#   chainsel_protomers_n_res    = 90 + 90 + <binder length>

# 2. now ringfit has one chain per protomer, and the numbering the reference uses
sapia run     ringfit <run_dir> -t table1 \
                  -i chainsel_protomers_path \
                  --ref-structure inputs/slyb_assembly.cif \
                  --target-chains A,B --ref-target-chains A,B \
                  --resnum-match resnum
sapia collect ringfit <run_dir> -t table1
```

Why each piece:

- **`--merge-groups`, not `--rename-to`.** `--rename-to A,A,B,B` is a hard error by
  design — an accidental collapse must never be silent. The group spelling is the
  explicit request.
- **The binder needs its own group.** `--merge-groups` replaces `--chains` entirely, so
  any chain you do not name is **dropped**. `E:E` keeps the binder (ringfit's
  `--binder-chains auto` then picks it up as the non-target chain).
- **`--renumber-from`, not `--renumber`.** Boltz numbers each chain 1..N, so both
  segments start at 1 and the default (preserve) would collide — the run would come back
  `error: merged chain A ... has duplicate residue numbers`. `--renumber` would fix the
  collision but give 1..90, which only `ringfit --resnum-match ordinal` could use, and
  ordinal is wrong here: the reference chain also carries residues 60–106, so the i-th
  residues do not correspond. `--renumber-from` restores 19–59 / 107–155, and
  `--resnum-match resnum` pairs equal numbers and tolerates the gap.
- **The `-<last>` is a checksum.** `C:107-155` asserts 49 residues; if Boltz returned a
  different count that design errors instead of shifting the numbering by a few.
- **Order is the group's order.** `A+C:A` writes 19–59 then 107–155. Writing `C+A:A`
  would write 107–155 first — legal, but then the chain's numbers run backwards at the
  seam. Write the segments in sequence order.

Sanity-check one design before scaling: `chainsel_protomers_n_chains == 3`,
`chainsel_protomers_n_res == 180 + binder`, then after ringfit read
`ringfit_resnum_offset` (should be **0** — the merged chain now speaks the reference's
numbering) and `ringfit_seq_match_frac` (should be ~1.0). A non-zero offset or a poor
match means the merge or the ranges were wrong, and `bridge_ratio` would be
wrong-but-plausible.

> **gemmi is pinned to 0.7.5** in `modal_image.py` (Modal). Split mode's reader handling
> and every measurement in this skill are 0.7.5's, and the golden byte-identical test
> asserts the version. On a gemmi upgrade, re-check the reader behaviour (the canary
> test) and regenerate the golden md5s **deliberately, from the pre-split baseline
> worker** (`git show 3753f78:tools/chainsel/...`) run under the new gemmi — never from
> the current code. On vib the workspace venv gets whatever gemmi prosapia's
> `>=0.7.5,<0.8` resolves to; it is not pinned there.
>
> **Not yet run.** The tool was authored and its worker exercised locally on synthetic
> structures; it has never been submitted to Modal. Treat the first run as a shakedown —
> two designs, then read `chainsel_status` and `chainsel_n_res` before scaling.
>
> **`--split-groups` (added 2026-10-07) has not run on Modal, on vib or on a real
> RFdiffusion file.** It is tested only on synthetic structures written to mimic the
> described layout (`tests/test_chainsel.py`). Its first real run is a verification:
> read `chain_n_res` and `chain_mean_b` on every row and open one split file.
