---
name: verifying-a-tool
description: >-
  The checklist a new or edited prosapia tool must pass before a campaign uses its columns: a code review that runs nothing, a real-data run under its own label, an independent recompute of the headline quantity, a deliberately bad input that must error, and a regression check against the pre-run table. Load when commissioning, editing or verifying a tool. A passing test suite is not verification.
---

# Verifying a tool

A tool's columns become gates, and a gate decides which designs get GPU hours, so a tool is done only when something other than what built it has checked it. The author wrote the tests from the same understanding as the code: they catch arithmetic, not a wrong premise, a wrong chain convention in the real outputs, or a column that is precisely computed and answers a different question. "The tests pass" is not a verification claim; say what the tests prove and what stays unknown until real data.

## The checklist, in cost order

1. **Review the code (`tool-reviewer`) before anything runs.** It reads `spec.py`, the premise docstring, builder, worker, collector, task script and skill, and hunts the defects real data reveals expensively: a plausible number on an inferred residue or chain mapping, a silent default where an error belongs, a trust metric emitted only on the success path, `0` where `NA` belongs, a name that claims more than the code measures. It returns findings, the invariant to use in check 2, and the assumptions that need real data; work from that list.
2. **Run it on real data under its own label.** A handful of designs from the real run_dir, never synthetic input and never the full set, with `-l verification` so columns and output dir stay out of the campaign's leaves. Require back: the collected columns as written; an invariant against a number the tool did not compute (a parent-table length, a residue identity at a known position, a count of raw result files, `n_res == target_len + binder_len`); the trust columns; the `<leaf>_status` counts, including any `missing` that should have been submitted. **Read the trust columns before the science columns**: a beautiful metric with `seq_match_frac` 0.6 is a metric about the wrong residues.
3. **Recompute the headline quantity independently.** A separate script from the raw inputs (structure files, raw tool outputs), not the tool's code path, compared row by row. Agreement of the tool with itself proves nothing.
4. **Break it on purpose.** A deliberately bad input must give an error status and `NA`, never a plausible number; check the table, not the exit code. A tool returning `0.0` for an undefined measurement passes a `< threshold` filter and looks like the best design. Pick the failure the premise excludes:

| premise | the bad input |
| --- | --- |
| two adjacent protomers from an assembly | a single-chain structure |
| a predicted complex with a matching reference | a reference with a different chain count |
| a binder plus a target | a monomer with the target chain absent |
| a structure with a ligand belt | the same structure with ligands stripped |
| matching per-chain atom counts | a truncated chain |

For a tool that places or moves one structure relative to another, also re-run it with an input re-expressed in an arbitrary frame (non-trivial rotation, translation of 100+ A, applied to one file and then to the whole complex) and require every frame-independent column to be identical; a column that changes with where the file sits in space measures in the wrong frame. Routing noise from a grid is a tolerance stated as a number; a different answer is not.
5. **Regression.** Back up the table before the run; afterwards every pre-existing column must be identical to the backup (an additive change touches nothing else). An edit that alters a default, a cutoff, a sign convention or units makes old and new rows non-comparable: say which columns change meaning, and give the edited tool a new `-l` label rather than mixing two behaviours under one leaf.
6. **Use invariants that are not tautologies.** A quantity defined as a sum cannot verify itself, and a partition whose parts and total come from the same loop passes by construction. Say in the report which checks are real (independent path, different inputs) and which are true by construction.
7. **Check direction and sign on a known example.** Choose a case whose answer is known in advance and confirm each column moves the way its documented sign says; a sign error passes every other check.
8. **Compare inputs of different provenance only like-for-like.** Predicted versus relaxed coordinates, or two recipes, give different numbers for the same design; split the table by provenance before ranking or correlating.

## Re-review a fix, and re-review it on the same function

**A fix is new code and inherits none of the review that preceded it.** Measured on `hbcheck`, 2026-10-17: pass 1 found a collector whose early returns carried no `data`, so a failed rerun kept the previous run's numbers beside its failure status. The fix for it swept columns by leaf-name prefix — and because a labelled leaf is `<tool>_<label>`, the unlabelled leaf's prefix is a **strict prefix of every labelled leaf's columns**. The first production collect would have nulled the `verification` leaf's science *and its status*, destroying the evidence the tool was trusted on, on the very workflow this page prescribes. Pass 3, scoped to that one function, found a third defect beside it: a sanitiser that rewrote `H|H` into the valid-but-unwritten `HH` instead of refusing the cell.

Three lessons, in order of how much they cost:

- **Column ownership is the collector's hardest problem.** Anything that decides which columns to blank, null, carry or overwrite is where a plausible, test-invisible defect lives. Derive the set of names a leaf legitimately owns; never accept "everything after my prefix".
- **A fix in that code gets its own pass, narrowly scoped to the function.** A general re-review spreads attention; naming the one function concentrates it. If a third defect appears in the same place, that is the signal the approach is reaching past what it can do safely — simplify rather than patch again.
- **Prefer refusing to repairing.** When a cell is unusable, error. A sanitiser that strips characters until the value parses hands back a number nobody wrote, which is the same failure class as a silent default.

None of the three was reachable by a test suite: each needed a reader asking what the code owns and what it would do on the *next* command.

## What "verified" means

The code was read by something that did not write it; the columns were produced on real designs from the real run_dir; an independent recompute and an independent invariant held; the trust columns said the mapping was right; a violation of the premise produced an error; and nothing that existed changed. It does not mean the tool is right for the next target: a tool verified on a trimmed oligomer is unverified on a single-chain target, and reusing it there is a scope decision.

## Recording it

Verification is a campaign step: it goes to the `tracker` with its run_dir, its `verification` label, the invariant and whether it held. If it fails, report the raw evidence and hand it back to whoever built it; do not patch from the verification seat, because the fix needs the agent that restates the premise, and the checklist again.

## After a verified edit

Tell the user to commit it the moment it is verified (or copy it somewhere safe), say in the hand-back whether it is committed, and edit the tool's skill only together with the commit: a working-tree edit that a skill already documents is a trap the next session cannot see.
