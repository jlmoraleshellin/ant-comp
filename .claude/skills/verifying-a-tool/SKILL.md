---
name: verifying-a-tool
description: >-
  The three checks a new or edited prosapia tool must pass before a campaign uses its columns: a code review that runs nothing, a real-data run under a verification label, and a deliberately bad input that must error. Load when commissioning, editing or verifying a tool. A passing test suite is not verification.
---

# Verifying a tool

A tool's columns become gates, and a gate decides which designs get GPU hours. So a tool is not "done" when it runs — it is done when something **other than the thing that built it** has checked it.

## Why the test suite does not count

The agent that wrote the tool wrote the tests, from the same understanding of the problem. Tests catch arithmetic and they document intent, both worth having. They cannot catch a misunderstanding of the **premise**, a wrong chain convention in the real outputs, or a column that is precisely computed and answers a different question.

*The precedent in this workspace:* a tool shipped with **78 of 78 checks passing** — a binder translated by a known vector, a rotated complex agreeing to 0.0006, an independent quaternion Kabsch written in plain numpy agreeing to 1e-3 — while it was still unknown whether the real Boltz and rfd3 outputs carry the chain labels and per-chain CA counts it assumed. The arithmetic was right. Whether it was arithmetic about the right atoms was untested.

So: **"the tests pass" is not a verification claim.** Say what the tests prove and what remains unknown until real data.

## The three checks, in cost order

### 1. Review the code — free, finds the most

`tool-reviewer` reads `spec.py`, the premise docstring, the manifest builder, the worker, the collector, the task script and the skill. It runs nothing, costs no GPU, and it is looking for exactly the class of defect that real data will otherwise reveal expensively: a plausible number on an inferred residue or chain mapping, a silent default where an error belongs, a trust metric that is only emitted on the success path, `0` written where `NA` belongs, a column whose name claims more than the code measures.

It returns findings by severity plus **the invariant it recommends** for check 2 and **the assumptions that cannot be checked without real data**. Do step 2 against that list rather than inventing your own.

### 2. Run it on real data, under a `verification` label

Not on synthetic input, and not on the full set. A handful of designs from the actual run_dir, with `-l verification` so the columns and the output directory are quarantined from the campaign's own leaves.

Order from the worker, and require all four back:

1. **the collected columns** for those designs, as actually written to the table;
2. **an invariant checked against a number the tool did not compute** — a length from a parent table, a residue identity at a known position, a count of raw per-design result files, an arithmetic identity like `n_res == target_len + binder_len`;
3. **the trust columns**, confirming the mapping was right: sequence-identity fraction at or near 1.0 where it should be, a numbering offset of 0 where it should be, an alignment RMSD in the range the premise implies;
4. **the `<leaf>_status` counts**, plus whether any row came back `missing` that should have been submitted.

**Read the trust columns before the science columns.** A beautiful metric with `seq_match_frac` of 0.6 is a metric about the wrong residues, and the order you read them in is what stops you believing it.

### 3. Break it on purpose

**A deliberately bad input must produce an error status, not a plausible number.** This is the check that distinguishes a tool from a function, and it is the one most often skipped because the tool "obviously works".

Pick a failure the premise is supposed to exclude, and give it exactly that:

| premise | the bad input |
| --- | --- |
| two adjacent protomers from an assembly | a single-chain structure |
| a predicted complex with a matching reference | a reference with a different chain count |
| a binder plus a target | a monomer with the target chain absent |
| a structure with a ligand belt | the same structure with ligands stripped |
| matching per-chain atom counts | a truncated chain |

Then check the table, not the exit code: the row must carry an **error status and `NA`**, not `0`, not an empty string, and not a number. A tool that returns `0.0` for an undefined measurement will pass a `< threshold` filter and look like the best design in the batch.

## What "verified" means, and what it does not

Verified means: the code was read by something that did not write it, the columns were produced on real designs from this run_dir, an independent invariant held, the trust columns said the mapping was right, and a violation of the premise produced an error rather than a value.

It does **not** mean the tool is right for the next target. The premise is still a premise: a tool verified on a trimmed oligomer is unverified on a single-chain target, and reusing it there is a scope decision, not a verified one.

## Recording it

The verification is a campaign step like any other: it goes to the `tracker` with its run_dir, its `verification` label, the invariant and whether it held. A tool used in a campaign whose ledger has no verification step is a tool whose columns cannot be defended later.

If verification **fails**, report the raw evidence and hand it back to whoever built it. Do not patch the tool from the verification seat — the agent that fixes it should be the one that re-states the premise, and the fix then needs the three checks again.

## Editing an existing tool

The same three checks, with one addition: **say which existing columns change meaning.** If an edit alters a default, a cutoff, a sign convention or a column's units, every row already collected under the old behaviour is now non-comparable with the new ones — and that fact belongs in the campaign log as a property of the column, not as a changelog entry nobody reads. Give the edited tool a new `-l` label rather than silently mixing two behaviours under one leaf.
