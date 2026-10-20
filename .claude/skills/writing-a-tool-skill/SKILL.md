---
name: writing-a-tool-skill
description: >-
  The house shape for a prosapia tool's SKILL.md: the frontmatter contract (an unquoted colon silently strips the description), the required sections in order, and the pre-publish self-check. Load when writing or upgrading a tool skill, or when a skill's documented columns and its collector have drifted apart. Pairs with `authoring-a-tool`.
---

# Writing a tool skill

No directory under `tools/` carries a `.md`: the skill at `.claude/skills/<tool>/SKILL.md` is the tool's documentation, so a stale skill is an undocumented tool and a missing section is a question the next agent answers by guessing. `authoring-a-tool` covers the code; this covers the file that makes the code usable.

**A tool skill holds the facts of the tool**: its premise, how to invoke it, its columns with meaning and sign, its traps with their signatures, and the invariants a worker can check. It holds no campaign workflow (which designs, which thresholds, which order: that is the campaign log and the thinker), and it is edited only together with the verified code: never document a column, flag or behaviour that the committed source does not have.

## The frontmatter contract

```yaml
---
name: <tool>
description: >-
  <the question it answers>. Load before composing a <tool> run, and before reading <the column most likely to be misread>: <the one-clause trap>. Covers <the two or three things an agent cannot guess>.
---
```

- **Always the `>-` block scalar.** A plain value breaks at the first `: `; when the YAML fails to parse the skill still loads with no fields set and the description falls back to the first markdown line, hiding every trigger word. `>-` needs no quoting or escaping.
- **Say when to load it, and name the trap.** "How to run X" does not warn an agent that it is about to misread a column; the description is the only text read before the decision to load.
- **Keep it tight.** `description` plus `when_to_use` is truncated at 1,536 characters in the listing, and the whole listing has a character budget beyond which Claude Code drops the descriptions of the least-invoked skills, so a long one can cost another skill its description. Detail belongs in the body.

## The sections, in this order

1. **Identity, one line.** What it measures or builds, its `action` (`create` mints a child table, `update` annotates in place), whether `-t` is required, which executors it runs on.
2. **Premise, and when its numbers are meaningless.** The conditions the tool assumes (chains, what must be present, which frame), **what it does when they fail** (most tools here return a number anyway: say so, or say that it errors), and which columns are honest under which flags (a table of regimes naming the honest and the lying column in each). A premise only in the module docstring is read by nobody.
3. **Trust columns, read before the science.** Anything that aligns, maps, renumbers or matches needs trust columns (identity, numbering offset, alignment RMSD, count agreement); say to read them first and give the expected value ("about 1.0 against its own backbone"), not just the name. Include the generator conventions that make them necessary (rfd3 renumbers every chain from 1).
4. **The question and the threshold that belongs to it.** If one tool serves two gates under different flags, separate them and attach each threshold to its own gate; a bar correct for one mode and inverted for another is carried across. Say which mode a conventional bar is conventional for.
5. **Invocation and flags.** A verified command line, then every flag with its **default**; bold a default that is wrong for the common case. Include resources and whether the builder forces `gpus_per_task = 0`, and label every executor-specific statement ("Modal only").
6. **Traps, each with its signature.** What you will see, not a worry: "`Submitting N designs` comes back far under the row count because `--col-b` resolved to nothing on those rows" is a diagnosis; "rows can be dropped" is not.
7. **Columns collected**, as a table of name, meaning, **sign and direction** (which way is better, which sign is wanted, units), and how to read it. The third and fourth columns are what make the table worth having. Mark trust columns, say which is `_status` and that `OK` is the only proof of success, and note the on-disk layout if a caller will open raw output. **The list must match the collector name for name**, and the sign and direction must be verified against the source: a documented column the code does not write sends every caller to a filter that selects nothing, and a recalled sign is a silent error.
8. **A filter expression** a caller would write (`def apply_filter(df): return df[(df["<leaf>_status"] == "OK") & (df["<leaf>_<metric>"] < 2.0)]`). If the column set cannot support one, it has not been thought through. If a threshold needs a campaign property (a fixed target length), say so and point to the campaign log.
9. **Invariants a worker can check.** The checks to run on a first real leaf under a verification label (a hand reproduction of one number, an identity between two columns, a count against raw files, a break-it input that must error), marking which are real and which are tautologies (see `verifying-a-tool`).
10. **What it is blind to.** What a caller will wrongly believe the tool ruled out (geometry is blind to chemistry, fold scores to pose, confidence to whether it is about the right thing), naming the tool that answers each gap. Describe each column's blindness here, not only the tool's.

## Self-check before publishing

1. The frontmatter parses: `python -c "import yaml,sys; yaml.safe_load(open(sys.argv[1]).read().split('---')[1])" SKILL.md`.
2. The description says when to load and names the trap.
3. The premise says what happens when it is violated.
4. Trust columns are named, ordered first, with expected values.
5. Every flag has its default; every wrong-for-the-common-case default is bolded.
6. Every column has a how-to-read entry and a sign; the table matches `collect_<tool>.py`, which you read rather than recalled.
7. The filter expression would run against those names and dtypes.
8. Every executor-specific statement is labelled.
9. The blind-to section names the tool that answers each gap.
10. No campaign workflow, run names or thresholds-as-decisions have crept in.

## When updating rather than writing

- **Say which columns change meaning.** An edit to a default, cutoff, sign convention or unit makes rows collected before it non-comparable; record that in the skill as a note on the column and in the campaign log.
- **When a tool is removed, hunt the references.** Removing `tools/<name>/` and `.claude/skills/<name>/` is two of three: a tool named inside another skill's prose survives both, so `grep -rn '<name>' .claude/skills/` before calling the removal done.
- **Edit the skill with the committed code, not before it.** A skill that documents an edit living only in a working tree misleads every session that does not have that tree; if a review copy of `.claude/` exists, write both.
