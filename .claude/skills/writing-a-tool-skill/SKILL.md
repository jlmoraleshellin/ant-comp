---
name: writing-a-tool-skill
description: >-
  The house shape for a prosapia tool's SKILL.md: the frontmatter contract (an unquoted colon silently strips the description), the required sections in order, and the pre-publish self-check. Load when writing or upgrading a tool skill, or when a skill's documented columns and its collector have drifted apart. Pairs with `authoring-a-tool`.
---

# Writing a tool skill

**None of the directories under `tools/` contains a `.md` file.** The skill at
`.claude/skills/<tool>/SKILL.md` *is* the tool's documentation. A stale skill is an undocumented
tool, and a missing section is a question the next agent answers by guessing.

`authoring-a-tool` covers the code. This covers the file that makes the code usable.

## The frontmatter contract

```yaml
---
name: <tool>
description: >-
  <what it measures>. Load before composing a <tool> run, and before reading <the column most
  likely to be misread>: <the one-clause trap>. Covers <the two or three things an agent cannot
  guess>.
---
```

**Use the `>-` block scalar, always.** A plain unquoted value breaks the moment the text contains
`: ` — and when the YAML fails to parse, **Claude Code still loads the skill but with no fields set**,
so the description falls back to the first non-empty markdown line. Measured in this repo: `boltz`
advertised itself to the model as the seven characters `# boltz`, and `rpxdock` as its `#` heading —
in both cases every trigger keyword the author wrote was invisible in the skill listing, while the
file itself looked fine to a human reading it. `>-` removes the whole class: no quoting, no escaping,
colons safe.

Three more rules for the description:

- **Say when to load it, not just what the tool is.** "How to run X" does not tell an agent that it
  is about to misread a column. Put the key use case first.
- **Keep it tight.** `description` plus `when_to_use` is truncated at **1,536 characters** in the
  listing, and the whole listing has a character budget (1% of the context window by default) beyond
  which Claude Code **drops descriptions from the skills you invoke least**. A long description does
  not just cost context; past the budget it can cost another skill its description entirely. Detail
  belongs in the body, which loads only when the skill is invoked.
- **Name the trap in the description.** It is the one sentence that is read before the decision to
  load, so it is the only place a warning reaches an agent that was not going to read the file.

## The sections, in this order

The order is the order an agent needs them in, which is not the order they are easiest to write in.

### 1. One line of identity

What it compares/measures/builds, its **`action`** (`create` mints a child table, `update` annotates
in place), and whether `-t` is required.

### 2. Premise, and when its numbers are meaningless

**The section most often missing and the one that prevents the expensive failures.** State:

- the conditions the tool assumes — how many chains, what must be present, what frame the input is in;
- **what it does when those conditions do not hold.** Almost every tool here returns a number
  anyway. Say so explicitly: *"it always returns a value; there is no input for which it reports that
  the comparison is meaningless."*
- which of its own columns are honest **under which flags**. If the answer depends on a flag's value,
  give it as a table of regimes with the honest column and the lying column named in each.

A tool whose premise is only in the module docstring is a tool whose premise nobody reads.

### 3. Trust columns — read these before the science

If the tool aligns, maps, renumbers or matches anything, it must have trust columns (sequence
identity, residue-numbering offset, alignment RMSD, count agreement), and this section must tell the
reader to **read them first**. Give the expected value, not just the name: "~1.0 when comparing a
design to its own backbone" is actionable; "sequence identity" is not.

Include the generator conventions that make trust columns necessary — e.g. **rfd3 renumbers every
output chain from 1**, so a target numbered 18–155 comes back as 1–138.

### 4. The question(s) it answers, and the threshold that belongs to each

If one tool serves two gates under different flags, **separate them and attach each threshold to its
own gate.** A threshold that is correct for one mode and inverted for another is worse than no
threshold, because it will be carried across. Where a bar is conventional, say which mode it is
conventional *for*.

### 5. Invocation and flags

A verified command line, then every flag with its **default** — defaults are what a caller gets when
they do not think about it, so a default that is wrong for the common case is itself a trap and
should be bolded. Include resources and whether the manifest builder forces `gpus_per_task = 0`.

Mark anything **executor-specific**. A statement about Volumes, images, `modal-shell`, `--gpu-type`,
image build times or `.exit` files is Modal-only, and a worker on a cluster will read this file. Label
it rather than deleting it: *"(Modal only)"*.

### 6. Gotchas, each with its signature

Not a list of warnings — a list of **what you will see** when each one happens. "Rows can be dropped
silently" is a worry; "`Submitting N designs` comes back well under the table's row count because
`--col-b` resolved to nothing on those rows and all their ancestors" is a diagnosis.

### 7. What it collects

A table: column, meaning, **how to read it**. The third column is the one that makes the section
worth having. Mark trust columns as such. Say which column is `_status` and that `OK` is the only
proof of success. Note on-disk layout if a caller will ever open the raw output.

**This list must match the collector's actual output, name for name.** A documented column the code
does not write sends every caller to a filter that selects nothing.

### 8. At least one filter expression

```python
def apply_filter(df):
    return df[(df["<leaf>_status"] == "OK") & (df["<leaf>_<metric>"] < 2.0)]
```

A tool exists so its numbers can gate a decision. If the skill cannot show the `-f` module a caller
would write, the column set has not been thought through — this is the same question
`authoring-a-tool` requires answered before the tool is built, and the answer belongs here where the
caller will look for it.

If a threshold needs a value that is a property of the campaign rather than of the table (a fixed
target length, a reference count), say so and tell the reader to record the number in the campaign
log.

### 9. What it is blind to

What a caller will wrongly believe this tool has ruled out. Geometry tools are blind to chemistry;
fold scores are blind to pose; confidence is blind to whether the thing it is confident about is the
thing you cared about. Name the tool that *does* answer each one, so the section routes rather than
just warning.

## Self-check before you publish

Run through this; it is the same list `tool-reviewer` applies to the pair.

1. Does the frontmatter parse? `python -c "import yaml,sys; yaml.safe_load(open(sys.argv[1]).read().split('---')[1])" SKILL.md`
2. Does the description say **when** to load, and name the trap?
3. Does the premise section say what happens when the premise is violated?
4. Are the trust columns named, ordered first, and given expected values?
5. Is every flag's default present, and every wrong-for-the-common-case default bolded?
6. Does the column table have a "how to read it" entry for every row?
7. **Does that table match the collector?** Read `collect_<tool>.py`, not your memory of it.
8. Is there a filter expression that would actually run against those names and dtypes?
9. Is every Modal-only statement labelled?
10. Does the blind-to section name the tool that answers each gap?

## When you are updating rather than writing

- **Say which columns change meaning.** An edit to a default, a cutoff, a sign convention or a unit
  makes every row collected under the old behaviour non-comparable with the new ones. That belongs in
  the campaign log as a property of the column, and in the skill as a dated note.
- **If a tool was removed, hunt the references.** Removing `tools/<name>/` and
  `.claude/skills/<name>/` is two of three: a tool named inside *another*
  skill's prose survives both. Measured: `framefit` was deleted in full and remained the
  recommended tool in `binder-campaign` for the pose gate. `grep -rn '<name>' .claude/skills/` before
  calling a removal done.
- **Keep the mirror in step.** If the repo carries a review copy of `.claude/` (`agent-context/`),
  write both. A skill that exists in only one of the two will be edited in the wrong copy later.
