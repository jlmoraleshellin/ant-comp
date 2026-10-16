"""Submit tasks running AtomiUM on designs, spawning a new child table.

AtomiUM is a ProteinMPNN-like sequence designer (a PyG re-implementation with
noise-conditioned weights), and it ships ProteinMPNN's own helper scripts under
``proteinMPNN_helper_scripts/`` (on Modal; off Modal the activation script supplies
the helper directory, see ``atomium.sh``). Its jsonl interface is the same one
(``--jsonl_path``, ``--chain_id_jsonl``, ``--fixed_positions_jsonl``,
``--tied_positions_jsonl``, ``--bias_AA_jsonl``), so this builder is deliberately
the proteinmpnn builder with the same two mini-languages and the same batching:
designs whose params are identical -- the signature ``(chains, fixed_positions,
tie_mode, tied_positions)`` -- are grouped into a single batched ``atomium.py``
call (model loaded once), and groups are bin-packed onto tasks up to
``--designs-per-task`` designs each.

The pipeline (``parse_multiple_chains`` -> optional helper steps -> ``atomium.py``)
runs inline in ``atomium.sh``, once per group.

Chain mini-language (``--chains-to-design``): the chains you design, in ProteinMPNN
order. ``:`` is an inclusive letter range, ``,`` separates:

    --chains-to-design A:C,E     # -> "A B C E"
    --chains-to-design A,C        # -> "A C" (non-contiguous)
    --chains-to-design ''         # design all chains (default)

Position mini-language (``--fixed-positions`` / ``--tied-positions``): ``/`` breaks
chains (its groups map one-to-one onto ``--chains-to-design``, in order), ``,``
separates fragments within a chain, ``start:end`` expands inclusively, a single
position passes through; table column expressions live in ``{...}`` islands (resolved
up the lineage with ``+ - * //`` arithmetic), everything else is a literal integer:

    --fixed-positions 9:23/10,11,18:20,22    # -> "9 10 ... 23, 10 11 18 19 20 22"
    --tied-positions  1:8/1:8                # -> "1 2 3 4 5 6 7 8, 1 2 3 4 5 6 7 8"

Positions are 1-indexed within each parsed chain (every chain is renumbered to
1..L), NOT original PDB numbering.

Symmetry (``--symmetry auto|N``): ties all chains via ``make_tied_positions_dict
--homooligomer 1``, and makes the position mini-language describe a single
asymmetric unit (a lone group is broadcast across the designed chains). ``auto``
takes the order from the input's polymer chain count. Mutually exclusive with
``--tied-positions``.

What differs from the proteinmpnn tool, all verified against the checkout:

  * ``--model-noise`` picks the weights (``model_noised_=_<noise>.pt``) and replaces
    ProteinMPNN's model-name/backbone-noise knobs. AtomiUM has no ``--ca_only``,
    ``--use_soluble_model``, ``--save_score`` or ``--save_probs``.
  * This tool exposes no ``--batch-size``. CORRECTED 2026-10-14: the commit this
    note used to cite for ``atomium.py``'s own ``--batch_size`` being LIVE
    (dee5a4d67ae02293cee59f86c86085c8fa11879d) is an **atomiumph** commit, not this
    tool's pinned checkout (``ATOMIUM_COMMIT = def1d6a...`` in ``modal_image.py``).
    Its status at the pinned commit is therefore UNVERIFIED -- it was DEAD at
    292e29c and LIVE somewhere later. Either way it affects CHUNKING ONLY: the
    sampling count is ``num_seq_per_target * len(temperatures)`` regardless, and
    this tool exposes no flag for it.
  * ``--sampling-temp`` may carry SEVERAL temperatures, and each one produces a full
    ``--num-seq-per-target`` set (so N temps => N x num_seq sequences per design).
    It therefore travels as its own manifest field and is quoted in the task script;
    folding it into the unquoted argv tail would split "0.1 0.2" into two arguments.
  * ``--seed 0`` means "pick a random seed" (AtomiUM treats 0 as falsy).

WHAT THE FASTA CONTAINS (measured, 2026-10-14 -- do not re-derive from memory)
--------------------------------------------------------------------------------
A 2-design probe on real output (run_dir ``20261009_130629_5M2J_vhh_redmotif_lowtemp_scaleup``,
child table ``table1_atomium_probe``, ``--chains-to-design D`` with 17 pinned
positions) settled the question this tool's skill used to get wrong:

  * the generated sequence is the FULL designed chain, pinned residues INCLUDED
    (``scaleup0_diff_1``: chain-D length 144, generated length 144; the 17 pinned
    characters read ``TSNYWYENTNLIKRSPF`` in both the native chain and every
    generated sequence);
  * ``_init_sequence`` is byte-identical to the native designed chain, NOT to the
    native minus the pinned positions.

So nothing is reinserted anywhere: ``atomium_sequence`` is the binder. The
collector's job is to PROVE that per row rather than assume it -- see
``collect_atomium.py``'s trust columns.

Like ProteinMPNN's, the FASTA covers the DESIGNED chains only. **Above ONE
designed chain this is assumption, not measurement**: the probe had a single
designed chain, so the ``/`` join and the order of the chains within it are read
from ProteinMPNN's printing convention and have never been checked against a real
multi-chain AtomiUM run. The collector is written to survive being wrong about the
order (its length and echo checks compare multisets) and reports
``init_order_matches`` so the first real multi-chain run settles it from the table.

Usage:
    sapia run atomium outputs/<run> --table table0 --table-label atomium \\
        --model-noise n05 --num-seq-per-target 8

    # scan two noise levels as two runs into two child tables
    sapia run atomium outputs/<run> -t table0 --table-label n03 --model-noise n03
    sapia run atomium outputs/<run> -t table0 --table-label n05 --model-noise n05
"""

import hashlib
import json
import shutil
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
from typing import cast

import pandas as pd

from prosapia.core import (
    CommonArgs,
    LookupFn,
    ManifestCtx,
)
from prosapia.core.executors import volume_path
from prosapia.utils import (
    ensure_pdb,
    expand_chain_spec,
    parse_positions,
    polymer_chain_lengths,
    polymer_chain_names,
    resolve_template,
)

# The weight files shipped in the checkout: model_weights/model_noised_=_<noise>.pt.
# Guarded up front because a typo would otherwise burn a GPU container per design
# before torch.load raises.
MODEL_NOISES = ("n00", "n01", "n02", "n03", "n04", "n05", "n05_v2", "n06", "n07")


class AtomiumArgs(CommonArgs):
    chains_to_design: str
    chains_to_design_column: str
    fixed_positions: str
    fixed_positions_column: str
    tied_positions: str
    symmetry: str | None
    bias_aa: str
    set: list[str]
    model_noise: str
    num_seq_per_target: int
    sampling_temp: str
    seed: int
    designs_per_task: int


def _parse_chains(spec: str) -> str:
    """Expand the chain mini-language into a space-separated ``--chain_list`` string."""
    return " ".join(expand_chain_spec(spec))


def _parse_positions(
    spec: str,
    lookup: LookupFn,
    name: str,
    *,
    broadcast_to: int | None = None,
    chain_lengths: list[int] | None = None,
    flag: str = "--fixed-positions",
) -> str:
    """Expand the position mini-language into a ``--position_list``.

    Per-chain groups joined by ``", "``, positions within a group by spaces. Order is
    preserved and NOT de-duplicated (tied positions are index-parallel). Empty spec
    -> "". ``broadcast_to`` (set from ``--symmetry``) replicates a spec describing a
    SINGLE asymmetric unit across that many chains; any other group count passes
    through untouched.

    ``chain_lengths`` is the residue count of each DESIGNED chain, in the order the
    groups map onto them. Passing it refuses, at submit time, two specs that are
    otherwise silently wrong:

      * a position outside ``1..length``. ``parse_chain_positions`` bounds-checks
        ONLY when it is handed a length -- with ``length is None`` its own docstring
        leaves bounds "to the caller" (``utils/positions.py``), which is how an
        out-of-range position used to reach the task;
      * a group count that does not match the designed chains.
        ``make_fixed_positions_dict.py`` zips ``--chain_list`` with the ``,``-joined
        groups, so an EXTRA group is dropped without a word (``--chains-to-design D
        --fixed-positions '49/12'`` pins only 49) and a MISSING one dies mid-task
        with an IndexError.

    It also enables open-ended ranges (``10:``, ``:50``, ``:``), which raise without
    a length.
    """
    # Resolved here so the group count can be counted before parsing; parse_positions
    # resolves again, which is a no-op on the brace-free string this produces.
    resolved = resolve_template(spec, lookup, name).strip().strip("[]").strip()
    if not resolved:
        return ""
    lengths: list[int] | None = None
    if chain_lengths is not None:
        n_groups = len(resolved.split("/"))
        if broadcast_to and n_groups == 1:
            # One group describing a single asymmetric unit, replicated below across
            # `broadcast_to` chains: it has to be in bounds for EVERY one of them.
            lengths = [min(chain_lengths)]
        elif n_groups != len(chain_lengths):
            raise ValueError(
                f"design {name!r}: {flag} has {n_groups} '/'-separated chain "
                f"group(s) but {len(chain_lengths)} chain(s) are designed. Groups "
                f"map one-to-one onto the designed chains, in order; the helper "
                f"zips them, so an extra group is dropped SILENTLY (pinning fewer "
                f"positions than you wrote) and a missing one kills the task. "
                f"Spec: {resolved!r}"
            )
        else:
            lengths = list(chain_lengths)
    groups = parse_positions(resolved, lookup, name, lengths=lengths)
    if broadcast_to and len(groups) == 1:
        groups = groups * broadcast_to
    return ", ".join(" ".join(str(p) for p in group) for group in groups)


def _check_symmetry_spec(spec: str | None) -> None:
    """Validate ``--symmetry`` up front, before anything is submitted."""
    if spec is None or spec == "auto":
        return
    try:
        order = int(spec)
    except ValueError:
        raise ValueError(
            f"--symmetry: expected 'auto' or an integer chain count, got {spec!r} "
            f"(write the number of chains in the oligomer, e.g. --symmetry 12)"
        )
    if order < 2:
        raise ValueError(f"--symmetry: order must be at least 2, got {order}")


def _check_model_noise(noise: str) -> None:
    """Reject a noise level with no weight file, before any container starts."""
    if noise not in MODEL_NOISES:
        raise ValueError(
            f"--model-noise: {noise!r} has no weights in the AtomiUM checkout; "
            f"choose one of {', '.join(MODEL_NOISES)}"
        )


def _check_sampling_temp(spec: str) -> None:
    """Every ``--sampling-temp`` token must be a float; AtomiUM splits on whitespace."""
    tokens = spec.split()
    if not tokens:
        raise ValueError("--sampling-temp: at least one temperature is required")
    for token in tokens:
        try:
            float(token)
        except ValueError:
            raise ValueError(
                f"--sampling-temp: {token!r} is not a number (expected one or more "
                f"space-separated temperatures, e.g. '0.1' or '0.1 0.2')"
            )


def _cell_text(value: object) -> str:
    """One table cell as a stripped string; every missing form becomes ``""``.

    ``pd.isna`` first, as ``tools/proteinmpnn/run_proteinmpnn.py`` does: stringifying
    before testing turns ``pd.NA`` into the literal ``"<NA>"``, which then sails
    through as a spec and dies much later with a misdirecting message. The string
    forms are kept as a second net for a TSV that round-tripped a null as text.
    """
    if isinstance(value, str):
        text = value.strip()
    elif value is None or pd.isna(value):
        return ""
    else:
        text = str(value).strip()
    return "" if text.lower() in ("nan", "none", "<na>") else text


def _explicit_chains(spec: str, source: str, design_name: str = "") -> str:
    """Expand an EXPLICIT chain spec, refusing one that expands to NO chains.

    ``"[]"``, ``","`` and ``"[ ]"`` are all non-empty strings that
    ``expand_chain_spec`` turns into ``[]``. An empty ``--chain_list`` does not mean
    "no chains" anywhere downstream -- it means "design every chain": ``atomium.sh``
    writes no ``assigned_pdbs.jsonl`` at all, and still calls
    ``make_fixed_positions_dict --chain_list ""``, which pins NOTHING. The task exits
    0 and the row collects ``status OK`` with ``n_fixed == 0`` and ``fixed_held`` NA
    -- a motif-graft design whose motif was never held, indistinguishable in the
    table from a legitimate no-pin run. So an explicit source that expands to nothing
    is refused here.
    """
    chains = _parse_chains(spec)
    if not chains:
        who = f"design {design_name!r}: " if design_name else ""
        raise ValueError(
            f"{who}{source} is {spec!r}, which expands to NO chains. An empty chain "
            f"list means 'design EVERY chain' to the helper scripts and pins "
            f"nothing, at exit 0 with status OK -- so it is refused rather than "
            f"submitted. Write the chain letter(s), or omit the flag entirely if "
            f"designing every chain is really what you meant."
        )
    return chains


def _grp_digest(
    chains: str, fixed_pl: str, tie_mode: str, tied_pl: str, names: list[str]
) -> str:
    """The subgroup directory's name, derived from WHAT IS IN IT, not a counter.

    Unchanged since the staging fix, payload included: changing it would rename
    every grp_ dir already on disk and make this submit look like a cross-subgroup
    restage of designs that never moved.
    """
    payload = "\n".join([chains, fixed_pl, tie_mode, tied_pl, *sorted(names)])
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _text_digest(text: str) -> str:
    """A short, stable content digest of a file's own bytes (64 bits of sha256)."""
    return hashlib.sha256(text.encode()).hexdigest()[:16]


_CHAIN_LENGTH_CACHE: dict[str, dict[str, int]] = {}


def _chain_lengths_by_name(pdb_src: Path) -> dict[str, int]:
    """``{chain name: residue count}`` for the first model's polymer chains.

    Cached per path: the same structure is asked about twice (chain validation and
    position bounds), and a wide table would otherwise re-read every file.

    This is the SUBMITTER's view of the structure (gemmi, via prosapia's own
    helpers). The task's view is ``parse_multiple_chains.py``'s, which counts the
    residues carrying the backbone atoms it needs. The two agree on an ordinary
    design PDB; where they could not, this is a GUARD and not a guarantee -- a
    refusal here is real, but a spec that passes here can still be wrong at task
    time, which is what the collector's ``len_matches_parsed`` exists to catch.
    """
    key = str(pdb_src)
    if key not in _CHAIN_LENGTH_CACHE:
        _CHAIN_LENGTH_CACHE[key] = dict(
            zip(polymer_chain_names(pdb_src), polymer_chain_lengths(pdb_src))
        )
    return _CHAIN_LENGTH_CACHE[key]


def _validate_chains(chains: str, pdb_src: Path, design_name: str, source: str) -> None:
    """Refuse a designed-chain letter this design's own structure does not have.

    Nothing downstream treats it as an error: ``assign_fixed_chains.py`` copies
    ``--chain_list`` verbatim, so a chain the file lacks designs NOTHING, every real
    chain is marked fixed, the task exits 0 and the row collects ``OK`` carrying the
    input sequence back. The letter is the author chain id as it appears in the PDB
    -- the same one ``parse_multiple_chains.py`` reads.
    """
    present = _chain_lengths_by_name(pdb_src)
    missing = [c for c in chains.split() if c not in present]
    if missing:
        raise ValueError(
            f"design {design_name!r}: {source} names chain(s) "
            f"{', '.join(missing)}, which are not polymer chains of {pdb_src} "
            f"(present: {', '.join(present) or '<none>'}). A chain letter the "
            f"structure does not have designs nothing and still collects OK, so "
            f"the submit is refused instead."
        )


def _column_specs(ready, column: str, flag: str) -> dict[str, str]:
    """Per-design spec strings read verbatim from ``column``, for every ready design.

    Validates the WHOLE batch in one pass and raises ONCE, naming every offending
    design, so the table is repaired in a single edit instead of one resubmit per
    bad row. A missing value is REFUSED, never defaulted: an empty spec would sail
    through as ``""``, produce no jsonl, and silently redesign the whole chain at
    exit 0.
    """
    if column not in ready.columns:
        available = ", ".join(sorted(str(c) for c in ready.columns)) or "<none>"
        raise ValueError(
            f"{flag} {column!r} is not a column of the source table. "
            f"Available columns: {available}"
        )
    specs: dict[str, str] = {}
    missing: list[str] = []
    for name in sorted(cast(str, n) for n in ready.index):
        text = _cell_text(ready.at[name, column])
        if not text:
            missing.append(name)
        else:
            specs[name] = text
    if missing:
        shown = ", ".join(missing[:10])
        more = f" (and {len(missing) - 10} more)" if len(missing) > 10 else ""
        raise ValueError(
            f"{flag} column {column!r} is missing or empty for "
            f"{len(missing)} of {len(ready.index)} ready design(s): {shown}{more}. "
            f"Refusing the submit rather than redesigning those chains in full. "
            f"Use -f/--filter to exclude them, e.g. "
            f"'df[df[\"{column}\"].notna() & (df[\"{column}\"] != \"\")]'."
        )
    return specs


def _resolve_restage(
    out_dir: Path, planned_dir_by_design: dict[str, Path], force: bool
) -> None:
    """Refuse (or, under ``--force``, clean up) a design staged elsewhere here.

    PORTED VERBATIM IN BEHAVIOUR from ``tools/proteinmpnn/run_proteinmpnn.py``
    (cross-tool imports are impossible: ``tool_registry._load`` roots each tool
    folder in its own synthetic package). Keep the two in step.

    With content-addressed names a design lands in a DIFFERENT directory than last
    time only if its signature or its subgroup membership changed since an earlier
    submit to this same table+leaf. The earlier directory keeps its outputs, and the
    collector indexes every ``*/seqs/*.fa`` under ``out_dir`` by FASTA stem -- so two
    directories holding one design means two FASTAs and whichever was read last wins.
    A silent wrong answer, so it is handled here, at submit time.

      * ``--force`` and every design in the stale directory is restaged by THIS
        submit -> its outputs are wholly superseded, so it is removed (loudly).
      * anything else -> raise, naming every design and directory. Deleting would
        destroy a result the caller did not ask to redo.
    """
    members_by_dir: dict[Path, set[str]] = {}
    dir_by_design: dict[str, Path] = {}
    for inputs_dir in sorted(out_dir.glob("grp_*/inputs")):
        grp_dir = inputs_dir.parent
        for link in sorted(inputs_dir.glob("*.pdb")):
            dir_by_design[link.stem] = grp_dir
            members_by_dir.setdefault(grp_dir, set()).add(link.stem)

    clashes = [
        (name, dir_by_design[name], planned_dir_by_design[name])
        for name in sorted(planned_dir_by_design)
        if name in dir_by_design
        and dir_by_design[name] != planned_dir_by_design[name]
    ]
    if not clashes:
        return

    restaged = set(planned_dir_by_design)
    removable = {
        old
        for _, old, _ in clashes
        if force and members_by_dir.get(old, set()) <= restaged
    }
    blocking = [(n, old, new) for n, old, new in clashes if old not in removable]

    if blocking:
        shown = "; ".join(
            f"{name}: already staged in {old.name}, this submit would use {new.name}"
            for name, old, new in blocking[:5]
        )
        more = f" (and {len(blocking) - 5} more)" if len(blocking) > 5 else ""
        hint = (
            "Re-run with --force to replace them"
            if not force
            else "Those dirs also hold designs this submit is not redoing, so they "
            "cannot be replaced wholesale"
        )
        raise ValueError(
            f"{len(blocking)} design(s) are already staged under this output dir "
            f"({out_dir}) in a DIFFERENT subgroup than this submit would use: "
            f"{shown}{more}. The same design was submitted to this table+leaf "
            f"before with a different chains/fixed-positions/tied signature or in "
            f"a different batch; keeping both would leave two FASTAs for one design "
            f"and the collector would pick one arbitrarily. {hint}, or use a fresh "
            f"-l/--dir-label for this variant, or delete the stale grp_ dir(s) "
            f"named above."
        )

    for old in sorted(removable):
        print(
            f"atomium: --force, removing superseded subgroup {old} "
            f"({len(members_by_dir[old])} design(s), all restaged by this submit)"
        )
        shutil.rmtree(old)


def _stage_inputs(inputs_dir: Path, members: list["Member"]) -> None:
    """Stage exactly ``members`` into ``inputs_dir`` -- no more, no less.

    PORTED from ``tools/proteinmpnn/run_proteinmpnn.py``; keep the two in step.
    ``_symlink_member`` (still the thing that creates each link) only ever APPENDS:
    it skips a link that exists and never removes one, so a changed target was
    silently ignored and a stray entry silently joined the batched call. Here the
    directory is made to MATCH the member list:

      * an entry that is not an expected member is removed (with content-addressed
        names this cannot arise from a legitimate submit; it guards a hand-edited or
        interrupted directory, and it prints when it fires);
      * a symlink pointing somewhere other than this member's staged PDB is
        repointed, so a re-run whose input path moved does not silently redesign the
        old structure.
    """
    inputs_dir.mkdir(parents=True, exist_ok=True)
    expected = {f"{name}.pdb": Path(pdb_src) for name, pdb_src in members}

    for entry in sorted(inputs_dir.iterdir()):
        if entry.name not in expected and (entry.is_symlink() or entry.is_file()):
            print(
                f"atomium: removing stale staged input {entry} "
                f"(not a member of this subgroup)"
            )
            entry.unlink()

    for name, pdb_src in members:
        link = inputs_dir / f"{name}.pdb"
        if link.is_symlink():
            if link.readlink() == Path(pdb_src):
                continue
            print(
                f"atomium: repointing staged input {link} "
                f"({link.readlink()} -> {pdb_src})"
            )
            link.unlink()
        elif link.exists():
            link.unlink()
        _symlink_member(Path(pdb_src), inputs_dir, name)


def _resolve_chains(args: AtomiumArgs, pdb_src: Path) -> str:
    """One design's designed-chain list, as a ``--chain_list`` string.

    ``--chains-to-design`` wins when given. Otherwise ``--symmetry`` derives the list
    from the structure itself -- the first ``order`` polymer chains, where ``auto``
    takes every one. With neither flag, "" leaves every chain designed.
    """
    if args.chains_to_design:
        return _explicit_chains(args.chains_to_design, "--chains-to-design")
    if args.symmetry is None:
        return ""
    names = polymer_chain_names(pdb_src)
    order = len(names) if args.symmetry == "auto" else int(args.symmetry)
    return " ".join(names[:order])


def _write_bias_jsonl(spec: str, out_dir: Path) -> str:
    """Write a global AA-bias jsonl from an ``AA:bias`` spec; return its path.

    ``make_bias_AA.py`` is purely generative (no structure input), so the dict is
    written directly and the AtomiUM env is not needed on the submit node.
    """
    spec = spec.strip()
    if not spec:
        return ""
    bias: dict[str, float] = {}
    for token in spec.split():
        if ":" not in token:
            raise ValueError(
                f"--bias-aa: malformed pair {token!r} (expected 'AA:bias', "
                f"e.g. 'D:1.39')"
            )
        aa, value = token.split(":", 1)
        try:
            bias[aa] = float(value)
        except ValueError:
            raise ValueError(f"--bias-aa: bias for {aa!r} is not a number: {value!r}")
    # Content-addressed, for the reason grp_ dirs and task files are: a fixed name
    # is REWRITTEN by the next submit to the same leaf, so a second arm with a
    # different --bias-aa silently replaces the file the first arm's rows were
    # designed with, and the audit (`cat` it) then reports the wrong bias for
    # already-collected rows. Its path travels in the task file, which is itself
    # content-addressed, so nothing has to agree on the name.
    payload = json.dumps(bias, sort_keys=True) + "\n"
    bias_path = out_dir / f"bias_AA_{_text_digest(payload)}.jsonl"
    bias_path.write_text(payload)
    return str(bias_path)


def _atomium_extra(args: AtomiumArgs, bias_path: str) -> str:
    """The run-wide atomium.py argv tail (same for every design).

    Whitespace-free by construction: the task script expands this unquoted so each
    token becomes its own argv entry. ``--sampling_temp`` is NOT here -- it may hold
    several space-separated temperatures and travels as its own quoted field.
    """
    parts = [
        "--model_noise",
        args.model_noise,
        "--num_seq_per_target",
        str(args.num_seq_per_target),
        "--seed",
        str(args.seed),
    ]
    if bias_path:
        parts += ["--bias_AA_jsonl", bias_path]
    parts.extend(args.set)
    return " ".join(parts)


# A design ready to run, and the (chains, fixed_positions, tie_mode, tied_positions)
# signature that decides which designs can share one batched atomium.py call.
Member = tuple[str, Path]  # (design_name, staged_pdb_source)
Signature = tuple[str, str, str, str]
Subgroup = tuple[Signature, list[Member]]


def _pack_subgroups(chunks: list[Subgroup], max_per_task: int) -> list[list[Subgroup]]:
    """First-fit-pack subgroup chunks (each <= max_per_task) into tasks.

    Each task holds at most ``max_per_task`` designs total, possibly spanning several
    signatures; every subgroup stays intact so batching is preserved.
    """
    tasks: list[list[Subgroup]] = []
    remaining: list[int] = []
    for chunk in chunks:
        size = len(chunk[1])
        for t, cap in enumerate(remaining):
            if size <= cap:
                tasks[t].append(chunk)
                remaining[t] -= size
                break
        else:
            tasks.append([chunk])
            remaining.append(max_per_task - size)
    return tasks


def _symlink_member(pdb_src: Path, inputs_dir: Path, name: str) -> None:
    """Symlink a staged design into ``inputs_dir`` as ``<name>.pdb`` (no copy).

    The stem is what collect keys lineage off: AtomiUM names each output FASTA after
    the parsed input, so ``<name>.pdb`` -> ``seqs/<name>.fa`` -> parent row ``name``.
    """
    link = inputs_dir / f"{name}.pdb"
    if not link.exists() and not link.is_symlink():
        link.symlink_to(pdb_src)


def build_atomium_manifest(ctx: ManifestCtx[AtomiumArgs]) -> list[tuple[str, ...]]:
    # A present input --input-column (ctx.ready) already implies the upstream step
    # succeeded (collect writes a _path only on status OK), so no extra prefilter.
    ready = ctx.ready

    # Run-wide pieces (all structure-independent): argv tail and the bias jsonl.
    bias_path = _write_bias_jsonl(ctx.args.bias_aa, ctx.out_dir)
    atomium_extra = _atomium_extra(ctx.args, bias_path)

    # Guardrails, all cheap and all before submission.
    _check_model_noise(ctx.args.model_noise)
    _check_sampling_temp(ctx.args.sampling_temp)
    _check_symmetry_spec(ctx.args.symmetry)
    # Exactly one source per concept. The column used to win silently over the
    # flag, which is the kind of precedence nobody notices until the designs exist.
    if ctx.args.fixed_positions_column and ctx.args.fixed_positions:
        raise ValueError(
            "--fixed-positions (one spec for the whole run) and "
            "--fixed-positions-column (a per-design spec read from a table column) "
            "are mutually exclusive; pass exactly one. Omit both to redesign every "
            "position, which is the default."
        )
    if ctx.args.chains_to_design_column and ctx.args.chains_to_design:
        raise ValueError(
            "--chains-to-design (one chain list for the whole run) and "
            "--chains-to-design-column (a per-design chain list read from a table "
            "column) are mutually exclusive; pass exactly one."
        )
    if ctx.args.symmetry is not None and ctx.args.tied_positions:
        raise ValueError(
            "--symmetry (homo-oligomer tie) and --tied-positions (explicit tie) are "
            "mutually exclusive"
        )
    if (
        (ctx.args.fixed_positions or ctx.args.fixed_positions_column or ctx.args.tied_positions)
        and not ctx.args.chains_to_design
        and not ctx.args.chains_to_design_column
        and ctx.args.symmetry is None
    ):
        raise ValueError(
            "--fixed-positions/--fixed-positions-column/--tied-positions require "
            "--chains-to-design or --chains-to-design-column (their groups map "
            "one-to-one onto those chains), or --symmetry to derive the chains and "
            "broadcast one group across them"
        )
    # That guardrail tests the FLAGS, which is only a cheap pre-check: what matters
    # is the RESOLVED chain list, which is per-design and therefore re-checked inside
    # the loop below. The run-wide flag's expansion can be checked once, here, so its
    # failure is reported once instead of per design.
    if ctx.args.chains_to_design:
        _explicit_chains(ctx.args.chains_to_design, "--chains-to-design")

    # Pre-validate column specs: one raise naming every offending design, before
    # the design loop and therefore before anything is staged or queued.
    fixed_specs = (
        _column_specs(ready, ctx.args.fixed_positions_column, "--fixed-positions-column")
        if ctx.args.fixed_positions_column
        else {}
    )
    chain_specs = (
        _column_specs(
            ready, ctx.args.chains_to_design_column, "--chains-to-design-column"
        )
        if ctx.args.chains_to_design_column
        else {}
    )

    # Stage each design (CIF->PDB via the shared cache; PDBs returned as-is) and
    # compute its signature. Positions may resolve per-design via lineage, so grouping
    # by signature lets same-param designs share one batched atomium.py call.
    # Per-design failures are COLLECTED, not raised on the first one: on a wide
    # table a single wrong chain letter would otherwise be a resubmit-per-bad-row
    # loop. One raise at the end, naming every offender, matching _column_specs.
    groups: dict[Signature, list[Member]] = defaultdict(list)
    design_errors: list[str] = []
    for design_name in sorted(cast(str, n) for n in ready.index):
        input_path = Path(str(ready.at[design_name, ctx.args.input_column]))
        pdb_src = volume_path(ensure_pdb(input_path, ctx.args.run_dir))
        try:
            # Resolve chains: the column and the flag are mutually exclusive (above).
            if ctx.args.chains_to_design_column:
                chains_source = (
                    f"--chains-to-design-column {ctx.args.chains_to_design_column!r}"
                )
                chains = _explicit_chains(
                    chain_specs[design_name], chains_source, design_name
                )
            else:
                chains = _resolve_chains(ctx.args, pdb_src)
                chains_source = (
                    "--chains-to-design" if ctx.args.chains_to_design else "--symmetry"
                )
            # The designed chains must exist in THIS design's structure, and their
            # lengths are what bounds-check the position specs below. "" here means
            # "design every chain" and is only reachable with no chain flag at all;
            # an explicit source that expands to nothing was refused above.
            chain_lengths: list[int] | None = None
            if chains:
                _validate_chains(chains, pdb_src, design_name, chains_source)
                by_name = _chain_lengths_by_name(pdb_src)
                chain_lengths = [by_name[c] for c in chains.split()]
            # Under --symmetry a lone position group describes one asymmetric unit and
            # is replicated across the designed chains (one group per chain).
            broadcast = len(chains.split()) if ctx.args.symmetry is not None else None
            # Resolve fixed positions: the column and the flag are mutually exclusive.
            fixed_spec = (
                fixed_specs[design_name]
                if ctx.args.fixed_positions_column
                else ctx.args.fixed_positions
            )
            # The guardrail above tested the flags; this tests the RESOLVED list, which
            # is what the helper scripts actually receive. Positions with no designed
            # chain pin nothing and exit 0.
            if (fixed_spec or ctx.args.tied_positions) and not chains:
                raise ValueError(
                    f"design {design_name!r}: positions were given but the resolved "
                    f"designed-chain list is EMPTY ({chains_source}). The helpers "
                    f"read an empty --chain_list as 'every chain' and pin nothing, "
                    f"at exit 0 with status OK."
                )
            fixed_pl = _parse_positions(
                fixed_spec,
                ctx.lookup,
                design_name,
                broadcast_to=broadcast,
                chain_lengths=chain_lengths,
                flag=(
                    f"--fixed-positions-column {ctx.args.fixed_positions_column!r}"
                    if ctx.args.fixed_positions_column
                    else "--fixed-positions"
                ),
            )
            if ctx.args.symmetry is not None:
                tie_mode, tied_pl = "homo", ""
            elif ctx.args.tied_positions:
                tie_mode, tied_pl = (
                    "explicit",
                    _parse_positions(
                        ctx.args.tied_positions,
                        ctx.lookup,
                        design_name,
                        chain_lengths=chain_lengths,
                        flag="--tied-positions",
                    ),
                )
            else:
                tie_mode, tied_pl = "", ""
        except ValueError as exc:
            # Every collected line must name its design: most of these messages
            # already do, but a helper raising from deeper (e.g. a structure with no
            # polymer chain at all) does not, and an unattributable line in a batch
            # of 600 is useless.
            message = str(exc)
            if not message.startswith("design "):
                message = f"design {design_name!r}: {message}"
            design_errors.append(message)
            continue
        sig: Signature = (chains, fixed_pl, tie_mode, tied_pl)
        groups[sig].append((design_name, pdb_src))

    if design_errors:
        shown = "\n".join(f"  - {m}" for m in design_errors[:10])
        more = (
            f"\n  ... and {len(design_errors) - 10} more"
            if len(design_errors) > 10
            else ""
        )
        raise ValueError(
            f"{len(design_errors)} of {len(ready.index)} ready design(s) failed "
            f"validation; NOTHING was submitted and nothing was staged:\n"
            f"{shown}{more}\n"
            f"Fix the upstream rows, or drop them with a -f/--filter module."
        )

    # Split each signature group into subgroups of <= X, then bin-pack subgroups into
    # tasks with a total budget of X designs each.
    x = ctx.args.designs_per_task
    chunks: list[Subgroup] = []
    for sig, members in groups.items():
        for i in range(0, len(members), x):
            chunks.append((sig, members[i : i + x]))
    tasks = _pack_subgroups(chunks, x)

    # Materialize: one grp_<digest>/ output dir per subgroup (inputs/ holds
    # symlinks), one sub-manifest per task (a row per subgroup), and a 1-field
    # top-level row per task. Plan every directory name FIRST and refuse the whole
    # submit before touching the filesystem if a design is already staged elsewhere
    # under this leaf.
    plan: list[list[tuple[Path, Signature, list[Member]]]] = []
    planned_dir_by_design: dict[str, Path] = {}
    for subgroups in tasks:
        task_plan: list[tuple[Path, Signature, list[Member]]] = []
        for sig, members in subgroups:
            digest = _grp_digest(*sig, [m[0] for m in members])
            grp_dir = ctx.out_dir / f"grp_{digest}"
            for name, _ in members:
                planned_dir_by_design[name] = grp_dir
            task_plan.append((grp_dir, sig, members))
        plan.append(task_plan)
    _resolve_restage(ctx.out_dir, planned_dir_by_design, bool(ctx.args.force))

    tasks_dir = ctx.out_dir / "atomium_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[tuple[str, ...]] = []
    for task_plan in plan:
        sub_rows: list[tuple[str, ...]] = []
        for grp_dir, (chains_sig, fixed_pl, tie_mode, tied_pl), members in task_plan:
            _stage_inputs(grp_dir / "inputs", members)
            sub_rows.append(
                (
                    str(grp_dir),
                    chains_sig,
                    fixed_pl,
                    tie_mode,
                    tied_pl,
                    ctx.args.sampling_temp,
                    atomium_extra,
                )
            )
        text = "".join("\t".join(row) + "\n" for row in sub_rows)
        # Content-addressed for the same reason grp_ dirs are: `task_0.tsv` restarted
        # at 0 on every submit, and a task reads its sub-manifest at CONTAINER START,
        # which can be long after the submit returned. A content-addressed file is
        # only ever rewritten with identical bytes, so no in-flight task can have one
        # changed under it.
        task_file = tasks_dir / f"task_{_text_digest(text)}.tsv"
        task_file.write_text(text)
        manifest_rows.append((str(task_file),))

    return manifest_rows


def add_run_atomium_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--model-noise",
        type=str,
        default="n05",
        metavar="N00..N07",
        help="Noise level of the AtomiUM weights to design with, i.e. the sigma the "
        "model was trained to denoise: n00=0.0A, n01=0.1A ... n07=0.7A (plus the "
        "n05_v2 retrain). Higher tolerates rougher backbones; lower trusts the input "
        "geometry more. Default n05. Selects model_weights/model_noised_=<noise>.pt, "
        "so an unknown value is rejected at submit time.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=10,
        help="Max designs per task (default 10). Designs are auto-grouped by identical "
        "params (chains, fixed/tied positions, symmetry) so a group shares one batched "
        "atomium.py call, and groups are bin-packed up to this budget. Lower it for "
        "more parallelism, raise it for fewer tasks.",
    )
    parser.add_argument(
        "--chains-to-design",
        type=str,
        default="",
        help="Chains to design, in ProteinMPNN order, as a chain mini-language: ':' is "
        "an inclusive letter range and ',' separates (e.g. 'A:C,E' -> 'A B C E'). Passed "
        "as --chain_list to assign_fixed_chains and, for the position flags, to "
        "make_{fixed,tied}_positions_dict. Empty (default): design all chains, or the "
        "structure's first N when --symmetry says so.",
    )
    parser.add_argument(
        "--chains-to-design-column",
        type=str,
        default="",
        help="Column to read per-design chain lists from. MUTUALLY EXCLUSIVE with "
        "--chains-to-design (passing both raises). Same chain mini-language. A "
        "chain letter the design's own structure does not have refuses the submit "
        "(it would otherwise design nothing and still collect OK). Missing or "
        "empty values refuse the whole submit.",
    )
    parser.add_argument(
        "--fixed-positions",
        type=str,
        default="",
        help="Positions to keep FIXED (not redesigned), as a position mini-language: "
        "'/' breaks chains (groups map one-to-one onto --chains-to-design, in order), "
        "',' separates fragments within a chain, 'start:end' is an inclusive range; table "
        "column expressions go in {...} islands (resolved up the lineage with + - * // "
        "arithmetic), everything else is a literal integer. 1-indexed within each parsed "
        "chain. E.g. '9:23/10,11,18:20,22'. One group per designed chain is expected, "
        "unless --symmetry broadcasts a single group across them. Empty (default): "
        "redesign everything.",
    )
    parser.add_argument(
        "--fixed-positions-column",
        type=str,
        default="",
        help="Column to read per-design fixed-position specs from. MUTUALLY "
        "EXCLUSIVE with --fixed-positions (passing both raises). The column value "
        "must be in the position mini-language (e.g. '11:18,45:54'), 1-indexed "
        "within the designed chain. Positions are bounds-checked against that "
        "design's own chain length, and the number of '/'-separated groups must "
        "equal the number of designed chains. Missing or empty values refuse the "
        "whole submit.",
    )
    parser.add_argument(
        "--tied-positions",
        type=str,
        default="",
        help="Positions to TIE across chains for symmetric design, same position "
        "mini-language as --fixed-positions (groups map to --chains-to-design; tied "
        "groups must be equal length and are tied index-parallel). E.g. '1:8/1:8'. "
        "Mutually exclusive with --symmetry. Empty (default): no explicit tie.",
    )
    parser.add_argument(
        "--symmetry",
        type=str,
        default=None,
        metavar="AUTO|N",
        help="Homo-oligomer convenience. Ties all chains via make_tied_positions_dict "
        "--homooligomer 1, AND lets --fixed-positions describe a single asymmetric unit: "
        "a lone position group is broadcast across the designed chains. 'auto' takes the "
        "order from the input structure's polymer chain count; otherwise pass that count "
        "as a plain integer. With no --chains-to-design, the designed chains are the "
        "structure's first N. Omitted by default. Mutually exclusive with "
        "--tied-positions.",
    )
    parser.add_argument(
        "--bias-aa",
        type=str,
        default="",
        help="Global AA composition bias as space-separated AA:bias pairs "
        "(e.g. 'D:1.39 E:1.39 K:1.39'), applied via --bias_AA_jsonl. "
        "Empty (default): no bias.",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="FLAG",
        help="Raw flag forwarded verbatim to atomium.py (repeatable); also accepts "
        "pre-made jsonl paths. E.g. --set '--omit_AAs AC', --set '--max_length 5000', "
        "--set '--pssm_jsonl /path/pssm.jsonl'. Tokens are split on whitespace, so a "
        "value containing a space cannot be passed this way.",
    )
    parser.add_argument("--num-seq-per-target", type=int, default=2)
    parser.add_argument(
        "--sampling-temp",
        type=str,
        default="0.1",
        help="One or more space-separated sampling temperatures (e.g. '0.1' or "
        "'0.1 0.2'). AtomiUM generates a full --num-seq-per-target set at EACH "
        "temperature, so N temperatures give N x num_seq sequences per design. "
        "Default '0.1'.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=37,
        help="Random seed (default 37). 0 means pick a random seed per task.",
    )
