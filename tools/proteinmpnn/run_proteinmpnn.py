"""Fork of the bundled proteinmpnn manifest builder: PER-DESIGN fixed positions.

WHAT THIS FORK ADDS
-------------------
Two optional flags that read a *per-row table column* instead of one run-wide
CLI string:

    --fixed-positions-column  <column>    per-design --fixed-positions spec
    --chains-to-design-column <column>    per-design --chains-to-design spec

...and a fix to a staging bug inherited from upstream that those flags make
near-certain to hit (see SUBGROUP DIRECTORY NAMES below).

Everything else -- the task script, the Modal image, the collector, every
existing flag and the output layout -- is the bundled tool's, reused unchanged
via ``get_builtin("proteinmpnn").with_overrides(...)`` (see ``spec.py``). The
tool keeps the name ``proteinmpnn``, so it SHADOWS the built-in and every column
it writes (``proteinmpnn_sequence``, ``proteinmpnn_path``, ...) is byte-identical
to before. Downstream steps that key off those columns are unaffected.

WHY
---
``-i/--input-column`` already varies per row: the builder does
``ready.at[design_name, ctx.args.input_column]``. ``--fixed-positions`` does
not -- it is one parsed string broadcast to every design. Its only per-row
escape is a ``{column}`` island, and ``prosapia.utils.expr`` restricts those to
INTEGER arithmetic (``_resolve_column`` does ``float(val)`` and raises on
anything non-numeric). A composite multi-fragment spec like
``"11:18,45:54,87:92"`` therefore cannot be expressed as ``{...}`` islands
unless the upstream tool emits six separate integer boundary columns.

``motifrmsd`` emits ONE composite string column
(``motifrmsd_design_fixed_positions``, e.g. ``"11:18,45:54,87:92"``) holding the
grafted motif's segments in the design's own ProteinMPNN numbering, plus the
design-side chain letter in ``motifrmsd_fixed_positions_chain``. This fork is
the missing per-row input mode that lets those columns reach
``--fixed-positions`` / ``--chains-to-design``.

PREMISE, AND WHERE IT IS MEANINGLESS
------------------------------------
The column value is used VERBATIM as that design's spec. It must already be a
complete position mini-language string in **ProteinMPNN position space** --
1-indexed within each chain as ProteinMPNN renumbers it (1..L in file order),
NOT original PDB numbering. Nothing here checks that: a spec computed against
the wrong numbering fixes the wrong residues and the run still exits 0. Whoever
produced the column owns that guarantee (``motifrmsd`` does: it refuses to emit
``design_fixed_positions`` unless the design chain is numbered 1..N in file
order).

The spec's ``/``-separated groups map one-to-one onto the designed chains, in
order, exactly as for ``--fixed-positions``. A single-chain spec (no ``/``) is
therefore only correct when exactly one chain is designed.

REUSE, NOT REIMPLEMENTATION
---------------------------
Every parsing helper is imported from the built-in module and called unchanged:
``_parse_positions`` (which is the thin formatter over
``prosapia.utils.parse_positions`` -> ``parse_chain_positions``),
``_parse_chains``, ``_resolve_chains``, ``_write_bias_jsonl``, ``_mpnn_extra``,
``_check_symmetry_spec``, ``_pack_subgroups``, ``_symlink_member``. No colon /
comma / slash parsing is duplicated here.

A consequence of reusing ``_parse_positions`` unchanged: a column value that
contains ``{...}`` is STILL run through ``resolve_template`` and resolved as a
lineage expression. On a brace-free string (what ``motifrmsd`` emits) that is a
no-op, which is the whole reason the reuse is safe. A brace that is not a valid
integer column expression raises loudly; it never resolves to a wrong number
silently.

SUBGROUP DIRECTORY NAMES: THE BUG THIS FORK ALSO FIXES
------------------------------------------------------
UPSTREAM-INHERITED, not introduced here. Upstream names each staged subgroup
with a counter that restarts at 0 on every submit::

    gid = 0
    ...
    grp_dir = ctx.out_dir / f"grp_{gid}"

Two separate ``sapia run`` calls against the same table and the same leaf (one
``-l/--dir-label``, the shard pattern a wide campaign needs and the pattern
``running-a-step`` prescribes: submit, collect, submit) therefore BOTH create a
``grp_0``. It is the same physical directory, and upstream's ``_symlink_member``
appends without clearing -- so submit 2's ``grp_0/inputs/`` ends up holding
submit 1's PDB as well, and the single ``protein_mpnn_run`` call for that
subgroup processes BOTH under submit 2's spec.

Measured on Modal, 2026-10-09::

    File ".../protein_mpnn_utils.py", line 310, in tied_featurize
        fixed_position_mask[np.array(fixed_pos_list)-1] = 0.0
    IndexError: index 105 is out of bounds for axis 0 with size 105

The crash was the lucky case: it required the two designs' chain lengths to
differ. With equal lengths the same collision silently redesigns one design's
backbone under the OTHER design's fixed-positions spec, exit 0, status OK.

Upstream rarely notices because one run-wide spec means ONE signature and
usually one subgroup; this fork's per-design specs make every design its own
subgroup, so the very first subgroup of any second submit collides.

The fix, in the marked region at the end of the builder: subgroup directories
and task sub-manifests are CONTENT-ADDRESSED (``grp_<digest>`` over the
signature plus the member design names; ``task_<digest>.tsv`` over the file's
own bytes) instead of counter-numbered, staging makes a directory match its
member list exactly rather than appending to it, and a design found already
staged under the same leaf in a DIFFERENT subgroup refuses the submit (or, with
``--force`` and when that whole subgroup is being redone, replaces it). Nothing
downstream depends on the old names: ``proteinmpnn.sh`` reads the grp dir from
field 1 of its sub-manifest, and the bundled collector globs ``*/seqs/*.fa``
under the out_dir and keys rows by FASTA stem, never by ``grp_<n>``.

Residual risk this does NOT fix, because it lives in prosapia, not here: the
top-level manifest (``.manifests/<table>_<leaf>_manifest.txt``) and the task
logs (``<script>_<id>.out/.err/.exit``, ``<script>_modal.json``) are still
written per table+leaf with task ids restarting at 1, so a submit launched while
an earlier submit's tasks are still queued or running will clobber them. The
run-wide ``bias_AA.jsonl`` (upstream's, written once per submit into out_dir) is
overwritten the same way. All three are safe only sequentially: submit, wait,
collect, then submit the next shard -- which is already the rule.

THE DIFF AGAINST UPSTREAM
-------------------------
``build_proteinmpnn_manifest`` below is a copy of
``prosapia.tools.proteinmpnn.run_proteinmpnn.build_proteinmpnn_manifest`` with
three marked regions added (the guards, the per-row lookup inside the design
loop, and the content-addressed staging that replaces upstream's ``gid``
counter). Everything outside those markers is upstream's code, kept so the fork
is easy to re-diff when prosapia moves.
"""

import hashlib
import shutil
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
from typing import cast

import pandas as pd

from prosapia.core import ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.tools.proteinmpnn.run_proteinmpnn import (
    Member,
    ProteinMPNNArgs,
    Signature,
    Subgroup,
    _check_symmetry_spec,
    _mpnn_extra,
    _pack_subgroups,
    _parse_chains,
    _parse_positions,
    _resolve_chains,
    _symlink_member,
    _write_bias_jsonl,
    add_run_proteinmpnn_args as _add_builtin_run_args,
)
from prosapia.utils import ensure_pdb


class ProteinMPNNColumnArgs(ProteinMPNNArgs):
    """Upstream's args plus the two per-row column flags this fork adds."""

    fixed_positions_column: str
    chains_to_design_column: str


def _column_specs(
    ready: pd.DataFrame, column: str, flag: str
) -> dict[str, str]:
    """Per-design spec strings read verbatim from ``column``, for every ready design.

    Validates the WHOLE batch in one pass and raises ONCE, naming every offending
    design, so a caller repairs the table in a single edit instead of discovering
    the bad rows one resubmit at a time.

    A missing value is **refused, never defaulted**. ``NaN``, ``None``, the empty
    string and a whitespace-only string all count as missing. This matters more
    than it looks: an empty spec would sail through ``_parse_positions`` as ``""``,
    produce no ``--fixed_positions_jsonl``, and silently REDESIGN THE WHOLE CHAIN
    at exit 0 -- the exact failure a per-design motif graft cannot survive.

    No per-design isolation is possible here. ``build_manifest_fn`` runs once, in
    the submit process, before anything is queued (``core/base_run.py``), so one
    bad row refuses the whole submit. That is deliberate -- a half-submitted batch
    is worse -- and the remedy is a ``-f/--filter`` module, which the error names.
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
        value = ready.at[name, column]
        if isinstance(value, str):
            text = value.strip()
        elif value is None or pd.isna(value):
            text = ""
        else:
            text = str(value).strip()
        if not text:
            missing.append(name)
        else:
            specs[name] = text

    if missing:
        shown = ", ".join(missing[:10])
        more = f" (and {len(missing) - 10} more)" if len(missing) > 10 else ""
        raise ValueError(
            f"{flag} {column!r}: {len(missing)} of {len(ready.index)} ready "
            f"design(s) have no value in that column: {shown}{more}. "
            f"Refusing the submit rather than redesigning those chains in full "
            f"(an empty spec fixes nothing and exits 0). Either fix the upstream "
            f"tool's rows, or drop them with a -f/--filter module such as "
            f"'df[df[\"{column}\"].notna() & (df[\"{column}\"] != \"\")]'."
        )
    return specs


_DIGEST_LEN = 16


def _digest(*parts: str) -> str:
    """A short, stable content digest (64 bits of sha256, hex)."""
    h = hashlib.sha256()
    for part in parts:
        h.update(part.encode())
        h.update(b"\x1f")  # unambiguous field separator
    return h.hexdigest()[:_DIGEST_LEN]


def _subgroup_dir(out_dir: Path, sig: Signature, members: list[Member]) -> Path:
    """The staging dir for one subgroup, named by WHAT IS IN IT, not by a counter.

    FORK DIFF, and the reason it exists. Upstream names subgroups with a counter
    that restarts at 0 on every submit (``gid = 0`` ... ``grp_{gid}``). Two
    separate ``sapia run`` calls against the same table+leaf -- the shard pattern
    the campaign needs, and the one ``running-a-step`` prescribes (submit,
    collect, submit) -- therefore BOTH write a ``grp_0``, into the same physical
    directory, and ``_symlink_member`` appends rather than clears. The second
    submit's single ``protein_mpnn_run`` call then processes the first submit's
    PDB too, under the second submit's spec. Measured 2026-10-09:
    ``IndexError: index 105 is out of bounds for axis 0 with size 105`` in
    ``tied_featurize``. That crash was luck -- it needed the two chain lengths to
    differ. At equal lengths the same collision applies the WRONG design's
    fixed-positions spec to the other design's backbone and exits 0.

    The name is therefore derived from the subgroup's content: its signature plus
    its member design NAMES. Properties this buys, each of which the counter
    lacked:

      * two submits can never collide unless they stage literally the same
        designs with the same signature -- in which case sharing the directory is
        correct, not a collision;
      * a ``--force`` rerun of the same designs reuses the same directory and
        overwrites its outputs, instead of accumulating a duplicate subgroup and
        leaving two FASTAs per design for the collector to pick between;
      * no scan of "what is already on disk" is needed to pick the name, so it is
        not racy against a concurrent submit.

    Member PATHS are deliberately NOT in the digest: a design whose input path
    changed must reuse its directory (``_stage_inputs`` repoints the symlink), or
    its stale FASTA would be left behind in an orphan directory for the collector
    to find. For the same reason the run-wide ``mpnn_extra`` is not in the digest
    either -- a rerun with different sampling params must OVERWRITE the old
    outputs, not sit beside them.
    """
    chains, fixed_pl, tie_mode, tied_pl = sig
    names = sorted(name for name, _ in members)
    return out_dir / f"grp_{_digest(chains, fixed_pl, tie_mode, tied_pl, *names)}"


def _resolve_restage(
    out_dir: Path, planned_dir_by_design: dict[str, Path], force: bool
) -> None:
    """Refuse (or, under ``--force``, clean up) a design staged elsewhere here.

    With content-addressed names a design lands in a DIFFERENT directory than
    last time only if its signature or its subgroup membership changed since an
    earlier submit to this same table+leaf. The earlier directory keeps its
    outputs, and the bundled collector globs every ``*/seqs/*.fa`` under
    ``out_dir`` keying rows by FASTA stem -- so two directories holding that
    design means two FASTAs and the collector keeps whichever it read last. A
    silent wrong answer, so it is handled here, at submit time.

    Two cases, deliberately different:

      * ``--force`` and every design in the stale directory is being restaged by
        THIS submit -> the old outputs are wholly superseded, so the stale
        directory is removed (loudly) and the run proceeds. This is the ordinary
        "re-run a superset of what I ran before" case.
      * anything else (no ``--force``, or the stale directory also holds designs
        this submit is not redoing) -> raise, naming every design and directory.
        Deleting would destroy a result the caller did not ask to redo.
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
            else "Those dirs also hold designs this submit is not redoing, so "
            "they cannot be replaced wholesale"
        )
        raise ValueError(
            f"{len(blocking)} design(s) are already staged under this output "
            f"dir ({out_dir}) in a DIFFERENT subgroup than this submit would "
            f"use: {shown}{more}. The same design was submitted to this "
            f"table+leaf before with a different chains/fixed-positions/tied "
            f"signature; keeping both would leave two FASTAs for one design and "
            f"the collector would pick one arbitrarily. {hint}, or use a fresh "
            f"-l/--dir-label for this variant, or delete the stale grp_ dir(s) "
            f"named above."
        )

    for old in sorted(removable):
        print(
            f"proteinmpnn: --force, removing superseded subgroup {old} "
            f"({len(members_by_dir[old])} design(s), all restaged by this submit)"
        )
        shutil.rmtree(old)


def _stage_inputs(inputs_dir: Path, members: list[Member]) -> None:
    """Stage exactly ``members`` into ``inputs_dir`` -- no more, no less.

    FORK DIFF. Upstream's ``_symlink_member`` (still the thing that creates each
    link, reused unchanged) only ever APPENDS: it skips a link that exists and
    never removes one. That is what turned a directory-name collision into two
    designs in one ``protein_mpnn_run`` call. Here the directory is made to match
    the member list exactly:

      * an entry that is not an expected member is removed (with
        content-addressed names this cannot arise from a legitimate submit; it is
        an invariant guard against a hand-edited or interrupted directory, and it
        prints when it fires);
      * a symlink pointing somewhere other than this member's staged PDB is
        repointed, so a re-run whose input path moved does not silently redesign
        the old structure.
    """
    inputs_dir.mkdir(parents=True, exist_ok=True)
    expected = {f"{name}.pdb": Path(pdb_src) for name, pdb_src in members}

    for entry in sorted(inputs_dir.iterdir()):
        if entry.name not in expected and (entry.is_symlink() or entry.is_file()):
            print(
                f"proteinmpnn: removing stale staged input {entry} "
                f"(not a member of this subgroup)"
            )
            entry.unlink()

    for name, pdb_src in members:
        link = inputs_dir / f"{name}.pdb"
        if link.is_symlink():
            if link.readlink() == Path(pdb_src):
                continue
            print(
                f"proteinmpnn: repointing staged input {link} "
                f"({link.readlink()} -> {pdb_src})"
            )
            link.unlink()
        elif link.exists():
            link.unlink()
        _symlink_member(Path(pdb_src), inputs_dir, name)


def build_proteinmpnn_manifest(
    ctx: ManifestCtx[ProteinMPNNColumnArgs],
) -> list[tuple[str, ...]]:
    # A present input --input-column (ctx.ready) already implies the upstream
    # diffusion succeeded (collect writes a _path only on status OK), so no extra
    # input-status prefilter is needed here.
    ready = ctx.ready

    # Run-wide pieces (all structure-independent): argv tail and the bias jsonl.
    bias_path = _write_bias_jsonl(ctx.args.bias_aa, ctx.out_dir)
    mpnn_extra = _mpnn_extra(ctx.args, bias_path)

    # --- FORK DIFF START: per-row column flags -------------------------------
    fixed_col = (getattr(ctx.args, "fixed_positions_column", "") or "").strip()
    chains_col = (getattr(ctx.args, "chains_to_design_column", "") or "").strip()

    # Exactly one source per concept. Silent precedence between a run-wide string
    # and a per-row column is the kind of thing nobody notices until the designs
    # are already made.
    if fixed_col and ctx.args.fixed_positions:
        raise ValueError(
            "--fixed-positions (one spec for the whole run) and "
            "--fixed-positions-column (a per-design spec read from a table "
            "column) are mutually exclusive; pass exactly one. Omit both to "
            "redesign every position, which is the default."
        )
    if chains_col and ctx.args.chains_to_design:
        raise ValueError(
            "--chains-to-design (one chain list for the whole run) and "
            "--chains-to-design-column (a per-design chain list read from a "
            "table column) are mutually exclusive; pass exactly one."
        )
    # --- FORK DIFF END --------------------------------------------------------

    # Guardrails: tying is either the homo-oligomer shortcut or explicit, not both;
    # per-chain positions need a chain list to map their groups onto, which
    # --symmetry can supply from the structure instead.
    _check_symmetry_spec(ctx.args.symmetry)
    if ctx.args.symmetry is not None and ctx.args.tied_positions:
        raise ValueError(
            "--symmetry (homo-oligomer tie) and --tied-positions (explicit tie) are "
            "mutually exclusive"
        )
    if (
        # FORK DIFF: `or fixed_col` -- a per-row spec needs the chain list just as
        # much as a run-wide one, and `or chains_col` -- a per-row chain list
        # satisfies that need. Without this the new flag would bypass the
        # guardrail and ProteinMPNN would silently redesign everything.
        (ctx.args.fixed_positions or fixed_col or ctx.args.tied_positions)
        and not ctx.args.chains_to_design
        and not chains_col
        and ctx.args.symmetry is None
    ):
        raise ValueError(
            "--fixed-positions/--fixed-positions-column/--tied-positions require "
            "--chains-to-design or --chains-to-design-column (their groups map "
            "one-to-one onto those chains), or --symmetry to derive the chains "
            "and broadcast one group across them"
        )

    # --- FORK DIFF START: validate the whole batch before the design loop ------
    # One raise naming every offending design, rather than aborting on the first.
    fixed_specs = _column_specs(ready, fixed_col, "--fixed-positions-column") if fixed_col else {}
    chain_specs = (
        _column_specs(ready, chains_col, "--chains-to-design-column") if chains_col else {}
    )
    # --- FORK DIFF END --------------------------------------------------------

    # Stage each design (CIF->PDB via the shared cache; PDBs returned as-is) and
    # compute its signature. Positions may resolve per-design via lineage, so grouping
    # by signature lets same-param designs share one batched protein_mpnn_run call.
    groups: dict[Signature, list[Member]] = defaultdict(list)
    for design_name in sorted(cast(str, n) for n in ready.index):
        input_path = Path(str(ready.at[design_name, ctx.args.input_column]))
        pdb_src = volume_path(ensure_pdb(input_path, ctx.args.run_dir))
        # FORK DIFF: a per-row chain list, read with the same `ready.at[...]`
        # pattern -i already uses; otherwise upstream's flag/symmetry resolution.
        if chains_col:
            chains = _parse_chains(chain_specs[design_name])
        else:
            chains = _resolve_chains(ctx.args, pdb_src)
        # Under --symmetry a lone position group describes one asymmetric unit and
        # is replicated across the designed chains (ProteinMPNN wants one group per
        # chain in --chain_list).
        broadcast = len(chains.split()) if ctx.args.symmetry is not None else None
        # FORK DIFF: the row's own spec, fed STRAIGHT into the unchanged
        # _parse_positions. Brace-free values make resolve_template a no-op, so
        # this is the same code path a literal --fixed-positions takes.
        fixed_spec = (
            fixed_specs[design_name] if fixed_col else ctx.args.fixed_positions
        )
        fixed_pl = _parse_positions(
            fixed_spec,
            ctx.lookup,
            design_name,
            broadcast_to=broadcast,
        )
        if ctx.args.symmetry is not None:
            tie_mode, tied_pl = "homo", ""
        elif ctx.args.tied_positions:
            tie_mode, tied_pl = (
                "explicit",
                _parse_positions(ctx.args.tied_positions, ctx.lookup, design_name),
            )
        else:
            tie_mode, tied_pl = "", ""
        sig: Signature = (chains, fixed_pl, tie_mode, tied_pl)
        groups[sig].append((design_name, pdb_src))

    # Split each signature group into subgroups of <= X, then bin-pack subgroups
    # into array tasks with a total budget of X designs each.
    x = ctx.args.designs_per_task
    chunks: list[Subgroup] = []
    for sig, members in groups.items():
        for i in range(0, len(members), x):
            chunks.append((sig, members[i : i + x]))
    tasks = _pack_subgroups(chunks, x)

    # --- FORK DIFF START: collision-free, content-addressed staging ----------
    # Upstream materialized here with a counter (`gid = 0`; `grp_{gid}`;
    # `task_{t}.tsv`) that restarts at 0 on EVERY submit, so a second submit to
    # the same table+leaf reused the first one's directories and appended into
    # them. Names are now derived from content: see _subgroup_dir. Everything
    # else -- the sub-manifest's six fields and their order, one task file per
    # task, a 1-field top-level row per task -- is unchanged, so proteinmpnn.sh
    # and the bundled collector read exactly what they read before.

    # Plan every directory name first, and refuse the whole submit before
    # touching the filesystem if a design is already staged elsewhere here.
    plan: list[list[tuple[Path, Signature, list[Member]]]] = []
    planned_dir_by_design: dict[str, Path] = {}
    for subgroups in tasks:
        task_plan: list[tuple[Path, Signature, list[Member]]] = []
        for sig, members in subgroups:
            grp_dir = _subgroup_dir(ctx.out_dir, sig, members)
            for name, _ in members:
                planned_dir_by_design[name] = grp_dir
            task_plan.append((grp_dir, sig, members))
        plan.append(task_plan)
    _resolve_restage(ctx.out_dir, planned_dir_by_design, bool(ctx.args.force))

    tasks_dir = ctx.out_dir / "proteinmpnn_tasks"
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
                    mpnn_extra,
                )
            )
        text = "".join("\t".join(row) + "\n" for row in sub_rows)
        # Content-addressed too, and for the same reason: `task_0.tsv` collided
        # across submits exactly as `grp_0` did, and a task file is read by its
        # task at CONTAINER START, which can be long after the submit returned.
        # A content-addressed file is only ever rewritten with identical bytes,
        # so no in-flight task can have one changed under it.
        task_file = tasks_dir / f"task_{_digest(text)}.tsv"
        task_file.write_text(text)
        manifest_rows.append((str(task_file),))
    # --- FORK DIFF END --------------------------------------------------------

    return manifest_rows


def add_run_proteinmpnn_args(parser: ArgumentParser) -> None:
    """Upstream's run flags, plus this fork's two per-row column flags.

    Upstream's ``add_run_proteinmpnn_args`` is called first and unmodified, so
    every existing flag keeps its exact name, default and help text.
    """
    _add_builtin_run_args(parser)
    parser.add_argument(
        "--fixed-positions-column",
        type=str,
        default="",
        help="Read each design's --fixed-positions spec from this TABLE COLUMN "
        "instead of a single run-wide string. The value is used VERBATIM as that "
        "design's spec, so it must already be a complete position mini-language "
        "string in ProteinMPNN numbering (1-indexed per chain, 1..L in file "
        "order) -- e.g. '11:18,45:54,87:92', as produced by "
        "motifrmsd_design_fixed_positions. Mutually exclusive with "
        "--fixed-positions; still requires --chains-to-design(-column) or "
        "--symmetry. A design whose value is empty or NA REFUSES THE SUBMIT "
        "(it is never treated as 'redesign everything'). Empty (default): use "
        "--fixed-positions.",
    )
    parser.add_argument(
        "--chains-to-design-column",
        type=str,
        default="",
        help="Read each design's --chains-to-design list from this TABLE COLUMN "
        "instead of a single run-wide string (e.g. "
        "motifrmsd_fixed_positions_chain, which names the chain its "
        "design_fixed_positions spec addresses). Same chain mini-language as "
        "--chains-to-design. Mutually exclusive with --chains-to-design. A "
        "design whose value is empty or NA refuses the submit. Empty (default): "
        "use --chains-to-design / --symmetry.",
    )
