#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) asking, per design: DID THIS BINDER LAND ON THE
EPITOPE I CHOSE, AND WHICH TARGET RESIDUES DOES IT ACTUALLY TOUCH?

The decision it exists for: a de-novo binder campaign against a chosen epitope of a
target produces designs that are confident, well folded, and bound somewhere else.
"Somewhere else" is invisible in every confidence score and in every interface-size
metric -- a binder sitting on the opposite face of the target has a perfectly good
if_dG and a perfectly good CMS. This tool turns the epitope question into COLUMNS,
so it can be filtered:

    sapia run proteinmpnn ... -f 'epitope_status == "OK" and epitope_hotspot_recall >= 0.5'
    # and its inverse, to see what ignored the epitope entirely:
    #   -f 'epitope_status == "OK" and epitope_hotspot_recall == 0.0'
    # recall is binary per residue; when nearly every design grazes a hotspot it
    # discriminates nothing, and the WEIGHT is the question instead:
    #   -f 'epitope_his_status == "OK" and epitope_his_hotspot_contact_frac < 0.05'

PREMISE AND SCOPE. Two sets of chains in ONE structure file: ``--design-chains``
(the binder) and ``--target-chains`` (the target), plus a list of target residues
the binder was MEANT to touch (``--hotspots``). It measures heavy-atom CONTACT
between the two sides and nothing else:

  * NOT buried surface area, NOT shape complementarity -- use ``cms``.
  * NOT interface energy, NOT packing -- use ``pyrosetta`` (if_dG, if_dSASA).
  * NOT whether the structure it is scoring is TRUSTWORTHY. It scores the
    coordinates it is given. On a predicted complex it must be read AFTER the gate
    that says the target landed (usalign of the predicted target chains against the
    template) -- a prediction that placed the target wrongly will still produce a
    tidy hotspot_recall, of a fiction.
  * It makes no assembly, oligomer or membrane assumptions -- that is ``ringfit``,
    whose numbers (bridge_ratio, n_clash, lipid_clash) are undefined here. This
    tool is equally happy on a single-chain target, which ringfit is not.

``action: update`` -- where a binder landed is a property of a design that already
exists, so it annotates the table in place and mints no child. Run it twice with
different ``-l/--dir-label`` to score the designed backbone and its prediction onto
one table.

``default_input_column`` is the ``"not applicable"`` sentinel (as in cms, chainsel,
ssprofile): there is no honest default -- ``rfdiffusion3_path`` (designed complex),
``boltz_path`` and ``bindcraft2_path`` (predicted complexes) are all normal inputs.
The builder raises unless ``-i/--input-column`` is given, so it is effectively
required.

THE STRUCTURE COLUMN IS RESOLVED UP THE LINEAGE. The commonest invocation runs
epitope on a CHILD table (the sequence-design generation) with
``-i rfdiffusion3_path``, which lives in the PARENT. ``core/base_run.py`` hands an
update tool a plain ``read_frame(table)`` -- **not** ``join_lineage()`` -- so that
column is simply not in the frame; and ``filter_ready`` is defensive, leaving the
frame unchanged when the input column is missing rather than emptying it. So the
rows arrive "ready" with nothing to read. The builder therefore resolves each path
with ``ctx.lookup`` (row first, then up ``parent_name``/``parent_table``), exactly
as ``usalign`` resolves ``--col-b``. A design whose structure resolves NOWHERE is
**submitted anyway with an empty path** and lands in the table as
``error: could not resolve ...`` -- never dropped from the manifest, because a run
that silently submits fewer designs than the table holds is this workspace's most
expensive quiet failure.

``--sequence-column`` (optional, off by default) takes the design chain's residue
IDENTITIES from a table column instead of from the coordinates, resolved the same
way up the lineage. It exists because the sequence designer writes into a CHILD
table while the backbone that defines the interface GEOMETRY lives in the PARENT,
and that backbone's own residue names are generator artifacts rather than the
designed sequence. With the flag omitted, behaviour is exactly as before.

THE NUMBERING CONTRACT -- the reason this tool is worth its columns. A hotspot
residue number that does not exist in the target chain is an ERROR for that design,
never a recall of 0.0. Generators and predictors renumber chains from 1 (rfd3 does,
routinely), so hotspots copied off the reference structure can name residues that
are not there, or -- worse -- residues that ARE there but are different residues. A
mis-numbered hotspot silently scoring 0.0 is indistinguishable from a genuinely bad
binder, and would make a campaign discard good designs while believing its data.
Two defences, both mandatory: the error above, and ``epitope_hotspot_resnames`` --
the residue NAME found at each listed position, in the order listed, in the table
where it can be eyeballed against the reference (``A96=LYS,A99=GLU,...``).

Chain lists use the prosapia chain mini-language (``,`` separates, ``:`` is an
inclusive letter range): ``B``, ``A,B``, ``A:D``. Hotspots use the workspace's
chain-prefixed residue syntax, the same one bindcraft2 takes:
``A96,A99,A101,A155`` with ranges ``A96-99``, and ``{expr}`` islands resolved per
design against the table (``A{epitope_start}-{epitope_end}``).

Usage:
    sapia run epitope outputs/20260930_egfr_binder \
        --table table1 \
        --input-column boltz_path \
        --design-chains B --target-chains A \
        --hotspots A96,A99,A101,A155 \
        --dir-label pred

    # On the sequence-design generation: geometry from the PARENT's backbone,
    # residue identities from THIS table's designed sequence.
    sapia run epitope outputs/20260930_egfr_binder \
        --table table1 \
        --input-column rfdiffusion3_path \
        --sequence-column atomium_sequence \
        --design-chains B --target-chains A \
        --hotspots A96,A99,A101,A155 \
        --dir-label bb
"""

import math
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import cast

import pandas as pd

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.core.naming import status_column
from prosapia.utils import ensure_pdb, expand_chain_spec, resolve_template

# The sentinel default_input_column (see the module docstring).
NO_DEFAULT_COLUMN = "not applicable"

# 'A96' or 'A96A' (insertion code); ranges are plain integers only.
_SINGLE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z]?)$")
_RANGE = re.compile(r"^([A-Za-z])(\d+)-([A-Za-z]?)(\d+)$")


# The --seq-source value meaning "identities came from the coordinates".
SEQ_SOURCE_STRUCTURE = "structure"


class EpitopeArgs(CommonArgs):
    design_chains: str
    target_chains: str
    hotspots: str
    contact_cutoff: float
    designs_per_task: int
    clash_cutoff: float
    sequence_column: str | None


def parse_hotspots(spec: str) -> list[str]:
    """Chain-prefixed residue spec -> ordered, expanded residue labels.

    ``'A96,A99,A101-103'`` -> ``['A96', 'A99', 'A101', 'A102', 'A103']``. A single
    residue may carry an insertion code (``A96A``); a range may not (its endpoints
    must be plain integers, and it is expanded by number). Order is preserved --
    hotspot_resnames is read positionally against the reference structure, so the
    order the caller wrote is the order the table reports.

    Raises on anything malformed or duplicated: the hotspot list is the definition
    of the question being asked, and a silently dropped or doubled residue changes
    the denominator of hotspot_recall.
    """
    labels: list[str] = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if (m := _SINGLE.match(token)) is not None:
            chain, num, icode = m.group(1).upper(), m.group(2), m.group(3).upper()
            labels.append(f"{chain}{int(num)}{icode}")
            continue
        if (m := _RANGE.match(token)) is not None:
            chain, start, chain2, end = (
                m.group(1).upper(),
                int(m.group(2)),
                m.group(3).upper(),
                int(m.group(4)),
            )
            if chain2 and chain2 != chain:
                raise ValueError(
                    f"hotspot range {token!r} spans two chains ({chain} -> {chain2}); "
                    f"write one range per chain."
                )
            if end < start:
                raise ValueError(f"hotspot range {token!r} ends before it starts.")
            labels.extend(f"{chain}{n}" for n in range(start, end + 1))
            continue
        raise ValueError(
            f"cannot parse hotspot {token!r}: expected a chain-prefixed residue "
            f"('A96', or 'A96A' with an insertion code) or a range within one chain "
            f"('A96-99')."
        )

    if not labels:
        raise ValueError("--hotspots is empty: name the residues the binder is meant to touch.")
    duplicates = sorted({lab for lab in labels if labels.count(lab) > 1})
    if duplicates:
        raise ValueError(
            f"--hotspots lists {','.join(duplicates)} more than once; each residue "
            f"must appear once (it is the denominator of hotspot_recall)."
        )
    return labels


def _clean(value: object) -> str | None:
    """A table cell as a non-empty string, or None for null/blank/absent."""
    if value is None:
        return None
    try:
        if pd.isna(value): # type: ignore
            return None
    except (TypeError, ValueError):  # arrays and other non-scalars
        return None
    text = str(value).strip()
    return text or None


def _resolve(ctx: ManifestCtx[EpitopeArgs], name: str, column: str) -> str | None:
    """The row's own value for ``column``, else the nearest ancestor's, else None.

    Precedence is explicit -- the row first, the lineage only as a fallback -- and
    not left to ``ctx.lookup``'s internals, so the rule holds whatever the lookup is
    wired to.

    The lineage half is required, not a nicety: ``core/base_run.py`` hands an update
    tool a PLAIN ``read_frame(table)`` (no ``join_lineage``), so a child table
    produced by a sequence designer does not carry the parent's
    ``rfdiffusion3_path``. And ``filter_ready`` is defensive -- a missing input
    column leaves the frame unchanged rather than emptying it -- so without this the
    rows arrive "ready" and then have nothing to read. Same precedent as
    ``usalign --col-b``.
    """
    if column in ctx.df.columns and name in ctx.df.index:
        if (own := _clean(ctx.df.at[name, column])) is not None:
            return own
    return _clean(ctx.lookup(name, column))


def add_run_epitope_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--design-chains",
        type=str,
        required=True,
        help="The BINDER chain(s), in the chain mini-language (',' separates, ':' "
        "is an inclusive letter range): 'B', 'A,B', 'A:D'. Contacts are counted "
        "from the heavy atoms of these chains. A chain named here but ABSENT from "
        "a design is an error for that design, with every column NA -- never a "
        "smaller selection and never a recall of 0.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        required=True,
        help="The TARGET chain(s), same mini-language, e.g. 'A'. These define both "
        "the contact partner and the residue universe the epitope is reported over. "
        "A chain named here but absent from a design is an error for that design.",
    )
    parser.add_argument(
        "--hotspots",
        type=str,
        required=True,
        help="The target residues the binder was MEANT to touch -- the epitope you "
        "chose -- chain-prefixed and comma-joined, e.g. 'A96,A99,A101,A155', with "
        "ranges 'A96-99' and {expr} islands resolved per design ('A{ep_start}-"
        "{ep_end}'). Every chain named here must be in --target-chains. A residue "
        "number that does not exist in the target chain is an ERROR for that "
        "design, never a recall of 0.0 -- check epitope_hotspot_resnames whenever a "
        "recall surprises you.",
    )
    parser.add_argument(
        "--contact-cutoff",
        type=float,
        default=5.0,
        help="Heavy-atom distance (A) within which a design atom and a target atom "
        "count as a contact. Default 5.0. This is a contact test, not a buried-"
        "surface calculation: for interface AREA use cms or pyrosetta. Recorded as "
        "epitope_contact_cutoff, because a recall is only comparable against "
        "another recall at the same cutoff.",
    )
    parser.add_argument(
        "--clash-cutoff",
        type=float,
        default=2.5,
        help="Heavy-atom distance (A) below which a design atom and a target atom "
        "OVERLAP (added 2026-10-08; default 2.5 -- the same VALUE as ringfit's "
        "default, but epitope counts atom PAIRS while ringfit_n_clash counts binder "
        "atoms, so the counts are never comparable). Feeds epitope_n_clash, "
        "epitope_n_clash_design_res, epitope_n_clash_bb and epitope_n_clash_bbcb only "
        "and is echoed as "
        "epitope_clash_cutoff; no pre-existing column depends on it. A geometric "
        "overlap test between atoms of different chains: hydrogen-bond and salt-"
        "bridge distances (2.5-3.0 A) are NOT counted. Must be > 0.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=20,
        help="Designs scored per task (default 20). The work is a distance matrix, "
        "well under a second per design, so a task per design would be nearly all "
        "container start-up. Lower it for more parallelism.",
    )
    parser.add_argument(
        "--sequence-column",
        type=str,
        default=None,
        help="Optional table column holding the DESIGN chain's sequence as a plain "
        "one-chain, one-letter string (e.g. 'atomium_sequence', "
        "'proteinmpnn_sequence'). Resolved up the lineage like --input-column. When "
        "given, the design chain's residue IDENTITIES come from it -- mapped "
        "positionally onto the design chains' amino-acid residues in structure "
        "order -- instead of from the structure; the contact GEOMETRY still comes "
        "from the structure. Needed because a sequence designer writes into a CHILD "
        "table while the backbone that defines the geometry lives in the PARENT, and "
        "that backbone's own residue names are generator artifacts. A length "
        "mismatch against the structure, or a '/' in the string, is an ERROR for "
        "that design naming both numbers -- never a truncated or padded mapping. "
        "Omit it (the default) and the tool behaves exactly as before: identities "
        "from the structure. Recorded as epitope_seq_source.",
    )


def build_epitope_manifest(ctx: ManifestCtx[EpitopeArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0

    design = expand_chain_spec(ctx.args.design_chains)
    if not design:
        raise ValueError("--design-chains must name at least one chain (e.g. 'B').")
    if len(set(design)) != len(design):
        raise ValueError(f"--design-chains lists a chain twice: {design}.")

    target = expand_chain_spec(ctx.args.target_chains)
    if not target:
        raise ValueError("--target-chains must name at least one chain (e.g. 'A').")
    if len(set(target)) != len(target):
        raise ValueError(f"--target-chains lists a chain twice: {target}.")
    overlap = sorted(set(design) & set(target))
    if overlap:
        raise ValueError(
            f"chains {overlap} are in both --design-chains and --target-chains; an "
            f"interface needs two disjoint sides (and a chain cannot be its own "
            f"binding partner)."
        )

    if ctx.args.contact_cutoff <= 0:
        raise ValueError("--contact-cutoff must be > 0.")
    clash_cutoff = float(ctx.args.clash_cutoff)
    if not math.isfinite(clash_cutoff) or clash_cutoff <= 0:
        raise ValueError(f"--clash-cutoff must be a positive number, got {clash_cutoff}.")
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")

    # No honest default input column (see module docstring): refuse loudly rather
    # than submit against a column that resolves nowhere. Note the column does NOT
    # have to be in ctx.df -- see _resolve below.
    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"epitope has no default input column: pass -i/--input-column with the "
            f"structure column holding the binder/target COMPLEX. Structure columns "
            f"in table '{ctx.args.table}': {available or '(none)'} -- a PARENT "
            f"table's column (e.g. rfdiffusion3_path) is also valid and is resolved "
            f"up the lineage."
        )

    seq_column = ctx.args.sequence_column
    seq_source = SEQ_SOURCE_STRUCTURE if seq_column is None else seq_column

    target_set = set(target)

    # NOT ctx.ready. ctx.ready runs filter_ready(df, input_column), which DROPS any
    # row whose input cell is null or blank -- so a design whose structure lives on
    # a parent (the column is absent here) or whose cell is empty would vanish from
    # the manifest and the run would silently be smaller than the table. Every row
    # is a candidate instead, and one that resolves nowhere becomes an ERROR ROW in
    # the table rather than an absence. The resume skip below is the only filter,
    # and it is copied verbatim from ManifestCtx.ready so --force behaves the same.
    ready = ctx.df
    status_col = status_column(ctx.out_dir.name)
    if not ctx.args.force and status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    members: list[tuple[str, str, str, str]] = []
    n_unresolved = 0
    for name in ready.index:
        name = cast(str, name)
        # Row first, then up the lineage -- the usalign --col-b precedent. The
        # frame an update tool is handed is a PLAIN read of its own table (no
        # join_lineage), so a binder table produced by a sequence designer does NOT
        # carry the parent's rfdiffusion3_path and this fallback is what makes
        # `-i rfdiffusion3_path -t table1` work at all.
        raw = _resolve(ctx, name, column)

        # {expr} is resolved HERE, at build time, so a bad column or a malformed
        # hotspot kills the submit instead of producing N meaningless rows.
        try:
            hotspots = parse_hotspots(
                resolve_template(ctx.args.hotspots, ctx.lookup, name)
            )
        except ValueError as e:
            raise ValueError(f"--hotspots, for design {name}: {e}") from e
        stray = sorted({lab[0] for lab in hotspots} - target_set)
        if stray:
            raise ValueError(
                f"--hotspots names chain(s) {','.join(stray)} which are not in "
                f"--target-chains ({','.join(target)}); the epitope must live on the "
                f"target. (design {name})"
            )

        sequence = ""
        if seq_column is not None:
            # Same lineage resolution: the sequence normally lives on THIS row while
            # the structure lives on the parent. An unresolvable sequence is left
            # empty and becomes an 'error:' in the table, not a dropped row.
            sequence = (_resolve(ctx, name, seq_column) or "").replace(" ", "").upper()

        if raw is None:
            # Submitted with an EMPTY structure field on purpose: the worker turns it
            # into 'error: could not resolve ...' IN THE TABLE. Dropping the row here
            # would make the run silently smaller than the table -- the failure mode
            # that reads as "No designs to submit." and exits 0.
            n_unresolved += 1
            print(
                f"{name}: could not resolve {column} (row or lineage) -- submitting "
                f"anyway so it lands as an error, not as a missing row"
            )
            members.append((name, "", ",".join(hotspots), sequence))
            continue

        src = Path(raw)
        # CIF (and .cif.gz) -> PDB up front, cached under run_dir/.cif_to_pdb, so
        # the worker sees one format and one chain-naming convention. A path that
        # does not exist is passed through RAW (ensure_pdb would raise here), so the
        # worker records it as 'error: structure missing: ...' rather than the run
        # dropping the design.
        staged = ensure_pdb(src, ctx.args.run_dir) if src.exists() else src
        members.append(
            (name, str(volume_path(staged)), ",".join(hotspots), sequence)
        )

    if n_unresolved:
        print(
            f"epitope: {n_unresolved}/{len(members)} design(s) had no {column} on "
            f"their row or any ancestor; each will collect as 'error: could not "
            f"resolve ...'."
        )
    if members and n_unresolved == len(members):
        raise ValueError(
            f"epitope resolved no structure at all for any of the {len(members)} "
            f"ready designs via column {column!r} (searched table "
            f"'{ctx.args.table}' and its ancestors). Check -i/--input-column."
        )

    # One sub-manifest per task, and a top-level row per task pointing at it. The
    # hotspot list rides per DESIGN (not per task) so {expr} can differ per row,
    # and so the sub-manifest records exactly which residues were scored.
    tasks_dir = ctx.out_dir / "epitope_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    per_task = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), per_task)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, src, hotspots, sequence in members[i : i + per_task]:
                # The sequence rides LAST because it is the only field that may be
                # empty; the worker splits on tabs and pads, so a 3-field row from an
                # older manifest still reads.
                f.write(f"{name}\t{src}\t{hotspots}\t{sequence}\n")
        manifest_rows.append(
            (
                str(volume_path(task_file)),
                ",".join(design),
                ",".join(target),
                # Every field is non-empty, so none of them can vanish from `cut`.
                str(ctx.args.contact_cutoff),
                # Appended 2026-09-30; both always non-empty. epitope.sh cuts them
                # by index, so new fields go on the END, never in the middle.
                seq_source,
                column,
                # Appended 2026-10-08, always non-empty; the END again, so the six
                # fields above keep their indices (epitope.sh defaults a missing
                # 7th field to 2.5, so a 6-field manifest still runs).
                str(clash_cutoff),
            )
        )

    return manifest_rows
