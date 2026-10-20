#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) asking, per design: does one of the designed
histidines hydrogen-bond the target in the predicted complex?

Flags, scope limits, traps and how to read the columns:
``.claude/skills/hisbond/SKILL.md``.

``action: update``: whether a His bonds the target is a property of an existing design,
so the table is annotated in place and no child is minted. Run it twice with different
``-l/--dir-label`` to score two structure columns onto one table. CPU only (the builder
forces ``gpus_per_task = 0``).

The target side is written merged and per nitrogen (``..._nd1`` / ``..._ne2``), with the
same window and cone and no flag; an "ND1 contact" is ND1 inside the window and the cone
and nothing more. Each qualifying contact also gets a donor/acceptor role, read from the
tautomer when the file has hydrogens (relaxed PyRosetta PDB) and from the partner's
chemistry when it does not (AF3/Boltz CIF); ``role_basis`` says which. The pH-sensitive
contact is the one where the His N ACCEPTS, which ``any_switch_candidate`` gates on. The
same roles are read for INTRA-BINDER contacts (a core His breaking a bond inside its own
fold), with a per-centre verdict across both sides. No role changes a contact
definition.

``default_input_column`` is the ``"not applicable"`` sentinel (as in epitope, cms,
chainsel), so the builder raises unless ``-i/--input-column`` is given.

THE CENTRES COME IN AS PLAIN NUMBERS. ``--centers`` is a comma-separated list of
1-BASED positions along the binder chain in file order, with no default and no
knowledge of any other tool's columns. Per-design values reach it through prosapia's
``{expr}`` island layer, written by the caller::

    --centers '{atomiumph_centre1_idx1},{atomiumph_centre2_idx1}'

One island resolves to exactly one integer and the commas between them are literal
template text (``prosapia/utils/expr.py``), so a comma-joined cell cannot ride through
a single island and the caller names one column per centre. Resolution happens at
manifest-build time, once per design, so a row whose column is NA (fewer centres than
the template names) becomes an ERROR ROW rather than a silently shorter centre list.

THE POSITIONS ARE NEVER INHERITED FROM AN ANCESTOR ROW. ``--centers`` does NOT use
``ctx.lookup``: prosapia's lineage lookup returns the first NON-NULL value found
anywhere up the chain, which on a lineage carrying the column at two levels would hand
a design another row's histidine with status OK. ``_ColumnAwareLookup`` stops at the
first frame whose COLUMNS hold the name and returns that frame's cell as is, NaN
included, so a blank own cell raises and the design becomes a loud error row.
``centers_lineage_depth`` publishes where each design's positions came from, on every
row that carries a position list, so an already-collected table can still be audited.
It is the depth of the DEEPEST island: 0 is exact and means every island resolved from
the design's own table, but a template naming two producers can have islands at
different depths and only the deepest is reported. NA when no position resolved.

The structure column is resolved up the lineage, as in epitope. If this table HAS the
column, the row's own cell is the answer and a blank never falls back to an ancestor
(that would measure another design's structure); only an ABSENT column walks the
lineage. A design that resolves nothing is SUBMITTED ANYWAY and lands as ``error: could
not resolve ...``, never dropped: silently submitting fewer designs than the table holds
is this workspace's most expensive quiet failure.
"""

import math
from argparse import ArgumentParser
from pathlib import Path
from typing import cast

import pandas as pd

from prosapia.core import CommonArgs, DataManager, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.core.naming import PARENT_NAME, PARENT_TABLE, status_column
from prosapia.utils import expand_chain_spec
from prosapia.utils.expr import resolve_template

# The sentinel default_input_column (see the module docstring).
NO_DEFAULT_COLUMN = "not applicable"

# --expect-binder-len-col resolved nowhere: the worker must error, not skip the check.
UNRESOLVED_EXPECT = "__unresolved__"


class HisbondArgs(CommonArgs):
    binder_chain: str
    target_chains: str
    centers: str
    expect_binder_len_col: str | None
    min_dist: float
    max_dist: float
    max_angle: float
    env_radius: float
    designs_per_task: int


def _clean(value: object) -> str | None:
    """A table cell as a non-empty string, or None for null/blank/absent."""
    if value is None:
        return None
    try:
        if pd.isna(value):  # type: ignore[arg-type]
            return None
    except (TypeError, ValueError):  # arrays and other non-scalars
        return None
    text = str(value).strip()
    return text or None


def _resolve(ctx: ManifestCtx[HisbondArgs], name: str, column: str) -> str | None:
    """The row's value for ``column``, or None when it cannot be resolved.

    If this table HAS ``column`` the row's own cell is the answer; a blank never falls
    back to an ancestor (that would describe another design). Only an ABSENT column
    walks the lineage, nearest ancestor first, as in usalign's ``--col-b``.
    """
    if column in ctx.df.columns:
        if name not in ctx.df.index:
            return None
        return _clean(ctx.df.at[name, column])
    return _clean(ctx.lookup(name, column))


class _ColumnAwareLookup:
    """A ``resolve_template`` lookup that NEVER inherits a blank own cell.

    ``DataManager.lookup`` walks ``parent_table``/``parent_name`` and returns the FIRST
    NON-NULL value found anywhere up the chain, so a design whose own
    ``..._centre2_idx1`` is NA silently resolves to an ANCESTOR ROW'S position whenever
    the column exists at two levels (a second upstream leaf, an unlabelled re-run, a
    merge that lifts the column onto more than one generation). The row then collects
    ``OK`` with ``centers_his_ok == centers_n`` while describing a histidine this design
    never had, and the His-identity check only catches it when the inherited position
    happens not to be a His -- which, in a table of pH-switch binders, is the minority.

    This walk stops at the first frame whose COLUMNS contain the name and returns that
    frame's cell AS IS, NaN included, so ``expr.py`` raises and the design becomes an
    error row. ``depth`` is 0 for the row's own table, 1 for its parent, and so on. The
    caller reads ``max_depth`` after each design -- the DEEPEST island, not every
    island's -- and publishes it as ``centers_lineage_depth`` so an already-collected
    table can still be audited. 0 is exact (every island came from the design's own
    table); a non-zero value describes only the deepest of a mixed-producer template.

    Frames are cached for the whole build: one read per ancestor table, not per design.
    """

    def __init__(self, df: pd.DataFrame, run_dir: Path) -> None:
        self._df = df
        # Read-only: DataManager.__init__ only builds the backend closures and
        # __exit__ is scope-only, so a second handle on the run_dir takes no lock.
        self._read_frame = DataManager(Path(run_dir)).read_frame
        self._frames: dict[str, pd.DataFrame] = {}
        self.max_depth = 0

    def _frame(self, table: str) -> pd.DataFrame:
        if table not in self._frames:
            self._frames[table] = self._read_frame(table)
        return self._frames[table]

    def reset(self) -> None:
        self.max_depth = 0

    def __call__(self, name: str, column: str):
        cur_df, cur_name, depth = self._df, name, 0
        seen: set[tuple[str, str]] = set()
        while True:
            if column in cur_df.columns:
                # Found the column: this frame's cell is the answer, NaN included.
                self.max_depth = max(self.max_depth, depth)
                if cur_name not in cur_df.index:
                    return pd.NA
                value = cur_df.at[cur_name, column]
                if isinstance(value, pd.Series):
                    value = value.iloc[0]
                return value
            if cur_name not in cur_df.index:
                return pd.NA
            row = cur_df.loc[cur_name]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            parent_table = row.get(PARENT_TABLE) if PARENT_TABLE in row.index else None
            parent_name = row.get(PARENT_NAME) if PARENT_NAME in row.index else None
            if (
                parent_table is None
                or parent_name is None
                or pd.isna(parent_table)
                or pd.isna(parent_name)
            ):
                return pd.NA
            key = (str(parent_table), str(parent_name))
            if key in seen:  # cycle guard, as in DataManager.lookup
                return pd.NA
            seen.add(key)
            cur_df, cur_name, depth = self._frame(str(parent_table)), str(parent_name), depth + 1


def parse_centres(text: str) -> list[int]:
    """``'12,45'`` -> ``[12, 45]``: 1-BASED positions along the binder chain.

    Kept in sync BY HAND with ``hisbond_worker.py`` (no shared imports): the builder
    validates so a typo fails before anything is submitted, and the worker validates
    again because it is what reads the task file. Raises on anything malformed -- the
    position list defines the question, so guessing at it is never right.
    """
    positions: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            raise ValueError(
                f"empty position in centre list {text!r}; expected comma-separated "
                f"1-based positions, e.g. '12,45'"
            )
        try:
            value = int(part)
        except ValueError:
            raise ValueError(
                f"centre position {part!r} in {text!r} is not an integer; expected "
                f"comma-separated 1-based positions, e.g. '12,45'"
            ) from None
        if value < 1:
            raise ValueError(
                f"centre position {value} in {text!r} is not >= 1; positions are "
                f"1-BASED along the binder chain, so 0 is never one of them"
            )
        positions.append(value)
    if not positions:
        raise ValueError("the centre list is empty: there is nothing to measure")
    if len(set(positions)) != len(positions):
        dup = sorted({p for p in positions if positions.count(p) > 1})
        raise ValueError(
            f"centre position(s) {dup} appear twice in {text!r}; one residue measured "
            f"twice would be counted twice in every n_centres_* column"
        )
    return positions


def add_run_hisbond_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--binder-chain",
        type=str,
        required=True,
        help="The BINDER chain, exactly one chain ID (e.g. 'D'). NO DEFAULT: the "
        "binder is the last chain here, and which letter that is depends on the "
        "target, so a defaulted letter would measure the wrong molecule with "
        "status OK. The centre indices are counted along THIS chain, in file order. "
        "Check hisbond_binder_chain_n_res against the binder length you expect, or "
        "pass --expect-binder-len-col to make that an error instead of a column.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        required=True,
        help="The TARGET chain(s) the His must reach, in the chain mini-language "
        "(',' separates, ':' is an inclusive letter range): 'A,B,C', 'A:C'. NO "
        "DEFAULT. Must be disjoint from --binder-chain (checked at submit). A "
        "candidate partner is every N/O heavy atom of these chains IN A RESIDUE THAT "
        "IS NOT A TABULATED NON-AMINO-ACID -- waters, ions and tabulated ligands are "
        "excluded, the same filter the binder chain gets, and an UNTABULATED residue "
        "name passes. The atoms that filter removed are reported as "
        "hisbond_target_n_het_atoms, so an ignored ligand is visible rather than "
        "silent. A chain named here but absent from a design is an ERROR for that "
        "design, with every science column NA -- never a smaller target and never a "
        "count of 0.",
    )
    parser.add_argument(
        "--centers",
        type=str,
        required=True,
        help="The centres to measure: comma-separated 1-BASED POSITIONS ALONG THE "
        "BINDER CHAIN IN FILE ORDER (e.g. '12,45'), NOT residue numbers. NO DEFAULT "
        "-- there is no honest one, exactly as for --binder-chain and "
        "--target-chains. Per-design values come through prosapia's {expr} islands, "
        "which the CALLER writes: --centers "
        "'{atomiumph_centre1_idx1},{atomiumph_centre2_idx1}' resolves one integer per "
        "island per design (the commas are literal text), so one column is named per "
        "centre and this tool knows no other tool's column names. A design whose "
        "island resolves to nothing or to a non-integer is an ERROR ROW, never a "
        "shorter centre list. A POSITION IS NEVER INHERITED FROM AN ANCESTOR ROW: "
        "unlike prosapia's usual lineage lookup, the resolution stops at the first "
        "table that HAS the column, so a design whose own cell is blank errors "
        "instead of silently borrowing another design's histidine. Where the "
        "positions did come from is published as hisbond_centers_lineage_depth on "
        "every row that carries a position list (0 = this design's own table; with "
        "several islands it is the DEEPEST one). Positions must be >= 1 and distinct. THE RESIDUE FOUND "
        "AT EACH POSITION MUST BE A HISTIDINE or the row errors: that identity check "
        "is why positions are given rather than histidines auto-discovered, and it is "
        "what makes an acid_p/all_p table (Asp/Glu centres) error on every row.",
    )
    parser.add_argument(
        "--expect-binder-len-col",
        type=str,
        default=None,
        help="Optional table column holding the EXPECTED binder length -- an integer "
        "(a sequence-length column) or a plain one-letter sequence (a sequence "
        "column), resolved up the lineage. Default: none, i.e. the check is not made; "
        "no column name is assumed. When given, a binder chain whose residue count "
        "differs is an ERROR for that design: that is the check that catches a chain "
        "letter which EXISTS but is the wrong molecule, which nothing else here can "
        "see. Pass it whenever such a column exists -- the centre positions are "
        "counted along this chain, so a wrong chain puts every one of them on the "
        "wrong residue.",
    )
    parser.add_argument(
        "--min-dist",
        type=float,
        default=2.4,
        help="MINIMUM ND1/NE2-to-partner heavy-atom distance (A). Default 2.4. A "
        "partner closer than this is a steric clash, not a hydrogen bond -- a "
        "predictor can and does jam two heavy atoms to 2.07 A (MEASURED on a real "
        "design) -- so it is excluded from BOTH the qualifying and the distance-only "
        "counts, and reported separately as n_centres_too_close_target / "
        "n_centres_too_close_binder. hisbond_min_dist_target stays the "
        "UNCONDITIONAL minimum distance, so a clash is still readable there. Must be "
        ">= 0 and below --max-dist. Recorded as hisbond_min_dist.",
    )
    parser.add_argument(
        "--max-dist",
        type=float,
        default=3.5,
        help="Maximum ND1/NE2-to-partner heavy-atom distance (A) for a polar "
        "contact. Default 3.5 (the usual donor-acceptor bar). Recorded as "
        "hisbond_max_dist: a count is only comparable against another count at the "
        "same cutoff.",
    )
    parser.add_argument(
        "--max-angle",
        type=float,
        default=60.0,
        help="Maximum angle (degrees) between the imidazole N's external ring "
        "bisector -- where its H or lone pair points -- and the N-to-partner vector, "
        "for a contact to count as a qualifying H-BOND PROXY. Default 60. THE CONE IS "
        "A HEAVY-ATOM ORIENTATION PROXY, not a hydrogen-bond geometry: it says the "
        "partner sits where that nitrogen's substituent points and never looks at a "
        "hydrogen, so it cannot itself tell an H from a lone pair -- which is what "
        "keeps a count comparable between a prediction and a relaxed pose. THE "
        "DONOR/ACCEPTOR ROLE IS ANSWERED SEPARATELY by the role columns: from the "
        "TAUTOMER when the structure carries hydrogens (a relaxed PyRosetta PDB) and "
        "from the partner's own chemistry when it does not (an AF3/Boltz CIF); "
        "hisbond_role_basis says which. "
        "180 makes it a pure distance test. The distance-only counts are reported "
        "alongside (n_near_target in per_centre, and "
        "n_centres_polar_dist_only_target) so the cost of the proxy is visible. "
        "The SAME cone is applied to ND1 and to NE2 separately, and the per-nitrogen "
        "target columns (hisbond_n_centres_hb_target_nd1 / _ne2, "
        "hisbond_any_target_hb_proxy_nd1 / _ne2, hisbond_per_centre_nd1_target / "
        "_ne2_target) are always written. Recorded as hisbond_max_angle.",
    )
    parser.add_argument(
        "--env-radius",
        type=float,
        default=6.0,
        help="Radius (A) around each centre's imidazole RING CENTROID (the mean of "
        "CG, ND1, CD2, CE1, NE2) in which charged TARGET atoms are counted, as "
        "hisbond_centre<k>_n_acidic_env (ASP OD1/OD2, GLU OE1/OE2) and "
        "hisbond_centre<k>_n_basic_env (LYS NZ, ARG NE/NH1/NH2). Default 6.0, must be "
        "> 0. ATOMS ARE COUNTED, NOT RESIDUES: an Asp contributes 2 and an Arg 3. "
        "TARGET CHAINS ONLY, deliberately -- the binder's own charges are present in "
        "the free binder too, so only the target's differ between the complex and the "
        "apo binder. This is pKa CONTEXT, not a pKa: nothing here computes a shift. "
        "Recorded as hisbond_env_radius; two leaves with different radii must not be "
        "pooled.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=20,
        help="Designs measured per task (default 20). The work is milliseconds per "
        "design, so one task per design would be nearly all container start-up.",
    )


def build_hisbond_manifest(ctx: ManifestCtx[HisbondArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0

    binder = expand_chain_spec(ctx.args.binder_chain)
    if len(binder) != 1:
        raise ValueError(
            f"--binder-chain must name exactly one chain (e.g. 'D'), got "
            f"{ctx.args.binder_chain!r} -> {binder}. The centre indices are counted "
            f"along one chain; two chains have no single residue order."
        )
    target = expand_chain_spec(ctx.args.target_chains)
    if not target:
        raise ValueError("--target-chains must name at least one chain (e.g. 'A,B,C').")
    if len(set(target)) != len(target):
        raise ValueError(f"--target-chains lists a chain twice: {target}.")
    overlap = sorted(set(binder) & set(target))
    if overlap:
        raise ValueError(
            f"chain {overlap[0]} is in both --binder-chain and --target-chains; the "
            f"binder cannot be its own binding partner."
        )

    max_dist = float(ctx.args.max_dist)
    if not math.isfinite(max_dist) or max_dist <= 0:
        raise ValueError(f"--max-dist must be a positive number, got {max_dist}.")
    min_dist = float(ctx.args.min_dist)
    if not math.isfinite(min_dist) or min_dist < 0:
        raise ValueError(f"--min-dist must be a number >= 0, got {min_dist}.")
    if min_dist >= max_dist:
        raise ValueError(
            f"--min-dist ({min_dist}) must be below --max-dist ({max_dist}); "
            f"otherwise the distance window is empty and no contact can ever qualify."
        )
    max_angle = float(ctx.args.max_angle)
    if not math.isfinite(max_angle) or not 0 < max_angle <= 180:
        raise ValueError(
            f"--max-angle must be in (0, 180] degrees, got {max_angle}. Use 180 for a "
            f"pure distance test."
        )
    env_radius = float(ctx.args.env_radius)
    if not math.isfinite(env_radius) or env_radius <= 0:
        raise ValueError(
            f"--env-radius must be a positive number of angstroms, got {env_radius}. "
            f"A non-positive radius would count nothing and the environment columns "
            f"would read as an empty neighbourhood rather than as a refused question."
        )
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")

    centers_arg = ctx.args.centers.strip()
    if not centers_arg:
        raise ValueError(
            "--centers must list at least one 1-based binder position (e.g. '12,45' "
            "or '{some_position_column}')."
        )
    if any(ch in centers_arg for ch in "\t\n\r"):
        raise ValueError(
            f"--centers {centers_arg!r} contains a tab or newline; it rides the "
            f"tab-separated manifest and would split into two fields."
        )
    expect_col = (ctx.args.expect_binder_len_col or "").strip()

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"hisbond has no default input column: pass -i/--input-column with the "
            f"structure column holding the PREDICTED COMPLEX (binder + target). "
            f"Structure columns in table '{ctx.args.table}': {available or '(none)'} "
            f"-- a PARENT table's column is also valid and is resolved up the "
            f"lineage. A binder-only structure has no target and the run is "
            f"meaningless even though it would submit."
        )

    # Not ctx.ready: it drops blank-input rows, which must become error rows.
    ready = ctx.df
    status_col = status_column(ctx.out_dir.name)
    if not ctx.args.force and status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    members: list[tuple[str, str, str, str, str, str]] = []
    n_no_structure = 0
    n_no_centres = 0
    n_inherited = 0
    first_centres_error = ""
    # Never ctx.lookup for --centers: see _ColumnAwareLookup for why that is unsafe.
    centres_lookup = _ColumnAwareLookup(ctx.df, ctx.args.run_dir)
    for name in ready.index:
        name = cast(str, name)
        raw = _resolve(ctx, name, column)
        if raw is None:
            n_no_structure += 1
            print(
                f"{name}: could not resolve {column} (blank own cell, or absent from "
                f"this table and every ancestor) -- submitting anyway so the task "
                f"records an error"
            )
        # One island -> one integer; the commas stay literal. The lookup refuses to
        # inherit a blank own cell, so a design missing a centre ERRORS rather than
        # quietly borrowing an ancestor row's position.
        centres = ""
        centres_error = ""
        centres_depth = ""
        centres_lookup.reset()
        try:
            centres = ",".join(
                str(p)
                for p in parse_centres(
                    resolve_template(centers_arg, centres_lookup, name)
                )
            )
            centres_depth = str(centres_lookup.max_depth)
            if centres_lookup.max_depth > 0:
                n_inherited += 1
        except Exception as e:  # noqa: BLE001 - a per-design failure is an error row
            centres_error = " ".join(str(e).split())  # one line: it lands in a TSV cell
            first_centres_error = first_centres_error or centres_error
            n_no_centres += 1
            print(
                f"{name}: could not resolve --centers {centers_arg!r}: "
                f"{centres_error} -- submitting anyway so the task records an error"
            )
        expect_value = _resolve(ctx, name, expect_col) if expect_col else None
        if expect_col and expect_value is None:
            # Asked for but cannot be made: an error row, not a skipped check.
            print(
                f"{name}: --expect-binder-len-col {expect_col} resolves nowhere; the "
                f"binder-length check cannot be made and the design will error"
            )
            expect_value = UNRESOLVED_EXPECT

        # As-is, CIF included: gemmi's reader would merge repeated residue numbers.
        src = "" if raw is None else str(volume_path(Path(raw)))
        expect_value = expect_value or ""
        # The task row is a fixed-width TSV: a cell carrying a tab would shift every
        # later field. Refuse at submit rather than leaving a dead task behind.
        for label, field in (("the structure path", src), (f"{expect_col}", expect_value)):
            if any(ch in field for ch in "\t\n\r"):
                raise ValueError(
                    f"design {name}: {label} contains a tab or newline "
                    f"({field!r}); it rides a tab-separated task file and would "
                    f"shift every later field."
                )
        members.append((name, src, centres, centres_depth, expect_value, centres_error))

    if n_no_structure:
        print(
            f"hisbond: {n_no_structure}/{len(members)} design(s) could not resolve "
            f"{column}. Their tasks write 'error: could not resolve ...', but a row "
            f"whose OWN {column} cell is blank is skipped by collect (filter_ready) "
            f"and keeps an EMPTY hisbond_status -- compare non-empty statuses with "
            f"the row count."
        )
    if n_no_centres:
        print(
            f"hisbond: {n_no_centres}/{len(members)} design(s) could not resolve "
            f"--centers {centers_arg!r}; each will collect as 'error: ...'. A row "
            f"whose OWN table has the column but NOT a value for this design lands "
            f"here BY DESIGN: the position is never inherited from an ancestor row. "
            f"Check that every {{column}} island has a value on every row."
        )
    if n_inherited:
        print(
            f"hisbond: {n_inherited}/{len(members)} design(s) resolved at least one "
            f"--centers island from an ANCESTOR table (the column is absent from "
            f"'{ctx.args.table}' itself). That is legitimate lineage, but check "
            f"hisbond_centers_lineage_depth on the collected rows: it reports the "
            f"DEEPEST island, and a depth > 0 means the positions describe this "
            f"design only if the ancestor row is its own."
        )
    if members and n_no_structure == len(members):
        raise ValueError(
            f"hisbond resolved no structure at all for any of the {len(members)} "
            f"ready designs via column {column!r} (searched table "
            f"'{ctx.args.table}' and its ancestors). Check -i/--input-column."
        )
    if members and n_no_centres == len(members):
        raise ValueError(
            f"hisbond resolved --centers {centers_arg!r} for NONE of the "
            f"{len(members)} ready designs (searched table '{ctx.args.table}' and its "
            f"ancestors). First failure: {first_centres_error}. Without centres there "
            f"is nothing to measure, so this is a typo in --centers, not a per-design "
            f"gap."
        )

    # Delete stale results only after every guard that can refuse this submit.
    for name, _src, _centres, _depth, _expect, _err in members:
        (ctx.out_dir / f"{name}.tsv").unlink(missing_ok=True)

    tasks_dir = ctx.out_dir / "hisbond_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    per_task = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), per_task)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, src, centres, depth, expect_value, err in members[
                i : i + per_task
            ]:
                # Seven fields, always: the worker splits on tabs, it does not `cut`.
                f.write(
                    f"{name}\t{src}\t{centres}\t{depth}\t{expect_col}"
                    f"\t{expect_value}\t{err}\n"
                )
        manifest_rows.append(
            (
                # Every field non-empty, so `cut` in hisbond.sh keeps them all.
                str(volume_path(task_file)),
                binder[0],
                ",".join(target),
                str(min_dist),
                str(max_dist),
                str(max_angle),
                str(env_radius),
                column,
                # The --centers template verbatim, echoed as hisbond_centers_arg: the
                # provenance of every position (validated tab-free above).
                centers_arg,
            )
        )
    return manifest_rows
