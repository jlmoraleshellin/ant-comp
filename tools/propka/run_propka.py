#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) asking, per design: what pKa does PROPKA3 predict for
each designed histidine of the binder in the complex and in the free binder, and how big
is the shift? Flags, traps and how to read the columns: ``.claude/skills/propka/SKILL.md``.

``action: update``: the table is annotated in place and no child is minted. Run it twice
with different ``-l/--dir-label`` to score two structure columns onto one table. CPU only
(the builder forces ``gpus_per_task = 0``).

THE FREE STATE IS A RIGID UNBINDING: the binder chains' records cut out of the same file,
same coordinates, no relaxation (``free_state_model`` says so in every row). A file with
no chain outside the binder chains is an error row, unless ``--free-only`` asks for the
free state alone.

CENTRES are named by exactly one of ``--resnums`` (a spec, ``{expr}`` islands allowed, so
table columns work as ``'{col1},{col2}'``) or ``--auto-centres`` (every histidine of the
binder chains). Each must be a histidine, uniquely numbered in its chain.

``default_input_column`` is the ``"not applicable"`` sentinel (as in epitope, cms,
chainsel, hisbond), so the builder raises unless ``-i/--input-column`` is given.

CIF IN, PDB OUT, HERE: PROPKA 3.5.1 reads ``.pdb`` only and the worker's image has no
prosapia, so a CIF is converted here with ``ensure_pdb`` (gemmi, auth chain letters and
numbers kept, cached under ``run_dir/.cif_to_pdb``). The collected ``<leaf>_path`` is the
file PROPKA read.

The structure, the centres and the expected length resolve as in hisbond: if this table
HAS the column the row's own cell is the answer and a blank never falls back to an
ancestor (that would measure another design); only an ABSENT column walks the lineage. A
design that resolves nothing is SUBMITTED ANYWAY and lands as ``error: could not resolve
...``, never dropped: silently submitting fewer designs than the table holds is this
workspace's most expensive quiet failure.
"""

import shutil
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

# Every shared constant and parser comes from the worker, so a submit can never accept a
# list the worker would refuse. The import runs on this side only: the worker imports
# nothing but the standard library at module level (it also runs in the tool image,
# which has no prosapia) and must stay that way, because a module-level failure there
# breaks tool DISCOVERY for every sapia command.
from .propka_worker import (  # noqa: E402
    AUTO_CENTRES,
    DEFAULT_PH_VALUES as _DEFAULT_PH_TUPLE,
    DEFAULT_REFERENCE_PH as _DEFAULT_REFERENCE_PH,
    DEFAULT_RELEASE_FACTORS as _DEFAULT_FACTOR_TUPLE,
    MAX_CENTRE_COLUMNS,
    MAX_PH_VALUES,
    MAX_RELEASE_FACTORS,
    UNRESOLVED_EXPECT,
    num_token,
    parse_centres,
    parse_ph_values,
    parse_reference_ph,
    parse_release_factors,
)

DEFAULT_PH_VALUES = ",".join(num_token(v, True) for v in _DEFAULT_PH_TUPLE)
DEFAULT_RELEASE_FACTORS = ",".join(num_token(v, False) for v in _DEFAULT_FACTOR_TUPLE)
DEFAULT_REFERENCE_PH = num_token(_DEFAULT_REFERENCE_PH, True)


class PropkaArgs(CommonArgs):
    binder_chains: str
    resnums: str | None
    auto_centres: bool
    free_only: bool
    expect_binder_len_col: str | None
    designs_per_task: int
    keep_propka_input: bool
    ph_values: str
    reference_ph: str
    release_factors: str


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


def _resolve(ctx: ManifestCtx[PropkaArgs], name: str, column: str) -> str | None:
    """The row's value for ``column``, or None when it cannot be resolved.

    If this table HAS ``column`` the row's own cell is the answer; a blank never falls
    back to an ancestor (that would describe another design). Only an ABSENT column
    walks the lineage, nearest ancestor first.
    """
    if column in ctx.df.columns:
        if name not in ctx.df.index:
            return None
        return _clean(ctx.df.at[name, column])
    return _clean(ctx.lookup(name, column))


def add_run_propka_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--binder-chains",
        type=str,
        required=True,
        help="The BINDER chain(s), chain mini-language ('D', 'A,B', 'A:C'). NO "
        "DEFAULT: which letter the binder is depends on how the complex was built, and "
        "a defaulted letter would compute the free-state pKa of the WRONG MOLECULE with "
        "status OK. These chains together, alone, are the free state (a homo-oligomeric "
        "binder is the whole oligomer), and the centres are looked up in them. Check "
        "propka_binder_chains_n_res, or pass --expect-binder-len-col.",
    )
    parser.add_argument(
        "--resnums",
        type=str,
        default=None,
        help="CENTRES AS A SPEC: '45,52' (one binder chain) or 'A:45,B:45' (any), in "
        "auth_seq_id numbering, i.e. RESIDUE NUMBERS as a viewer shows them, NOT 1..L "
        "indices. '{expr}' islands are resolved per design up the lineage (integers, "
        "column names, + - * //; the nearest non-null value wins), so table columns "
        "are '{col1},{col2}' and a multi-chain binder 'A:{colA},B:{colB}'. A design "
        "whose island resolves nothing (an NA, a missing column) is an error row: a "
        "varying number of centres per design cannot be expressed, use --auto-centres. "
        "Exactly one of --resnums, --auto-centres is required.",
    )
    parser.add_argument(
        "--auto-centres",
        action="store_true",
        help="CENTRES = EVERY HISTIDINE of the binder chains, in file order. The "
        f"first {MAX_CENTRE_COLUMNS} get per-centre columns (centres_reported says how "
        "many did); every other column, binder_his_pka included, covers all of them.",
    )
    parser.add_argument(
        "--free-only",
        action="store_true",
        help="Run the FREE STATE ALONE: one PROPKA run on the binder chains, no partner "
        "chain needed (the file may hold only the binder). Everything that needs the "
        "complex (pka_complex, dpka, switch ratio, release pH, n_centres_coupled, ...) "
        "is NA by construction and free_only is 1. Without it, a file with no chain "
        "outside the binder chains is an error row.",
    )
    parser.add_argument(
        "--expect-binder-len-col",
        type=str,
        default=None,
        help="Optional table column holding the EXPECTED total residue count of the "
        "binder chains: an integer (e.g. 'atomiumph_seq_len') or a one-letter sequence "
        "('/' chain breaks allowed), resolved up the lineage. When given, a different "
        "count is an ERROR for that design: that is the check that catches a chain "
        "letter which EXISTS but is the wrong molecule. A check asked for and "
        "unresolvable is also an error row. Pass it whenever the column exists.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=10,
        help="Designs per task (default 10). Each design is TWO PROPKA runs (one under "
        "--free-only) and the complex run is the slow one: its coupling scan is O(n^2) "
        "in the titratable groups of the whole assembly. One measurement: both runs on "
        "a 6-chain, ~590-residue complex took about 2 s, so 10 designs is a short task.",
    )
    parser.add_argument(
        "--ph-values",
        type=str,
        default=DEFAULT_PH_VALUES,
        help=f"Comma-separated pH values at which to report the WEAKENING FACTOR F, the "
        f"factor by which binding is weakened relative to the high-pH limit, "
        f"F = prod_i (1 + 10^(pKa_free_i - pH)) / (1 + 10^(pKa_complex_i - pH)) over "
        f"the histidines. Default '{DEFAULT_PH_VALUES}'. Each between 0 and 14, at most "
        f"{MAX_PH_VALUES}, no duplicates. EACH VALUE ADDS TEN COLUMNS (five at the "
        f"reference pH), named with the pH: 6.5 -> f_ph6p5_{{designed,all}}, "
        f"f_ceiling_ph6p5_{{designed,all}}, switch_ph6p5_{{designed,all}}, "
        f"stab_dg_fold_ph6p5, stab_loss_ph6p5 and its _designed/_other shares. No "
        f"extra PROPKA run: arithmetic on the pKa values of the two runs.",
    )
    parser.add_argument(
        "--reference-ph",
        type=str,
        default=DEFAULT_REFERENCE_PH,
        help=f"The DENOMINATOR of the switch ratio: the pH the binder is meant to HOLD "
        f"at. Default {DEFAULT_REFERENCE_PH}. MUST be one of --ph-values. For every "
        f"OTHER pH the tool writes switch_ph<tag>_designed and switch_ph<tag>_all = "
        f"F(pH) / F(this pH), which under the same model is K_D(pH) / K_D(reference): "
        f"THE number a pH switch is judged on. The reference pH gets no ratio column.",
    )
    parser.add_argument(
        "--release-factors",
        type=str,
        default=DEFAULT_RELEASE_FACTORS,
        help=f"Comma-separated weakening factors X for the release-pH columns: "
        f"ph_release_x<X>_designed / _all = THE HIGHEST pH IN [3, 9] AT WHICH F REACHES "
        f"X. F is NOT monotonic in pH, so this does NOT say F >= X everywhere below it: "
        f"f_dips_below_x<X>_<set> is 1 where F falls back below X further down. GATE ON "
        f"f_ph<tag>_<set> AT THE pH YOU CARE ABOUT; use the release pH to rank. Default "
        f"'{DEFAULT_RELEASE_FACTORS}'. Each > 1, at most {MAX_RELEASE_FACTORS}, no "
        f"duplicates. Each adds FOUR columns (two release, two dips flags).",
    )
    parser.add_argument(
        "--keep-propka-input",
        action="store_true",
        help="Keep the temporary PDB files PROPKA read, as <out_dir>/<name>_propka/"
        "complex.pdb and free.pdb (what `propka3 <file>` must reproduce by hand; "
        "hydrogens already stripped). Off by default. THE INPUT FILE IS NEVER MODIFIED.",
    )


def build_propka_manifest(ctx: ManifestCtx[PropkaArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0

    binder_chains = expand_chain_spec(ctx.args.binder_chains)
    if not binder_chains:
        raise ValueError("--binder-chains must name at least one chain (e.g. 'D').")
    if len(set(binder_chains)) != len(binder_chains):
        raise ValueError(f"--binder-chains lists a chain twice: {binder_chains}.")

    literal = (ctx.args.resnums or "").strip()
    if bool(literal) + bool(ctx.args.auto_centres) != 1:
        raise ValueError(
            "name the centres with exactly one of --resnums (a spec, '{expr}' islands "
            "allowed) or --auto-centres (every histidine of the binder chains)."
        )
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")
    if literal:
        if any(ch in literal for ch in "\t\r\n"):
            raise ValueError("--resnums must not contain a tab or a newline.")
        centres_source = f"resnums:{literal}"
        if "{" not in literal:
            parse_centres(literal, binder_chains, centres_source)  # syntax, up front
    else:
        centres_source = AUTO_CENTRES

    # Refused HERE, before a single task is queued: a bad pH list would otherwise cost a
    # whole array of containers. The worker re-validates the same way.
    ph_values = parse_ph_values(ctx.args.ph_values)
    factors = parse_release_factors(ctx.args.release_factors)
    reference_ph = parse_reference_ph(ctx.args.reference_ph, ph_values)
    ph_values_token = ",".join(num_token(v, True) for v in ph_values)
    factors_token = ",".join(num_token(v, False) for v in factors)
    reference_token = num_token(reference_ph, True)

    expect_col = (ctx.args.expect_binder_len_col or "").strip()

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"propka has no default input column: pass -i/--input-column with the "
            f"structure column holding the COMPLEX (binder + target) whose pKa you "
            f"want. Structure columns in table '{ctx.args.table}': "
            f"{available or '(none)'} -- a PARENT table's column is also valid and is "
            f"resolved up the lineage. (A binder-only file needs --free-only.)"
        )

    # Not ctx.ready: it drops blank-input rows, which must become error rows.
    ready = ctx.df
    status_col = status_column(ctx.out_dir.name)
    if not ctx.args.force and status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    members: list[tuple[str, str, str, str]] = []
    n_no_structure = 0
    n_no_centres = 0
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

        if literal:
            cell = literal
            if "{" in literal:
                try:
                    cell = resolve_template(literal, ctx.lookup, name)
                except ValueError as exc:
                    print(f"{name}: {exc} -- submitting anyway so the task records an error")
                    cell = ""
        else:
            cell = AUTO_CENTRES
        if not cell.strip():
            n_no_centres += 1
            cell = ""
            print(
                f"{name}: no centre residue number resolved from {centres_source} -- "
                f"submitting anyway so the task records an error"
            )

        expect_value = _resolve(ctx, name, expect_col) if expect_col else None
        if expect_col and expect_value is None:
            # Asked for but cannot be made: an error row, not a skipped check.
            print(
                f"{name}: --expect-binder-len-col {expect_col} resolves nowhere; the "
                f"binder-length check cannot be made and the design will error"
            )
            expect_value = UNRESOLVED_EXPECT

        src = ""
        if raw is not None:
            path = Path(raw)
            # PROPKA reads PDB only; convert here, where prosapia (and gemmi) exist. A
            # path that does not exist is passed raw: the worker records the error.
            if path.exists():
                path = ensure_pdb(path, ctx.args.run_dir)
            src = str(volume_path(path))
        members.append((name, src, cell, expect_value or ""))

    if n_no_structure:
        print(
            f"propka: {n_no_structure}/{len(members)} design(s) could not resolve "
            f"{column}. Their tasks write 'error: could not resolve ...', but a row "
            f"whose OWN {column} cell is blank is skipped by collect (filter_ready) and "
            f"keeps an EMPTY propka_status -- compare non-empty statuses with the row "
            f"count."
        )
    if n_no_centres:
        print(
            f"propka: {n_no_centres}/{len(members)} design(s) resolved no centre "
            f"residue number from {centres_source}; each will collect as 'error: no "
            f"centre residue number ...'."
        )
    if members and n_no_structure == len(members):
        raise ValueError(
            f"propka resolved no structure at all for any of the {len(members)} ready "
            f"designs via column {column!r} (searched table '{ctx.args.table}' and its "
            f"ancestors). Check -i/--input-column."
        )
    if members and n_no_centres == len(members):
        raise ValueError(
            f"propka resolved no centre residue number at all for any of the "
            f"{len(members)} ready designs from {centres_source} (searched table "
            f"'{ctx.args.table}' and its ancestors). Check the centre flags."
        )

    # Only NOW delete the stale per-design results, after every guard that can refuse
    # this submit; a task killed after this point collects as 'missing', never as the
    # previous run's OK. The kept PROPKA inputs go too: a stale pair would be reproduced
    # against the NEW row's numbers.
    for name, _src, _cell, _expect in members:
        (ctx.out_dir / f"{name}.tsv").unlink(missing_ok=True)
        shutil.rmtree(ctx.out_dir / f"{name}_propka", ignore_errors=True)

    tasks_dir = ctx.out_dir / "propka_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    per_task = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), per_task)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, src, cell, expect_value in members[i : i + per_task]:
                # expect_col/expect_value ride last: the only fields that may be empty.
                f.write(f"{name}\t{src}\t{cell}\t{expect_col}\t{expect_value}\n")
        manifest_rows.append(
            (
                # Every top-level field is non-empty, so `cut` never loses one.
                str(volume_path(task_file)),
                ",".join(binder_chains),
                column,
                centres_source,
                "1" if ctx.args.keep_propka_input else "0",
                ph_values_token,
                factors_token,
                reference_token,
                "1" if ctx.args.free_only else "0",
            )
        )
    return manifest_rows
