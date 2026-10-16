#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) asking, per design: DID THE PREDICTED COMPLEX PUT
THE BINDER WHERE THE DESIGN PUT IT, allowing for the target's symmetry?

The measurement (full premise and scope limits in ``dockrmsd_worker.py``):
superpose the MODEL's TARGET chains onto the REFERENCE's TARGET chains, apply that
one transform to the model's BINDER, and report the binder RMSD with NO further
fitting -- minimised over every permutation of the target chains.

CORRESPONDENCE IS PER ROLE, and that is the whole point. The TARGET pairs by
author residue number or by gapped sequence alignment (``--target-match``); the
BINDER pairs by FILE-ORDER INDEX with equal length enforced, always, because the
two binders are the same backbone carrying DIFFERENT sequences. Applying one mode
to both roles is what disqualified DockQ here.

``action: update`` -- the binding-mode agreement is a property of a design that
already exists. CPU only: the builder forces ``gpus_per_task = 0``.

ROLES ARE FORCED BY FLAGS, NEVER INFERRED FROM LENGTH. ``--target-chains`` is the
receptor and ``--binder-chain`` the thing whose placement is measured, whatever
their relative sizes. (DockQ picks the receptor by length and therefore silently
inverts the measurement whenever a binder is longer than a target protomer.)

ROWS AND LINEAGE -- the usalign/posrmsd contract, so the three compose:

* ``--col-a`` is the MODEL (the prediction) and selects the rows: every row of
  ``-t`` whose ``--col-a`` cell is present, minus rows already ``OK`` under this
  comparison unless ``--force``. ``-i/--input-column`` is not used and is refused.
* ``--col-b`` is the REFERENCE (the design the model came from), resolved per
  design UP THE LINEAGE (the row's own value if set, else the nearest ancestor's).
* A row whose ``--col-b`` resolves to nothing is not submitted and not dropped
  silently: a ``missing_input:`` TSV is written here, so it collects as a named
  situation distinct from ``missing`` (no output) and ``error:``.

PREFIX CONTRACT (same as usalign/posrmsd). Every column is
``dockrmsd[_<dir-label>]_<prefix>_<metric>`` with ``<prefix>`` from
``--output-prefix`` or ``<col_a minus _path>_vs_<col_b minus _path>``; the
per-design files live in ``<out_dir>/<prefix>/``. ``sapia collect dockrmsd`` must
be given the same ``--col-a``/``--col-b`` or ``--output-prefix``.

Usage:
    sapia run dockrmsd outputs/RUN -t table2 \\
        --col-a boltz_path --col-b rfdiffusion3_path \\
        --binder-chain D --target-chains A,B,C \\
        --ref-binder-chain D --ref-target-chains A,B,C
"""

import csv
import math
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, cast

import pandas as pd
from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.utils import expand_chain_spec

# spec.py's default_input_column: matches no column (as for usalign/posrmsd), so
# filter_ready leaves the frame unchanged on both run and collect.
NO_INPUT_COLUMN = "not applicable"
# The prefix names a directory and is embedded in every column name.
PREFIX_RE = re.compile(r"^[A-Za-z0-9_]+$")
# k! permutations are enumerated; this caps k at 6 by default.
DEFAULT_MAX_PERMUTATIONS = 720


class DockRmsdArgs(CommonArgs):
    col_a: str
    col_b: str
    binder_chain: str
    target_chains: str
    ref_binder_chain: str
    ref_target_chains: str
    fit_atoms: str
    target_match: str
    min_target_res: int
    min_target_res_frac: float
    min_target_identity: float
    fit_tolerance: float
    max_target_fit_rmsd: float
    max_permutations: int
    pin_mapping: bool
    designs_per_task: int
    output_prefix: str | None


def add_run_dockrmsd_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--col-a",
        type=str,
        required=True,
        help="Column holding the MODEL complex (the prediction) -- the structure "
        "that is MOVED onto the reference. Selects which rows run; read from the "
        "-t table itself, not up the lineage.",
    )
    parser.add_argument(
        "--col-b",
        type=str,
        required=True,
        help="Column holding the REFERENCE complex (the design the model came "
        "from) -- the fixed frame. Resolved per design up the lineage, like "
        "usalign's --col-b.",
    )
    parser.add_argument(
        "--binder-chain",
        type=str,
        required=True,
        help="The binder chain in the MODEL, exactly one chain id (gemmi's chain "
        "name = auth_asym_id). NO DEFAULT: a defaulted letter would measure the "
        "placement of the wrong molecule.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        required=True,
        help="The target chains in the MODEL, chain mini-language ('A,B,C' or "
        "'A:C'). NO DEFAULT. These and only these atoms determine the transform.",
    )
    parser.add_argument(
        "--ref-binder-chain",
        type=str,
        required=True,
        help="The binder chain in the REFERENCE, exactly one chain id. NO DEFAULT "
        "-- generators renumber and relabel chains, so it often differs from "
        "--binder-chain.",
    )
    parser.add_argument(
        "--ref-target-chains",
        type=str,
        required=True,
        help="The target chains in the REFERENCE; same COUNT as --target-chains. "
        "NO DEFAULT. Which model chain maps onto which reference chain is SEARCHED "
        "over all permutations, not taken from the order given here.",
    )
    parser.add_argument(
        "--fit-atoms",
        choices=("ca", "backbone"),
        default="ca",
        help="Atoms used for BOTH the target fit and the binder RMSD: 'ca' "
        "(DEFAULT, one atom per residue) or 'backbone' (N, CA, C, O). Numbers from "
        "the two settings are not comparable.",
    )
    parser.add_argument(
        "--target-match",
        choices=("number", "sequence"),
        default="number",
        help="How TARGET residues correspond. 'number' (DEFAULT): intersect "
        "author residue numbers per chain pair -- right when the two files share "
        "a numbering convention. 'sequence': gapped global alignment -- right "
        "when one file MODELS RESIDUES THE OTHER LACKS (an unresolved loop), "
        "which NO renumbering can fix because both files may already be numbered "
        "1..N. **Pass 'sequence' whenever model and reference can differ in "
        "modelled extent**: on this campaign's pilot pair (model target 1-149, "
        "reference 1-139, a 10-residue loop modelled only in the prediction) "
        "'number' gives 0.698 identity and errors every row. The BINDER is always "
        "paired by file-order index and is NEVER sequence-aligned, whatever this "
        "flag says. No silent fallback: the mode is recorded as target_match.",
    )
    parser.add_argument(
        "--min-target-res",
        type=int,
        default=10,
        help="A chain permutation in which ANY SINGLE chain pair contributes "
        "fewer than this many residue pairs is EXCLUDED (default 10). Applied PER "
        "CHAIN PAIR, never to the pooled sum: pooling would admit a mapping in "
        "which one pair carries the whole correspondence and the others nothing, "
        "i.e. a single-protomer fit reported as a whole-target one. A row where "
        "no permutation clears it is an error.",
    )
    parser.add_argument(
        "--min-target-res-frac",
        type=float,
        default=0.75,
        help="If the winning permutation used a smaller fraction of the shorter "
        "target than this, the row is 'warn:' (default 0.75) -- numbers kept, "
        "excluded by == 'OK'. It means the two files correspond over part of the "
        "target only; compare target_resnum_range_model with "
        "target_resnum_range_ref.",
    )
    parser.add_argument(
        "--min-target-identity",
        type=float,
        default=0.8,
        help="A chain permutation whose shared target residues agree in identity "
        "at a lower fraction than this is EXCLUDED (default 0.8). On a homo-"
        "oligomer every permutation passes; on a hetero target only the right "
        "pairings do. A row where none clears it is an error.",
    )
    parser.add_argument(
        "--fit-tolerance",
        type=float,
        default=1.0,
        help="A permutation competes for the minimum only if its target fit RMSD "
        "is within this many angstroms of the BEST permutation's (default 1.0). "
        "This is what keeps the minimum over genuine symmetry copies: on a C3 "
        "trimer three permutations fit to ~1 A and three to ~22 A, and the latter "
        "must not be allowed to win on an accidental binder RMSD.",
    )
    parser.add_argument(
        "--max-target-fit-rmsd",
        type=float,
        default=5.0,
        help="Above this target fit RMSD the row is 'warn:' -- numbers kept, "
        "excluded by == 'OK' (default 5.0). The model did not reproduce the "
        "target, so the frame the binder was moved into is not the reference "
        "frame and the docking RMSD means nothing at any magnitude.",
    )
    parser.add_argument(
        "--max-permutations",
        type=int,
        default=DEFAULT_MAX_PERMUTATIONS,
        help=f"Refuse the submit if the target chain count k gives more than this "
        f"many permutations k! (default {DEFAULT_MAX_PERMUTATIONS}, i.e. k <= 6).",
    )
    parser.add_argument(
        "--pin-mapping",
        action="store_true",
        help="Evaluate ONLY the chain order as given, with no permutation search "
        "(default off). Use on a PSEUDO-symmetric target -- protomers homologous "
        "enough to clear --min-target-identity but NOT physically interchangeable "
        "-- where a minimum over mappings would report the best of several "
        "INEQUIVALENT sites, scoring a binder docked at the wrong protomer as "
        "though it hit the designed one. Not needed on a true homo-oligomer.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=25,
        help="Designs per task (default 25; the measurement is well under a second "
        "each, so the container start dominates).",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="Prefix for the results subdir and every column. Defaults to "
        "'<col_a minus _path>_vs_<col_b minus _path>'. collect must be given the "
        "same value (or the same --col-a/--col-b).",
    )


def resolve_prefix(col_a: str | None, col_b: str | None, output_prefix: str | None) -> str:
    """The comparison prefix, shared verbatim with the collector."""
    if output_prefix:
        prefix = output_prefix.strip()
    elif col_a and col_b:
        prefix = f"{col_a.strip().removesuffix('_path')}_vs_{col_b.strip().removesuffix('_path')}"
    else:
        raise ValueError("Provide --output-prefix, or both --col-a and --col-b.")
    if not PREFIX_RE.match(prefix):
        raise ValueError(
            f"comparison prefix {prefix!r} must be letters, digits and '_' only "
            f"(it names a directory and every column); pass --output-prefix"
        )
    return prefix


def _chains(flag: str, spec: str) -> list[str]:
    chains = expand_chain_spec(spec)
    if not chains:
        raise ValueError(f"{flag} must name at least one chain, got {spec!r}.")
    if len(set(chains)) != len(chains):
        raise ValueError(f"{flag} lists a chain twice: {chains}.")
    return chains


def _lookup_str(lookup: Any, name: str, column: str) -> str | None:
    """A lineage-resolved cell as a stripped string, or None when absent/blank."""
    value = lookup(name, column)
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return None
    text = str(value).strip()
    return text or None


def _write_status_tsv(path: Path, name: str, status: str) -> None:
    """A per-design result written at submit time (no task ran). The collector
    NAs every metric on it."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["name", "status", "sup_path"])
        writer.writerow([name, status.replace("\t", " ").replace("\n", " "), ""])


def _validate(args: DockRmsdArgs) -> tuple[list[str], list[str], str, str]:
    """Every submit-time refusal that does not need the table."""
    binder = _chains("--binder-chain", args.binder_chain)
    ref_binder = _chains("--ref-binder-chain", args.ref_binder_chain)
    for flag, got in (("--binder-chain", binder), ("--ref-binder-chain", ref_binder)):
        if len(got) != 1:
            raise ValueError(f"{flag} must name exactly one chain, got {got}.")
    targets = _chains("--target-chains", args.target_chains)
    ref_targets = _chains("--ref-target-chains", args.ref_target_chains)
    if binder[0] in targets:
        raise ValueError(
            f"--binder-chain {binder[0]} is also in --target-chains {targets}: the "
            f"binder would determine the transform it is then measured against."
        )
    if ref_binder[0] in ref_targets:
        raise ValueError(
            f"--ref-binder-chain {ref_binder[0]} is also in --ref-target-chains "
            f"{ref_targets}."
        )
    if len(targets) != len(ref_targets):
        raise ValueError(
            f"--target-chains has {len(targets)} chains and --ref-target-chains "
            f"{len(ref_targets)}: a chain permutation is a bijection, so the two "
            f"must name the same number of target chains."
        )
    n_perms = math.factorial(len(targets))
    if not args.pin_mapping and n_perms > args.max_permutations:
        raise ValueError(
            f"{len(targets)} target chains give {n_perms} permutations, over "
            f"--max-permutations ({args.max_permutations}). Raise the flag only if "
            f"you mean to enumerate that many fits."
        )
    if args.min_target_res < 3:
        raise ValueError("--min-target-res must be >= 3 (a Kabsch fit needs 3 points).")
    res_frac = float(args.min_target_res_frac)
    if not math.isfinite(res_frac) or not 0.0 <= res_frac <= 1.0:
        raise ValueError(f"--min-target-res-frac must be in [0, 1], got {res_frac}.")
    identity = float(args.min_target_identity)
    if not math.isfinite(identity) or not 0.0 <= identity <= 1.0:
        raise ValueError(f"--min-target-identity must be in [0, 1], got {identity}.")
    if not math.isfinite(args.fit_tolerance) or args.fit_tolerance < 0:
        raise ValueError("--fit-tolerance must be >= 0.")
    if not math.isfinite(args.max_target_fit_rmsd) or args.max_target_fit_rmsd <= 0:
        raise ValueError("--max-target-fit-rmsd must be > 0.")
    if args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")
    return targets, ref_targets, binder[0], ref_binder[0]


def build_dockrmsd_manifest(ctx: ManifestCtx[DockRmsdArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0

    args = ctx.args
    col_a = (args.col_a or "").strip()
    col_b = (args.col_b or "").strip()
    if not col_a or not col_b:
        raise ValueError("--col-a (model) and --col-b (reference) are both required.")
    if col_a == col_b:
        raise ValueError(
            f"--col-a and --col-b are both {col_a!r}: every design would be "
            f"compared with itself (docking RMSD 0 by construction)."
        )
    if (getattr(args, "input_column", NO_INPUT_COLUMN) or NO_INPUT_COLUMN) != NO_INPUT_COLUMN:
        raise ValueError(
            f"dockrmsd selects rows by --col-a; -i/--input-column "
            f"({args.input_column!r}) is not used and would make collect iterate a "
            f"different row set. Drop -i."
        )
    targets, ref_targets, binder, ref_binder = _validate(args)
    prefix = resolve_prefix(col_a, col_b, args.output_prefix)

    df = ctx.df
    if col_a not in df.columns:
        available = ", ".join(sorted(map(str, df.columns))) or "<none>"
        raise ValueError(
            f"--col-a {col_a!r} is not a column of the -t table (it selects the "
            f"rows, so it is not resolved up the lineage). Available: {available}"
        )

    cells = df[col_a]
    ready = df[cells.notna() & (cells.astype(str).str.strip() != "")]
    # The driver leaf-prefixes collected data, so this comparison's status lives
    # at <leaf>_<prefix>_status (as for usalign/posrmsd).
    status_col = f"{ctx.out_dir.name}_{prefix}_status"
    if args.force:
        print("Re-running dockrmsd for all selected designs (including OK).")
    elif status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    results_dir = ctx.out_dir / prefix

    plan: list[tuple[str, str, str | None]] = []
    for name in ready.index:
        name = cast(str, name)
        src_a = str(ready.at[name, col_a]).strip()
        plan.append((name, src_a, _lookup_str(ctx.lookup, name, col_b)))

    if plan and all(b is None for _, _, b in plan):
        raise ValueError(
            f"--col-b {col_b!r} has no value on any of the {len(plan)} selected "
            f"row(s) or their ancestors -- a typo, or the column lives on a table "
            f"that is not an ancestor of -t. Nothing submitted."
        )

    members: list[tuple[str, str, str]] = []
    n_missing = 0
    refs: set[str] = set()
    for name, src_a, src_b in plan:
        result_tsv = results_dir / f"{name}.tsv"
        if src_b is None:
            n_missing += 1
            status = f"missing_input: --col-b {col_b!r} has no value on this row or any ancestor"
            print(f"{name}: {status}")
            _write_status_tsv(result_tsv, name, status)
            continue
        # A stale result must not survive a task killed before it writes.
        result_tsv.unlink(missing_ok=True)
        (results_dir / f"{name}_sup.cif").unlink(missing_ok=True)
        refs.add(src_b)
        members.append((name, str(volume_path(Path(src_a))), str(volume_path(Path(src_b)))))

    if n_missing:
        print(
            f"WARNING: {n_missing} design(s) not submitted for a missing --col-b; "
            f"each has a missing_input row in {results_dir} and will collect as such."
        )
    print(
        f"dockrmsd: {len(members)} design(s) against {len(refs)} distinct reference "
        f"structure(s) from {col_b!r}; model binder {binder} vs target "
        f"{','.join(targets)}, reference binder {ref_binder} vs target "
        f"{','.join(ref_targets)}; "
        f"{1 if args.pin_mapping else math.factorial(len(targets))} mapping(s) per "
        f"design, fit atoms {args.fit_atoms}, target match {args.target_match}."
    )

    tasks_dir = ctx.out_dir / f"{prefix}_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    per_task = args.designs_per_task
    rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), per_task)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, model, ref in members[i : i + per_task]:
                f.write(f"{name}\t{model}\t{ref}\n")
        rows.append(
            (
                # Every field is non-empty, so none can vanish from `cut`.
                str(volume_path(task_file)),
                binder,
                ",".join(targets),
                ref_binder,
                ",".join(ref_targets),
                args.fit_atoms,
                args.target_match,
                str(int(args.min_target_res)),
                str(float(args.min_target_res_frac)),
                str(float(args.min_target_identity)),
                str(float(args.fit_tolerance)),
                str(float(args.max_target_fit_rmsd)),
                "1" if args.pin_mapping else "0",
                prefix,
            )
        )
    return rows
