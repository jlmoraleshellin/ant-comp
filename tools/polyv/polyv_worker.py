#!/usr/bin/env python3
"""
Poly-valine interface geometry for a batch of binder/target complexes (PyRosetta only).

This is the per-array-task step of the polyv tool. One task triages several designs,
because ``pyrosetta.init()`` costs ~3-4 s per container while a design costs seconds:
packing designs into a task pays the init once. Each design gets its own try/except,
so one bad structure never costs the rest of the batch.

PREMISE -- sequence-blind backbone triage: *does this backbone present a large,
well-fitted surface to a fixed target, before any sequence exists?* Per design:

  0. reseed Rosetta's RNG from a hash of the DESIGN NAME, so the stochastic packer
     gives one input file one answer regardless of its position in the batch, the
     batch size or the row order (the seed is collected as ``seed``);
  1. load the structure (.pdb or .cif) with PyRosetta;
  2. validate that every requested binder and target chain is present;
  3. drop any chain that is in neither list, and any NON-POLYMER residue named in
     --exclude-resnames (counted as ``n_res_dropped_chain`` /
     ``n_res_dropped_resname``);
  4. mutate every binder residue to VAL (GLY and PRO included unless --keep-gly-pro);
  5. repack ONLY the binder's valine side chains -- the target is completely frozen,
     no repack, no minimisation, no relax, anywhere;
  6. InterfaceAnalyzerMover (binder vs target) with pack_separated OFF, dSASA ON,
     interface SC ON, and interface energy / hbond-SASA energy / packstat OFF;
  7. Rosetta's ContactMolecularSurface filter in BOTH directions.

It is NOT an interface-quality score for a designed sequence. Its numbers are not
comparable with the ``cms`` tool's (cms-cuda, as-built sequence) nor with
``pyrosetta``'s ``if_dSASA`` (relaxed pose, separated repack).

No energy, dG or hydrogen-bond count is computed or reported. ref2015 is used only to
place valine rotamers during the repack.

``dsasa`` is InterfaceAnalyzerMover's delta SASA summed over BOTH partners, i.e. about
TWICE the single-side buried surface area usually quoted in the binder literature.

Error discipline (as in cms/chainsel): **a requested chain that is absent is an error,
not a smaller selection**; a requested --exclude-resnames that matches nothing is an
error for the same reason; a polymer residue may not be excluded at all (it would
leave a hole and change every remaining residue's surface); and an incomplete polyV
mutation is an error, not a quietly mixed-sequence pose -- ``n_mutated + n_kept_gly_pro == n_res_binder`` is enforced on
both the default and the --keep-gly-pro path. A dSASA or CMS of 0 is NOT an error --
it is an honest, alarming measurement. Anything that could not be computed is written
as NA (an empty cell), never as 0.

Writes, per design, ``<name>.pdb`` (the packed polyV pose) and ``<name>.tsv``
(one row: name, status, path, metrics...).

Normally invoked per array task by the polyv tool -- run ``sapia run polyv`` rather
than calling this directly.

Usage (standalone):
    printf 'design_0\\tcomplex.pdb\\n' > task.tsv
    python tools/polyv/polyv_worker.py --task-file task.tsv \\
        --binder-chains B --target-chains A --out-dir out/
"""

import argparse
import csv
import hashlib
import math
import sys
import time
from pathlib import Path
from typing import Any

# PyRosetta init. Rationale for the non-obvious flags:
#   -detect_disulf false  a binder backbone becomes poly-valine, so its cysteines are
#                         about to disappear; auto-detected disulfides would make the
#                         mutation fail. Nothing measured here is energetic, so this
#                         changes no reported number.
#   -load_PDB_components false  never reach for the chemical component dictionary.
INIT_FLAGS = (
    "-mute all -ignore_unrecognized_res -ignore_zero_occupancy false "
    "-detect_disulf false -load_PDB_components false"
)

# Metric columns of the per-design TSV, in order (bare names: collect_polyv.py hands
# them to the driver, which leaf-prefixes them to polyv_<name>).
METRIC_COLUMNS = [
    # measurements
    "cms_target",
    "cms_binder",
    "dsasa",
    "dsasa_per_res",
    "nres_int",
    "nres_int_binder",
    "nres_int_target",
    "sc",
    # trust metrics
    "n_mutated",
    "n_kept_gly_pro",
    "n_res_binder",
    "n_res_target",
    "n_res_dropped_chain",
    "n_res_dropped_resname",
    "chains_found",
    "n_input_restypes",
    "seed",
    "seconds",
]
RESULT_COLUMNS = ["name", "status", "path", *METRIC_COLUMNS]

# A chain ID that is blank in the file is rendered as this in `chains_found`. It is
# REFUSED in --binder-chains/--target-chains (both here and at submit time): it is also
# the separator in InterfaceAnalyzerMover's partner string, so a blank-chain selection
# could not be expressed there even if we wanted it.
BLANK_CHAIN = "_"


def split_list(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def _call_first(obj: Any, names: list[str], *args: Any, what: str) -> None:
    """Call the first method of `obj` named in `names`, or raise with the real API.

    PyRosetta's binding names are the one thing that cannot be checked from here, so a
    mismatch must fail loudly and name what is actually available -- never fall back to
    a default, which would silently change what the number means.
    """
    for name in names:
        fn = getattr(obj, name, None)
        if fn is not None:
            fn(*args)
            return
    available = sorted(m for m in dir(obj) if not m.startswith("_"))
    raise RuntimeError(
        f"PyRosetta API mismatch: {type(obj).__name__} has none of {names} ({what}). "
        f"Available methods: {', '.join(available)}"
    )


def chain_ids(pose: Any) -> list[str]:
    """Per-residue PDB chain ID, 1-indexed like the pose (index 0 is a placeholder)."""
    info = pose.pdb_info()
    if info is None:
        raise ValueError("pose has no PDB info (chain IDs unavailable)")
    out = [""]
    for i in range(1, pose.total_residue() + 1):
        c = info.chain(i)
        out.append(BLANK_CHAIN if c in ("", " ") else str(c))
    return out


def ordered_unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(values))


def check_contiguous(per_res_chain: list[str]) -> None:
    """Every chain must occupy one contiguous block of residues.

    A chain that reappears later means the file holds several models (or several
    copies under one letter); every count below would then be silently doubled.
    """
    seen: set[str] = set()
    previous = None
    for c in per_res_chain[1:]:
        if c != previous:
            if c in seen:
                raise ValueError(
                    f"chain {c} is not contiguous (several models, or a repeated "
                    f"chain ID, in one file)"
                )
            seen.add(c)
            previous = c


def design_seed(name: str) -> int:
    """A stable 31-bit RNG seed derived only from the design name.

    Rosetta's packer is stochastic and its RNG state carries ACROSS designs inside one
    task, so without this a design's score depends on its position in the batch and on
    --designs-per-task (measured: the same file scored 402.4235 and 402.1243 at two
    positions of one task). Seeding per design from the name makes one input file give
    one answer, whatever the batch size or row order. ``seed`` is a collected column so
    a number can be traced to the state that produced it.

    SHA-256 rather than ``hash()``: Python's string hash is salted per process.
    """
    digest = hashlib.sha256(name.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "big") & 0x7FFFFFFF


def seed_rosetta_rng(seed: int) -> None:
    """Reset Rosetta's global RNG. Called once per design, before anything stochastic."""
    from pyrosetta.rosetta.numeric import random as rosetta_random

    rosetta_random.rg().set_seed(int(seed))


def drop_unscored_residues(
    pose: Any, keep_chains: set[str], exclude_resnames: set[str]
) -> tuple[int, int]:
    """Delete, in place, every residue that is not part of the two scored sides.

    That is: every residue whose chain is in neither selection, plus every residue on a
    scored chain whose name is in --exclude-resnames. Returns
    ``(n_dropped_chain, n_dropped_resname)`` -- collected separately so that "a whole
    chain went" and "the exclusion flag bit" are distinguishable, and so that a zero in
    the second against a non-empty flag is visible.

    InterfaceAnalyzerMover partitions the pose into exactly two partners, so a third
    chain would be in neither; and selection is otherwise purely by chain letter, so a
    heteroatom residue Rosetta accepted on a named chain would contribute to dSASA,
    CMS and the residue counts unless it is named here.

    Two errors, both deliberate, both mirroring the absent-chain contract:

    * **a POLYMER residue may not be excluded.** Removing an interior residue of a
      polymer chain exposes backbone and neighbour atoms it previously occluded, so the
      SASA, CMS and SC of every REMAINING residue change for reasons that have nothing
      to do with the question. The result would be believable numbers on a mutilated
      structure. ``is_polymer()`` is the authoritative test (the submit-time name list
      in ``run_polyv.py`` is only fast feedback and cannot be exhaustive).
    * **a requested resname that matches nothing is an error**, not a smaller
      selection -- exactly as an absent chain is. Otherwise a typo (``LIG`` for
      ``UNL``, ``HEME`` for ``HEM``) or a residue Rosetta already dropped under
      ``-ignore_unrecognized_res`` leaves the thing you meant to exclude silently in
      the measurement, with status OK. Matches are counted only among residues on the
      SCORED chains, since a match on a chain that was dropped anyway excluded nothing.
    """
    per_res = chain_ids(pose)
    drop_chain: list[int] = []
    drop_resname: list[int] = []
    matched: dict[str, int] = {name: 0 for name in exclude_resnames}
    polymer_hits: set[str] = set()
    for i in range(1, pose.total_residue() + 1):
        if per_res[i] not in keep_chains:
            drop_chain.append(i)
            continue
        name3 = pose.residue(i).name3().strip().upper()
        if name3 in exclude_resnames:
            matched[name3] += 1
            if pose.residue(i).is_polymer():
                polymer_hits.add(name3)
            drop_resname.append(i)

    if polymer_hits:
        raise ValueError(
            f"--exclude-resnames names polymer residue {','.join(sorted(polymer_hits))}"
            f"; the flag is for non-polymer heteroatoms (a ligand, ion, lipid or "
            f"cofactor). Deleting an interior polymer residue exposes backbone its "
            f"neighbours occluded, so every remaining residue's SASA, CMS and SC would "
            f"change for reasons unrelated to the interface."
        )
    unmatched = sorted(name for name, count in matched.items() if count == 0)
    if unmatched:
        raise ValueError(
            f"--exclude-resnames {','.join(unmatched)} matched no residue on chains "
            f"{','.join(sorted(keep_chains))}; a requested resname that is absent is "
            f"an error, not a smaller selection"
        )

    drop = sorted(drop_chain + drop_resname)
    if not drop:
        return 0, 0
    counts = (len(drop_chain), len(drop_resname))
    blocks: list[tuple[int, int]] = []
    start = prev = drop[0]
    for i in drop[1:]:
        if i == prev + 1:
            prev = i
        else:
            blocks.append((start, prev))
            start = prev = i
    blocks.append((start, prev))
    for first, last in reversed(blocks):
        pose.delete_residue_range_slow(first, last)
    return counts


def mutate_to_val(
    pose: Any, residues: list[int], keep_gly_pro: bool
) -> tuple[int, int]:
    """Mutate the given residues to VAL. Returns ``(n_val, n_kept_gly_pro)``.

    Both counts are recounted from the POSE afterwards, not from bookkeeping, so the
    completeness invariant in ``run_one`` checks what is actually there. With
    ``keep_gly_pro`` off, ``n_kept_gly_pro`` is forced to 0, so a GLY or PRO that
    survived for any other reason cannot be absorbed into the gap.
    """
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

    for i in residues:
        name3 = pose.residue(i).name3()
        if keep_gly_pro and name3 in ("GLY", "PRO"):
            continue
        if name3 == "VAL":
            continue
        MutateResidue(i, "VAL").apply(pose)
    n_val = sum(1 for i in residues if pose.residue(i).name3() == "VAL")
    n_kept = (
        sum(1 for i in residues if pose.residue(i).name3() in ("GLY", "PRO"))
        if keep_gly_pro
        else 0
    )
    return n_val, n_kept


def pack_binder(pose: Any, sfxn: Any, binder_res: list[int], extra_rotamers: str) -> None:
    """Repack ONLY the binder's side chains. Everything else is frozen."""
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        ExtraRotamersGeneric,
        IncludeCurrent,
        OperateOnResidueSubset,
        PreventRepackingRLT,
        RestrictToRepacking,
    )
    from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

    binder_sel = ResidueIndexSelector(",".join(str(i) for i in binder_res))

    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())  # identities are already set; never design
    tf.push_back(IncludeCurrent())
    if extra_rotamers != "none":
        ex = ExtraRotamersGeneric()
        ex.ex1(True)
        if extra_rotamers == "ex1ex2":
            ex.ex2(True)
        tf.push_back(ex)
    # flip_subset=True -> prevent repacking of everything that is NOT the binder.
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), binder_sel, True))

    packer = PackRotamersMover(sfxn)
    packer.task_factory(tf)
    packer.apply(pose)


def interface_metrics(
    pose: Any,
    sfxn: Any,
    binder_chains: list[str],
    target_chains: list[str],
    binder_res: set[int],
    target_res: set[int],
) -> dict[str, Any]:
    """InterfaceAnalyzerMover, configured explicitly -- no defaults are relied on."""
    from pyrosetta.rosetta.core.pose import DockingPartners
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

    partners = DockingPartners.docking_partners_from_string(
        f"{''.join(binder_chains)}_{''.join(target_chains)}"
    )
    iam = InterfaceAnalyzerMover(partners)
    iam.set_scorefunction(sfxn)
    # Rigid-body separation, NOT a repacked unbound state.
    _call_first(iam, ["set_pack_separated"], False, what="pack_separated = False")
    _call_first(iam, ["set_pack_input"], False, what="pack_input = False")
    _call_first(iam, ["set_calc_dSASA"], True, what="compute dSASA")
    _call_first(
        iam, ["set_compute_separated_sasa"], True, what="compute separated SASA"
    )
    _call_first(iam, ["set_compute_interface_sc"], True, what="compute interface SC")
    # Explicitly off: not wanted here, and they cost time.
    _call_first(iam, ["set_compute_packstat"], False, what="packstat off")
    _call_first(iam, ["set_calc_hbond_sasaE"], False, what="hbond-SASA energy off")
    _call_first(
        iam, ["set_compute_interface_energy"], False, what="interface energy off"
    )
    _call_first(
        iam,
        ["set_compute_interface_delta_hbond_unsat"],
        False,
        what="delta unsat hbonds off",
    )
    for optional, value in (("set_use_jobname", False), ("set_skip_reporting", True)):
        fn = getattr(iam, optional, None)
        if fn is not None:
            fn(value)

    iam.apply(pose)

    dsasa = float(iam.get_interface_delta_sasa())
    n_reported = int(iam.get_num_interface_residues())
    getter = getattr(iam, "get_interface_set", None)
    if getter is None:
        raise RuntimeError(
            "PyRosetta API mismatch: InterfaceAnalyzerMover has no get_interface_set(); "
            "per-side interface residue counts cannot be computed honestly."
        )
    iface = {int(i) for i in getter()}
    if len(iface) != n_reported:
        raise ValueError(
            f"interface residue count mismatch: get_num_interface_residues()="
            f"{n_reported} vs |get_interface_set()|={len(iface)}"
        )
    n_binder = len(iface & binder_res)
    n_target = len(iface & target_res)
    if n_binder + n_target != len(iface):
        raise ValueError(
            f"interface residues outside both sides: {len(iface) - n_binder - n_target}"
        )
    # Cross-check the chain partition we computed against Rosetta's own. The partner
    # string is "<binder>_<target>", so side1 is the binder (verified). A mismatch
    # means the mover split the pose differently than the chain IDs say, which would
    # make every per-side number wrong: that is an error, not a trimmed fit.
    side1, side2 = int(iam.get_side1_nres()), int(iam.get_side2_nres())
    if (n_binder, n_target) != (side1, side2):
        raise ValueError(
            f"interface side mismatch: chain-based (binder {n_binder}, target "
            f"{n_target}) vs InterfaceAnalyzer (side1 {side1}, side2 {side2})"
        )

    sc_raw = float(iam.get_all_data().sc_value)
    # SC is defined on (0, 1]. Rosetta reports 0 (or a negative sentinel) when the SC
    # calculation failed -- a tiny or degenerate interface. That is NA, never 0.
    sc = sc_raw if math.isfinite(sc_raw) and 0.0 < sc_raw <= 1.0 else None

    return {
        "dsasa": round(dsasa, 4),
        "dsasa_per_res": round(dsasa / len(iface), 4) if iface else None,
        "nres_int": len(iface),
        "nres_int_binder": n_binder,
        "nres_int_target": n_target,
        "sc": None if sc is None else round(sc, 4),
    }


def contact_molecular_surface(
    pose: Any, surface_res: list[int], partner_res: list[int], distance_weight: float
) -> float | None:
    """Rosetta's ContactMolecularSurface for `surface_res` weighted toward `partner_res`.

    ``selector1`` is the side whose molecular surface is MEASURED, ``selector2`` the
    side it is weighted toward. Verified against the filter's own RosettaScripts
    schema ("Calculates the contact molecular surface area on the target defined by
    the target_selector, and the area is weighted by the closest distance between the
    target and the binder") and by running the same pose through an XML
    ``<ContactMolecularSurface target_selector=... binder_selector=.../>``: the XML
    value equalled ``selector1=target`` to 4 decimals and differed from the swapped
    call, so target_selector == selector1.

    Returns None when the filter raises. In practice that means the two selections do
    not touch at all (Rosetta throws a message-less RuntimeError from
    ContactMolecularSurfaceFilter.cc rather than returning 0), so the value is NA, not
    0 -- read ``dsasa`` and ``nres_int``, which are an honest 0 in that case.
    """
    from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
    from pyrosetta.rosetta.protocols.simple_filters import ContactMolecularSurfaceFilter

    f = ContactMolecularSurfaceFilter()
    _call_first(
        f,
        ["selector1", "set_selector1"],
        ResidueIndexSelector(",".join(str(i) for i in surface_res)),
        what="the side whose surface is measured",
    )
    _call_first(
        f,
        ["selector2", "set_selector2"],
        ResidueIndexSelector(",".join(str(i) for i in partner_res)),
        what="the side the surface is weighted toward",
    )
    _call_first(
        f, ["distance_weight", "set_distance_weight"], float(distance_weight),
        what="distance_weight",
    )
    compute = getattr(f, "compute", None) or getattr(f, "report_sm", None)
    if compute is None:
        raise RuntimeError(
            "PyRosetta API mismatch: ContactMolecularSurfaceFilter has neither "
            "compute() nor report_sm()."
        )
    try:
        return float(compute(pose))
    except RuntimeError:
        return None


def run_one(
    src: Path,
    name: str,
    binder_chains: list[str],
    target_chains: list[str],
    keep_gly_pro: bool,
    extra_rotamers: str,
    distance_weight: float,
    exclude_resnames: set[str],
    out_dir: Path,
) -> tuple[dict[str, Any], str]:
    """One design -> (metrics, packed-pose path). Raises on anything unusable."""
    import pyrosetta

    t0 = time.perf_counter()
    # BEFORE anything stochastic, and per design: Rosetta's RNG state otherwise
    # carries across the designs of a task, so a score would depend on batch position.
    seed = design_seed(name)
    seed_rosetta_rng(seed)

    if not src.exists():
        raise FileNotFoundError(f"structure missing: {src}")

    pose = pyrosetta.pose_from_file(str(src))
    if pose.total_residue() == 0:
        raise ValueError(f"no residues read from {src.name}")

    per_res = chain_ids(pose)
    present = ordered_unique(per_res[1:])
    check_contiguous(per_res)
    for chain_id in binder_chains + target_chains:
        if chain_id not in present:
            raise ValueError(
                f"chain {chain_id} not in {src.name} (have: {','.join(present)})"
            )
    chains_found = ",".join(present)

    # Everything that is in neither selection goes, so the pose is exactly two sides.
    n_res_dropped_chain, n_res_dropped_resname = drop_unscored_residues(
        pose, set(binder_chains) | set(target_chains), exclude_resnames
    )
    per_res = chain_ids(pose)
    binder_res = [i for i in range(1, pose.total_residue() + 1) if per_res[i] in binder_chains]
    target_res = [i for i in range(1, pose.total_residue() + 1) if per_res[i] in target_chains]
    if not binder_res or not target_res:
        raise ValueError(
            f"empty side after chain selection (binder {len(binder_res)} residues, "
            f"target {len(target_res)} residues)"
        )

    n_res_binder = len(binder_res)
    n_res_target = len(target_res)
    # Distinct residue names on the binder BEFORE the mutation. 1 == a uniform
    # backbone whatever letter the generator used (poly-GLY, poly-ALA, poly-VAL all
    # read as 1); a handful or more == the pose arrived carrying a real sequence,
    # which polyv is about to erase.
    n_input_restypes = len({pose.residue(i).name3() for i in binder_res})

    n_mutated, n_kept_gly_pro = mutate_to_val(pose, binder_res, keep_gly_pro)
    # Completeness invariant, enforced on BOTH paths. With --keep-gly-pro off,
    # n_kept_gly_pro is 0 by construction, so this reduces to n_mutated ==
    # n_res_binder; with it on, the gap must be exactly the GLY/PRO still there and
    # nothing else, so a residue that failed to mutate cannot hide in it.
    if n_mutated + n_kept_gly_pro != n_res_binder:
        unmutated = sorted(
            {
                pose.residue(i).name3()
                for i in binder_res
                if pose.residue(i).name3() != "VAL"
                and not (keep_gly_pro and pose.residue(i).name3() in ("GLY", "PRO"))
            }
        )
        raise ValueError(
            f"polyV mutation incomplete: {n_mutated} VAL + {n_kept_gly_pro} kept "
            f"GLY/PRO != {n_res_binder} binder residues; still present: "
            f"{','.join(unmutated) or '(none)'}"
        )

    sfxn = pyrosetta.create_score_function("ref2015")  # rotamer placement only
    pack_binder(pose, sfxn, binder_res, extra_rotamers)

    metrics = interface_metrics(
        pose.clone(),
        sfxn,
        binder_chains,
        target_chains,
        set(binder_res),
        set(target_res),
    )
    cms_target = contact_molecular_surface(pose, target_res, binder_res, distance_weight)
    cms_binder = contact_molecular_surface(pose, binder_res, target_res, distance_weight)
    metrics["cms_target"] = None if cms_target is None else round(cms_target, 4)
    metrics["cms_binder"] = None if cms_binder is None else round(cms_binder, 4)
    metrics.update(
        n_mutated=n_mutated,
        n_kept_gly_pro=n_kept_gly_pro,
        n_res_binder=n_res_binder,
        n_res_target=n_res_target,
        n_res_dropped_chain=n_res_dropped_chain,
        n_res_dropped_resname=n_res_dropped_resname,
        chains_found=chains_found,
        n_input_restypes=n_input_restypes,
        seed=seed,
    )

    out_pdb = out_dir / f"{name}.pdb"
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    pose.dump_pdb(str(out_pdb))

    metrics["seconds"] = round(time.perf_counter() - t0, 3)
    return metrics, str(out_pdb)


def write_result(out_dir: Path, name: str, status: str, path: str, metrics: dict) -> None:
    """One-row TSV. A metric that could not be computed is written as NA (an empty
    cell), never as 0 -- a dSASA or CMS of 0 is a real, alarming measurement."""
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / f"{name}.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(RESULT_COLUMNS)
        writer.writerow(
            [name, status, path]
            + ["" if metrics.get(c) is None else str(metrics[c]) for c in METRIC_COLUMNS]
        )


class Args(argparse.Namespace):
    task_file: Path
    binder_chains: str
    target_chains: str
    keep_gly_pro: bool
    extra_rotamers: str
    distance_weight: float
    exclude_resnames: str
    out_dir: Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Poly-valine interface geometry for a batch of complexes."
    )
    ap.add_argument(
        "--task-file",
        type=Path,
        required=True,
        help="Tab-separated 'name<TAB>structure' rows, one per design.",
    )
    ap.add_argument("--binder-chains", required=True, help="Comma-joined chain IDs.")
    ap.add_argument("--target-chains", required=True, help="Comma-joined chain IDs.")
    ap.add_argument("--keep-gly-pro", action="store_true")
    ap.add_argument(
        "--extra-rotamers", choices=("none", "ex1", "ex1ex2"), default="ex1"
    )
    ap.add_argument("--distance-weight", type=float, default=0.5)
    ap.add_argument("--exclude-resnames", default="", help="Comma-joined resnames.")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(namespace=Args())

    binder_chains = split_list(args.binder_chains)
    target_chains = split_list(args.target_chains)
    exclude_resnames = {r.upper() for r in split_list(args.exclude_resnames)}
    overlap = sorted(set(binder_chains) & set(target_chains))
    if overlap:
        raise SystemExit(
            f"chains {overlap} are in both --binder-chains and --target-chains; an "
            f"interface needs two disjoint sides."
        )
    if BLANK_CHAIN in binder_chains + target_chains:
        raise SystemExit(
            f"{BLANK_CHAIN!r} is polyv's rendering of a BLANK chain ID in "
            f"chains_found, not a selectable chain, and it is also the separator in "
            f"InterfaceAnalyzerMover's partner string. Give the structure real chain "
            f"IDs (see the chainsel tool) instead."
        )

    with open(args.task_file) as f:
        members = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    import pyrosetta

    pyrosetta.init(INIT_FLAGS)

    for name, src in members:
        try:
            metrics, path = run_one(
                Path(src),
                name,
                binder_chains,
                target_chains,
                args.keep_gly_pro,
                args.extra_rotamers,
                args.distance_weight,
                exclude_resnames,
                args.out_dir,
            )
            write_result(args.out_dir, name, "OK", path, metrics)
            print(
                f"{name}: cms_target={metrics['cms_target']} "
                f"cms_binder={metrics['cms_binder']} dsasa={metrics['dsasa']} "
                f"sc={metrics['sc']} nres_int={metrics['nres_int']} "
                f"seed={metrics['seed']} ({metrics['seconds']} s)"
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            # One line, and short: pybind errors list every overload.
            message = " ".join(str(e).split())
            # Rosetta's own exceptions are often just "File: <path>:<line>" with no
            # text at all; name the exception type so the status is diagnosable.
            if not message or message.startswith("File:"):
                message = f"{type(e).__name__} {message}".strip()
            message = message[:300]
            print(f"{name}: ERROR {message}", file=sys.stderr)
            (args.out_dir / f"{name}.pdb").unlink(missing_ok=True)
            write_result(args.out_dir, name, f"error: {message}", "", {})


if __name__ == "__main__":
    main()
