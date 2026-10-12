#!/usr/bin/env python3
"""
The per-array-task step of the dockq tool. It CALLS the published DockQ v2 package
(Mirabello & Wallner) and only selects, maps and reports its numbers: no part of Fnat,
iRMSD, LRMSD or the DockQ score is computed here. Premise, scope limits, traps and how
to read the columns: the tool's skill, .claude/skills/dockq/SKILL.md.

THE CHAIN MAPPING -- the contract this file must honour. DockQ maps model to reference
chains by SEQUENCE, and at its default allowed_mismatches=0 a designed binder matches
nothing. So this worker FIXES the binder explicitly (reference binder <- model binder,
as DockQ's `--mapping D*:D*` does) and lets DockQ's own group_chains() cluster ONLY the
target chains, at --allowed-mismatches. Every permutation DockQ's get_all_chain_maps()
yields is scored with run_on_all_native_interfaces(); one whose target pairs differ at
more than --allowed-mismatches is not a candidate. Every candidate must yield the SAME
set of native binder interfaces, or the design errors: that is the homo-oligomer
symmetry assumption, made checkable, and it keeps the binder sums comparable term for
term. DockQ's own objective is the SUM over ALL native interfaces, target-target ones
included; this worker instead picks the permutation with the highest sum over the
BINDER interfaces only (tie: higher full sum, then the first in DockQ's order), and
records what DockQ's own objective would have picked (`mapping_full_objective`,
`mapping_agrees`) and the dockq_best that choice gives (`dockq_best_full_objective`).
Because the kept mapping maximises the binder-interface SUM and not the best single
interface, dockq_best_full_objective may sit above or below dockq_best. Every
permutation's numbers are kept in `per_mapping`.

The default --residue-match number passes DockQ's own no_align=True: residues
correspond BY NUMBER and every position counts. It is guarded here: the residue-number
lists of every mapped chain pair must be IDENTICAL, or the design errors.
--residue-match sequence reproduces DockQ's default behaviour.

Normally run per array task by dockq.sh; run ``sapia run dockq`` rather than this.
"""

import argparse
import csv
import sys
from pathlib import Path

# CAPRI-style DockQ classes (the DockQ v2 README legend).
CLASS_BOUNDS = ((0.80, "high"), (0.49, "medium"), (0.23, "acceptable"))
ACCEPTABLE = 0.23

TRUST_COLUMNS = [
    "binder_chain",
    "target_chains",
    "ref_binder_chain",
    "ref_target_chains",
    "residue_match",
    "allowed_mismatches",
    "min_binder_identity",
    "min_native_contacts",
    "dockq_version",
    "ref_path",
    "model_chain_ids",
    "ref_chain_ids",
    "binder_n_res",
    "ref_binder_n_res",
    "binder_identity",
    "binder_n_res_compared",
    "target_max_mismatches",
    "target_max_mismatches_any_mapping",
    "mapping",
    "mapping_full_objective",
    "mapping_agrees",
    "n_mappings_tried",
    "per_mapping",
    "n_native_binder_interfaces",
    "n_interfaces_used",
]
SCIENCE_COLUMNS = [
    "dockq_best",
    "fnat_best",
    "irmsd_best",
    "lrmsd_best",
    "class_best",
    "best_ref_chain",
    "nat_best",
    "dockq_mean_natw",
    "n_interfaces_acceptable",
    "dockq_best_full_objective",
    "per_interface",
]

# DockQ functions known (2.1.3) to be functools.lru_cache'd on Biopython chain
# objects. Every cached callable in the module is cleared between designs anyway
# (clear_package_caches); this list only makes a renamed one visible in the log.
EXPECTED_CACHED = (
    "align_chains",
    "get_residue_distances",
    "get_aligned_residues",
    "subset_atoms",
    "run_on_chains",
    "list_atoms_per_residue",
)
METRIC_COLUMNS = TRUST_COLUMNS + SCIENCE_COLUMNS
RESULT_COLUMNS = ["name", "status", "path", *METRIC_COLUMNS]

# The echo of the question: kept on an error row, as in hisbond.
ECHO_COLUMNS = {
    "binder_chain",
    "target_chains",
    "ref_binder_chain",
    "ref_target_chains",
    "residue_match",
    "allowed_mismatches",
    "min_binder_identity",
    "min_native_contacts",
    "dockq_version",
}

_CIF_SUFFIXES = (".cif", ".mmcif", ".cif.gz", ".mmcif.gz")


class DockQCallError(RuntimeError):
    """The package refused (it calls sys.exit) or raised."""


def capri_class(score: float) -> str:
    for bound, label in CLASS_BOUNDS:
        if score >= bound:
            return label
    return "incorrect"


def call_package(what: str, fn, *args, **kwargs):
    """Call a DockQ function. DockQ reports several failures with sys.exit(), which is
    a BaseException and would escape `except Exception`: convert it into an error."""
    try:
        return fn(*args, **kwargs)
    except SystemExit as e:
        raise DockQCallError(f"DockQ {what} exited ({e.code!r}); see the task .err") from None


def load(dq, path: Path, chains: list[str]):
    """Read a structure with DockQ's own parsers, keeping only `chains`.

    For mmCIF this calls DockQ's MMCIFParser exactly as DockQ's load_PDB() falls back to
    it (auth_chains=True, no HETATM): load_PDB tries its PDB parser first, which on a
    Boltz mmCIF raises (measured: KeyError 0) and falls through to this same call.
    Calling it directly removes the dependence on that exception. Chain IDs are
    auth_asym_id, residue numbers auth_seq_id."""
    if not path.is_file():
        raise ValueError(f"structure missing: {path}")
    name = path.name.lower()
    if name.endswith(_CIF_SUFFIXES):
        import gzip

        from DockQ.parsers import MMCIFParser

        opener = gzip.open if name.endswith(".gz") else open
        with opener(path, "rt") as handle:
            model = MMCIFParser(QUIET=True).get_structure(
                "-", handle, chains=chains, parse_hetatms=False, auth_chains=True,
                model_number=0,
            )
        model.id = str(path)
        return model
    return dq.load_PDB(str(path), chains=chains)


def chain_info(structure, wanted: list[str], role: str, path: Path):
    found = [c.id for c in structure]
    missing = [c for c in wanted if c not in found]
    if missing:
        raise ValueError(
            f"{role} chain(s) {','.join(missing)} not found in {path.name} "
            f"(protein chains present: {','.join(found) or 'none'})"
        )
    for c in wanted:
        ch = structure[c]
        if ch.is_het is not None:
            raise ValueError(f"{role} chain {c} is a HETATM ligand chain, not a protein")
        if len(ch) == 0:
            raise ValueError(f"{role} chain {c} has no residues")
    return found


def numbering(chain) -> list[int]:
    return [int(r.id[1]) for r in chain]


def pair_identity_by_number(model_chain, ref_chain, label: str) -> tuple[int, int]:
    """(n_identical, n_positions) for a chain pair compared BY RESIDUE NUMBER. The
    numbering must be identical: no_align=True would otherwise align numbers with gaps
    and silently drop residues."""
    mnum, rnum = numbering(model_chain), numbering(ref_chain)
    if mnum != rnum:
        raise ValueError(
            f"{label}: residue numbering differs between model ({len(mnum)} res, "
            f"{mnum[0]}..{mnum[-1]}) and reference ({len(rnum)} res, {rnum[0]}..{rnum[-1]}); "
            f"--residue-match number needs identical numbering (an insertion/deletion "
            f"or a renumbered chain)"
        )
    same = sum(
        1 for a, b in zip(model_chain, ref_chain) if a.get_resname() == b.get_resname()
    )
    return same, len(rnum)


def pair_identity_by_sequence(dq, model_chain, ref_chain) -> tuple[int, int]:
    """(n_identical, n_ref_res) from DockQ's own alignment: the '|' positions are the
    ONLY residues DockQ compares in its default mode."""
    aln = dq.format_alignment(dq.align_chains(model_chain, ref_chain, use_numbering=False))
    return aln["matches"].count("|"), len(ref_chain)


def fmt_map(chain_map: dict[str, str], ref_order: list[str]) -> str:
    """`model>ref` pairs in reference-chain order, e.g. 'B>A;C>B;A>C;D>D'."""
    return ";".join(f"{chain_map[r]}>{r}" for r in ref_order)


def measure(dq, model_path: Path, ref_path: Path, a, metrics: dict) -> None:
    binder, targets = a.binder_chain, a.target_chains
    rbinder, rtargets = a.ref_binder_chain, a.ref_target_chains
    no_align = a.residue_match == "number"

    model = load(dq, model_path, targets + [binder])
    ref = load(dq, ref_path, rtargets + [rbinder])
    # Which chains the FILES hold (all protein chains, unfiltered): a trust fact.
    metrics["model_chain_ids"] = ",".join(c.id for c in load(dq, model_path, []))
    metrics["ref_chain_ids"] = ",".join(c.id for c in load(dq, ref_path, []))
    chain_info(model, targets + [binder], "model", model_path)
    chain_info(ref, rtargets + [rbinder], "reference", ref_path)

    metrics["binder_n_res"] = len(model[binder])
    metrics["ref_binder_n_res"] = len(ref[rbinder])

    # Binder: fixed by hand (sequence matching cannot admit a designed variant).
    if no_align:
        same, n = pair_identity_by_number(model[binder], ref[rbinder], "binder")
        metrics["binder_n_res_compared"] = n
    else:
        same, n = pair_identity_by_sequence(dq, model[binder], ref[rbinder])
        metrics["binder_n_res_compared"] = same
    identity = same / n
    metrics["binder_identity"] = identity
    if identity < a.min_binder_identity:
        raise ValueError(
            f"binder identity {identity:.3f} ({same}/{n}) is below --min-binder-identity "
            f"{a.min_binder_identity}: the model binder is not a variant of the "
            f"reference binder (wrong reference, or wrong chain letter)"
        )

    # Targets: DockQ's own sequence clustering, at --allowed-mismatches.
    clusters, reverse_map = call_package(
        "group_chains (target chains do not match the reference by sequence at "
        f"--allowed-mismatches {a.allowed_mismatches})",
        dq.group_chains, model, ref, targets, rtargets, a.allowed_mismatches,
    )
    maps = list(
        dq.get_all_chain_maps(clusters, {rbinder: binder}, reverse_map, targets, rtargets)
    )
    if not maps:
        raise ValueError("DockQ produced no candidate chain mapping for the target chains")
    metrics["n_mappings_tried"] = len(maps)

    rorder = rtargets + [rbinder]
    candidates: list[dict] = []
    per_mapping: list[str] = []
    # Written into metrics as it grows, so an error raised INSIDE the loop still
    # carries the permutations computed so far ('' = failed before the first one).
    metrics["per_mapping"] = ""
    any_mm = 0
    for idx, cm in enumerate(maps):
        if sorted(cm) != sorted(rorder) or len(set(cm.values())) != len(cm):
            raise ValueError(f"DockQ returned an incomplete chain map {cm}")
        # Target mismatches FOR THIS PERMUTATION; one that exceeds the allowance is
        # not a candidate (it is recorded in per_mapping as skipped).
        mm = 0
        for r in rtargets:
            m = cm[r]
            if no_align:
                same_t, n_t = pair_identity_by_number(model[m], ref[r], f"target {m}>{r}")
            else:
                same_t, n_t = pair_identity_by_sequence(dq, model[m], ref[r])
            mm = max(mm, n_t - same_t)
        any_mm = max(any_mm, mm)
        # ',' inside a mapping here, because ';' separates the per_mapping entries.
        mstr = fmt_map(cm, rorder).replace(";", ",")
        if mm > a.allowed_mismatches:
            per_mapping.append(f"{mstr}:target_mm={mm}:skipped")
            metrics["per_mapping"] = ";".join(per_mapping)
            continue
        result, full_sum = call_package(
            "run_on_all_native_interfaces",
            dq.run_on_all_native_interfaces, model, ref, chain_map=cm, no_align=no_align,
        )
        bkeys = frozenset(k for k in result if rbinder in k)
        bsum = sum(float(result[k]["DockQ"]) for k in bkeys)
        # Best binder DockQ over the interfaces that pass --min-native-contacts.
        bbest = max(
            (
                float(result[k]["DockQ"])
                for k in bkeys
                if int(result[k]["nat_total"]) >= a.min_native_contacts
            ),
            default=None,
        )
        candidates.append(
            {"idx": idx, "cm": cm, "result": result, "full_sum": float(full_sum),
             "bkeys": bkeys, "bsum": bsum, "bbest": bbest, "mm": mm}
        )
        per_mapping.append(
            f"{mstr}:target_mm={mm}:n_binder={len(bkeys)}"
            f":binder_sum={bsum:.6f}:full_sum={float(full_sum):.6f}"
        )
        metrics["per_mapping"] = ";".join(per_mapping)
    metrics["target_max_mismatches_any_mapping"] = any_mm
    if not candidates:
        raise ValueError(
            f"every target permutation differs from the reference at more than "
            f"--allowed-mismatches {a.allowed_mismatches} positions (worst {any_mm})"
        )
    # The symmetry assumption, made checkable: every candidate must score the SAME
    # set of native binder interfaces, so the binder sums compare equal numbers of
    # terms and no permutation can win on count.
    key_sets = {c["bkeys"] for c in candidates}
    if len(key_sets) != 1:
        raise ValueError(
            "target permutations disagree on which native binder interfaces exist: "
            + " vs ".join(",".join(sorted(s)) or "(none)" for s in key_sets)
        )

    # Kept: highest binder-interface sum; tie -> higher full sum; tie -> first in
    # DockQ's order. DockQ's own: first strict max of the all-interface sum.
    best = candidates[0]
    for c in candidates[1:]:
        if c["bsum"] > best["bsum"] or (
            c["bsum"] == best["bsum"] and c["full_sum"] > best["full_sum"]
        ):
            best = c
    full_best = candidates[0]
    for c in candidates[1:]:
        if c["full_sum"] > full_best["full_sum"]:
            full_best = c
    metrics["target_max_mismatches"] = best["mm"]

    cm, result = best["cm"], best["result"]
    metrics["mapping"] = fmt_map(cm, rorder)
    metrics["mapping_full_objective"] = fmt_map(full_best["cm"], rorder)
    metrics["mapping_agrees"] = best["idx"] == full_best["idx"]

    entries = []
    for key, info in result.items():
        if rbinder not in key:
            continue  # target-target interface: not this tool's question
        other = key.replace(rbinder, "", 1)
        if len(other) != 1 or other not in rtargets:
            raise ValueError(f"unexpected DockQ interface key {key!r}")
        for field in ("DockQ", "fnat", "iRMSD", "LRMSD", "nat_total"):
            if field not in info or info[field] is None:
                raise ValueError(f"DockQ gave no {field} for interface {key}")
        entries.append(
            {
                "ref": other,
                "model": cm[other],
                "dockq": float(info["DockQ"]),
                "fnat": float(info["fnat"]),
                "irmsd": float(info["iRMSD"]),
                "lrmsd": float(info["LRMSD"]),
                "nat": int(info["nat_total"]),
            }
        )
    metrics["n_native_binder_interfaces"] = len(entries)
    if not entries:
        raise ValueError(
            f"the reference binder {rbinder} makes no native contact (DockQ: no residue "
            f"pair within 5 A) with any of {','.join(rtargets)}: nothing to reproduce"
        )

    entries.sort(key=lambda e: e["ref"])
    for e in entries:
        e["used"] = e["nat"] >= a.min_native_contacts
    used = [e for e in entries if e["used"]]
    metrics["n_interfaces_used"] = len(used)
    if not used:
        raise ValueError(
            f"no native binder interface has >= --min-native-contacts "
            f"{a.min_native_contacts} native contacts (counts: "
            + ",".join(f"{e['ref']}={e['nat']}" for e in entries)
            + ")"
        )
    # Highest DockQ among the USED interfaces. An exact tie is broken ALPHABETICALLY
    # by reference chain -- an arbitrary but deterministic rule; the letter carries
    # no meaning.
    top = max(used, key=lambda e: (e["dockq"], -ord(e["ref"][0])))
    total_nat = sum(e["nat"] for e in used)
    metrics.update(
        {
            "dockq_best": top["dockq"],
            "fnat_best": top["fnat"],
            "irmsd_best": top["irmsd"],
            "lrmsd_best": top["lrmsd"],
            "class_best": capri_class(top["dockq"]),
            "best_ref_chain": top["ref"],
            "nat_best": top["nat"],
            # Weighted by each interface's NATIVE contact count (DockQ's nat_total).
            "dockq_mean_natw": sum(e["dockq"] * e["nat"] for e in used) / total_nat,
            "n_interfaces_acceptable": sum(1 for e in used if e["dockq"] >= ACCEPTABLE),
            # What dockq_best would be under DockQ's own all-interface mapping. The
            # kept mapping maximises the binder SUM, not the max, so this can be
            # above or below dockq_best.
            "dockq_best_full_objective": full_best["bbest"],
            # Full repr precision, so the collector re-derives dockq_best, the mean,
            # the acceptable count and the tie-break from exactly these numbers.
            # used=0 marks an interface below --min-native-contacts (kept, not used).
            "per_interface": ";".join(
                f"{e['ref']}:dockq={e['dockq']!r}:fnat={e['fnat']!r}"
                f":irmsd={e['irmsd']!r}:lrmsd={e['lrmsd']!r}"
                f":class={capri_class(e['dockq'])}:nat={e['nat']}:model={e['model']}"
                f":used={int(e['used'])}"
                for e in entries
            ),
        }
    )


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, float):
        # repr round-trips exactly, so the collector's arithmetic checks against
        # per_interface compare the worker's own numbers, not a re-rounding of them.
        return repr(value)
    return str(value)


def clear_package_caches(dq) -> list[str]:
    """Clear every functools cache in DockQ's module. DockQ lru_caches functions on
    Biopython chain objects, whose hash/equality use full_id -- and the structure id
    is the same placeholder for every file -- so without this a later design in the
    same process could be handed an earlier design's cached result."""
    cleared = []
    for attr in dir(dq):
        fn = getattr(dq, attr, None)
        clear = getattr(fn, "cache_clear", None)
        if callable(clear):
            clear()
            cleared.append(attr)
    return cleared


def write_result(out_dir: Path, name: str, status: str, path: str, metrics: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / f"{name}.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(RESULT_COLUMNS)
        w.writerow([name, status, path, *(_fmt(metrics.get(c)) for c in METRIC_COLUMNS)])


def split_list(value: str) -> list[str]:
    return [t.strip() for t in value.split(",") if t.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description="DockQ of a predicted binder vs a reference complex")
    ap.add_argument("--task-file", type=Path, required=True,
                    help="'name<TAB>model<TAB>reference' rows; either path may be empty.")
    ap.add_argument("--binder-chain", required=True)
    ap.add_argument("--target-chains", required=True, help="Comma-joined.")
    ap.add_argument("--ref-binder-chain", required=True)
    ap.add_argument("--ref-target-chains", required=True, help="Comma-joined.")
    ap.add_argument("--allowed-mismatches", type=int, default=0)
    ap.add_argument("--residue-match", choices=("number", "sequence"), default="number")
    ap.add_argument("--min-binder-identity", type=float, default=0.5)
    ap.add_argument("--min-native-contacts", type=int, default=1)
    ap.add_argument("--input-column", default="the input structure column")
    ap.add_argument("--ref-col", default="the reference structure column")
    ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args()
    a.target_chains = split_list(a.target_chains)
    a.ref_target_chains = split_list(a.ref_target_chains)
    a.binder_chain = a.binder_chain.strip()
    a.ref_binder_chain = a.ref_binder_chain.strip()
    if not a.target_chains or len(a.target_chains) != len(a.ref_target_chains):
        raise ValueError("--target-chains and --ref-target-chains must be non-empty and equal in count")
    if a.binder_chain in a.target_chains or a.ref_binder_chain in a.ref_target_chains:
        raise ValueError("a binder chain is also listed as a target chain")
    if a.allowed_mismatches < 0 or not 0.0 <= a.min_binder_identity <= 1.0:
        raise ValueError("--allowed-mismatches must be >= 0 and --min-binder-identity in [0, 1]")
    if a.min_native_contacts < 1:
        raise ValueError("--min-native-contacts must be >= 1")

    echo = {
        "binder_chain": a.binder_chain,
        "target_chains": ",".join(a.target_chains),
        "ref_binder_chain": a.ref_binder_chain,
        "ref_target_chains": ",".join(a.ref_target_chains),
        "residue_match": a.residue_match,
        "allowed_mismatches": a.allowed_mismatches,
        "min_binder_identity": a.min_binder_identity,
        "min_native_contacts": a.min_native_contacts,
    }

    with open(a.task_file) as f:
        members = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    # The package import is per-run, but a failure is per-ROW data, not a crash.
    dq = None
    import_error = None
    try:
        from DockQ import DockQ as dq
    except Exception as e:  # noqa: BLE001
        import_error = f"could not import the DockQ package: {type(e).__name__}: {e}"
    # The version lookup is metadata, not the measurement: its failure must not turn
    # good rows into errors, so it has its own guard and an explicit 'unknown'.
    from importlib.metadata import PackageNotFoundError, version

    echo["dockq_version"] = "unknown"
    for dist in ("DockQ", "dockq"):
        try:
            echo["dockq_version"] = version(dist)
            break
        except PackageNotFoundError:
            print(f"no installed distribution named {dist!r}", file=sys.stderr)
    if dq is not None:
        cleared = clear_package_caches(dq)
        absent = [n for n in EXPECTED_CACHED if n not in cleared]
        print(
            f"DockQ {echo['dockq_version']}: caches cleared between designs: "
            f"{','.join(cleared) or 'none'}"
            + (f" (expected but not found: {','.join(absent)})" if absent else "")
        )

    for member in members:
        member = member + [""] * (3 - len(member))
        name, src, ref = member[:3]
        metrics = dict(echo)
        if dq is not None:
            clear_package_caches(dq)  # no design may see another design's cache
        try:
            if import_error:
                raise RuntimeError(import_error)
            if not src.strip():
                raise ValueError(
                    f"could not resolve a model structure from column {a.input_column!r} "
                    f"(own row, or absent from this table and every ancestor)"
                )
            if not ref.strip():
                raise ValueError(
                    f"could not resolve a reference structure from column {a.ref_col!r} "
                    f"(own row, or absent from this table and every ancestor)"
                )
            metrics["ref_path"] = ref
            if Path(src).resolve() == Path(ref).resolve():
                raise ValueError("model and reference are the same file")
            measure(dq, Path(src), Path(ref), a, metrics)
            write_result(a.out_dir, name, "OK", src, metrics)
            print(
                f"{name}: mapping={metrics['mapping']} agrees={metrics['mapping_agrees']} "
                f"best={metrics['dockq_best']:.3f} ({metrics['class_best']}, "
                f"{metrics['best_ref_chain']}) | {metrics['per_interface']}"
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            msg = " ".join(str(e).split())
            print(f"{name}: ERROR {msg}", file=sys.stderr)
            kept = {k: v for k, v in metrics.items() if k in ECHO_COLUMNS or k in TRUST_COLUMNS}
            write_result(a.out_dir, name, f"error: {msg}", "", kept)


if __name__ == "__main__":
    main()
