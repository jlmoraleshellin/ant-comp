#!/usr/bin/env python3
"""
Collect dockq results into the table.

Reads the per-design <name>.tsv files dockq_worker.py wrote under
<run_dir>/<table>/dockq[_<label>]/ and merges them back as <prefix>_status,
<prefix>_path (the MODEL structure that was scored) and one <prefix>_<field> column
per metric. Every number is DockQ v2's own; this module only types, gates and
cross-checks them. How to read them: .claude/skills/dockq/SKILL.md.

  TRUST COLUMNS -- read these before any of the others:
    binder_chain / target_chains / ref_binder_chain / ref_target_chains
                          echo of the chain flags.
    residue_match         'number' (by residue number, all positions) or 'sequence'
                          (DockQ's alignment, substituted positions dropped).
    allowed_mismatches    DockQ's allowed_mismatches, applied to the target chains.
    min_binder_identity   the identity floor the binder had to clear.
    min_native_contacts   interfaces below this many native contacts are excluded
                          from best/mean/acceptable and from
                          dockq_best_full_objective (1 = DockQ's own rule).
    dockq_version         the DockQ package version ('unknown' if the installed
                          metadata could not be read; the numbers are unaffected).
    ref_path              the REFERENCE file this row was scored against.
    model_chain_ids / ref_chain_ids
                          every protein chain in each file, in file order.
    binder_n_res / ref_binder_n_res
                          residues of the binder chain in model and reference.
    binder_identity       fraction of binder positions identical to the reference.
    binder_n_res_compared binder residues DockQ actually compared.
    target_max_mismatches positions differing in the worst target pair OF THE KEPT
                          MAPPING (0 expected).
    target_max_mismatches_any_mapping
                          the same over EVERY permutation DockQ produced, including
                          ones skipped for exceeding --allowed-mismatches.
    mapping               the chain mapping used, model>reference, reference order.
    mapping_full_objective  the mapping DockQ's OWN objective (sum over ALL native
                          interfaces, target-target included) would have picked.
    mapping_agrees        True when the two are the same mapping.
    n_mappings_tried      target permutations DockQ produced (6 for a trimer).
    per_mapping           every permutation (empty: the row failed before the first
                          permutation was checked; a shorter list on an error row:
                          the permutations computed before the failure): <map with ','>:target_mm=
                          :n_binder=:binder_sum=:full_sum=, or :skipped when it exceeded
                          --allowed-mismatches. ';'-joined.
    n_native_binder_interfaces
                          reference target chains the reference binder contacts.
    n_interfaces_used     of those, the ones with >= min_native_contacts contacts.

  SCIENCE:
    dockq_best            highest DockQ over the USED native binder interfaces.
    fnat_best / irmsd_best / lrmsd_best
                          Fnat, iRMSD (A), LRMSD (A) of THAT SAME interface.
    class_best            CAPRI class of dockq_best.
    best_ref_chain        REFERENCE target chain of that interface (exact ties go
                          to the alphabetically first chain).
    nat_best              native contacts of that interface.
    dockq_mean_natw       mean DockQ over native binder interfaces, weighted by each
                          interface's native contact count (DockQ's nat_total).
    n_interfaces_acceptable  native binder interfaces with DockQ >= 0.23.
    dockq_best_full_objective  dockq_best under DockQ's OWN mapping choice; it can be
                          above or below dockq_best (the kept mapping maximises the
                          binder SUM, not the max).
    per_interface         ';'-joined, reference-chain order:
                          <ref>:dockq=:fnat=:irmsd=:lrmsd=:class=:nat=:model=<chain>:used=0|1
                          (full repr precision; used=0 = below min_native_contacts)

NA (empty), never 0, marks a value that was not measured: every 'error: ...' and
'missing' row has all SCIENCE columns NA. A native interface the model fails to
reproduce is NOT NA: DockQ scores it (Fnat 0.0, a small DockQ) and that value stays.

An OK row is CROSS-CHECKED before it is admitted: every per_interface entry must parse
into its eight fields with used == (nat >= min_native_contacts); over the USED entries,
dockq_best must be the largest DockQ, best_ref_chain and nat_best that entry's,
dockq_mean_natw = sum(dockq_i * nat_i) / sum(nat_i) and n_interfaces_acceptable the
count >= 0.23 (floats to 1e-6); the entry count must equal n_native_binder_interfaces
and the used count n_interfaces_used. A breach is an error row, never a half-trusted OK.

Statuses: 'OK'; 'error: ...'; 'missing' when the design has no TSV.
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

STR_COLUMNS = [
    "binder_chain",
    "target_chains",
    "ref_binder_chain",
    "ref_target_chains",
    "residue_match",
    "dockq_version",
    "ref_path",
    "model_chain_ids",
    "ref_chain_ids",
    "mapping",
    "mapping_full_objective",
    "per_mapping",
    "class_best",
    "best_ref_chain",
    "per_interface",
]
BOOL_COLUMNS = ["mapping_agrees"]
INT_COLUMNS = [
    "allowed_mismatches",
    "binder_n_res",
    "ref_binder_n_res",
    "binder_n_res_compared",
    "target_max_mismatches",
    "target_max_mismatches_any_mapping",
    "n_mappings_tried",
    "n_native_binder_interfaces",
    "n_interfaces_used",
    "min_native_contacts",
    "nat_best",
    "n_interfaces_acceptable",
]
FLOAT_COLUMNS = [
    "min_binder_identity",
    "binder_identity",
    "dockq_best",
    "fnat_best",
    "irmsd_best",
    "lrmsd_best",
    "dockq_mean_natw",
    "dockq_best_full_objective",
]
RESULT_COLUMNS = STR_COLUMNS + BOOL_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS

# The MEASUREMENT. Named explicitly: every one of these is emptied on an error row.
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

# What an OK row genuinely MUST carry: the science, plus the mapping and identity
# facts without which the science cannot be trusted. Descriptive metadata
# (dockq_version, chain-id lists, echoes) is not in the gate: its absence must never
# turn a good measurement into an error.
REQUIRED_ON_OK = SCIENCE_COLUMNS + [
    "mapping",
    "mapping_full_objective",
    "mapping_agrees",
    "binder_identity",
    "binder_n_res",
    "ref_binder_n_res",
    "binder_n_res_compared",
    "target_max_mismatches",
    "n_native_binder_interfaces",
    "n_interfaces_used",
    "min_native_contacts",
    "ref_path",
]

CLASSES = {"incorrect", "acceptable", "medium", "high"}
ENTRY_KEYS = ["dockq", "fnat", "irmsd", "lrmsd", "class", "nat", "model", "used"]
TOL = 1e-6


def _empty() -> dict[str, Any]:
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _as_bool(value: Any) -> bool:
    text = str(value).strip().lower()
    if text in ("true", "1"):
        return True
    if text in ("false", "0"):
        return False
    raise ValueError(f"{value!r} is not a boolean")


def parse_per_interface(text: str) -> list[dict[str, Any]]:
    """Each entry back into its reference chain plus eight typed key=value fields.
    Raises ValueError on any malformed entry."""
    entries = []
    for raw in text.split(";"):
        parts = raw.split(":")
        if len(parts) != 1 + len(ENTRY_KEYS):
            raise ValueError(f"per_interface entry {raw!r} does not have 9 ':' fields")
        ref, fields = parts[0], parts[1:]
        if len(ref) != 1:
            raise ValueError(f"per_interface entry {raw!r} has no single-letter chain")
        kv = {}
        for key, field in zip(ENTRY_KEYS, fields):
            k, sep, v = field.partition("=")
            if k != key or not sep or not v:
                raise ValueError(f"per_interface entry {raw!r}: expected {key}=<value>")
            kv[key] = v
        entry = {
            "ref": ref,
            "dockq": float(kv["dockq"]),
            "fnat": float(kv["fnat"]),
            "irmsd": float(kv["irmsd"]),
            "lrmsd": float(kv["lrmsd"]),
            "class": kv["class"],
            "nat": int(kv["nat"]),
            "model": kv["model"],
            "used": kv["used"],
        }
        if entry["used"] not in ("0", "1"):
            raise ValueError(f"per_interface entry {raw!r}: used must be 0 or 1")
        entry["used"] = entry["used"] == "1"
        if entry["class"] not in CLASSES:
            raise ValueError(f"per_interface entry {raw!r}: unknown class")
        if entry["nat"] < 1:
            raise ValueError(f"per_interface entry {raw!r}: a native interface has nat >= 1")
        entries.append(entry)
    return entries


def check_ok_row(data: dict[str, Any]) -> list[str]:
    """Arithmetic the collected columns must satisfy against per_interface."""
    try:
        entries = parse_per_interface(str(data["per_interface"]))
    except ValueError as e:
        return [str(e)]
    problems = []
    if len(entries) != data["n_native_binder_interfaces"]:
        problems.append(
            f"per_interface has {len(entries)} entries, "
            f"n_native_binder_interfaces={data['n_native_binder_interfaces']}"
        )
    bad_flag = [e["ref"] for e in entries if e["used"] != (e["nat"] >= data["min_native_contacts"])]
    if bad_flag:
        problems.append(f"used flag disagrees with min_native_contacts on {','.join(bad_flag)}")
    used = [e for e in entries if e["used"]]
    if len(used) != data["n_interfaces_used"]:
        problems.append(f"{len(used)} used entries, n_interfaces_used={data['n_interfaces_used']}")
    if not used:
        return problems + ["no used interface on an OK row"]
    top = max(used, key=lambda e: (e["dockq"], -ord(e["ref"])))
    if top["nat"] != data["nat_best"]:
        problems.append(f"nat_best {data['nat_best']} != entry {top['nat']}")
    if abs(data["dockq_best"] - top["dockq"]) > TOL:
        problems.append(f"dockq_best {data['dockq_best']} != max entry {top['dockq']}")
    if data["best_ref_chain"] != top["ref"]:
        problems.append(f"best_ref_chain {data['best_ref_chain']} != entry {top['ref']}")
    total_nat = sum(e["nat"] for e in used)
    weighted = sum(e["dockq"] * e["nat"] for e in used)
    # dockq_mean_natw * sum(nat) == sum(dockq_i * nat_i), compared per unit of
    # sum(nat) so the 1e-6 tolerance applies to the mean itself.
    if abs(data["dockq_mean_natw"] - weighted / total_nat) > TOL:
        problems.append(
            f"dockq_mean_natw {data['dockq_mean_natw']} != "
            f"sum(dockq*nat)/sum(nat) {weighted / total_nat}"
        )
    acceptable = sum(1 for e in used if e["dockq"] >= 0.23)
    if acceptable != data["n_interfaces_acceptable"]:
        problems.append(
            f"n_interfaces_acceptable {data['n_interfaces_acceptable']} != {acceptable}"
        )
    if data["class_best"] not in CLASSES:
        problems.append(f"unknown class {data['class_best']!r}")
    if not 0.0 <= data["dockq_best"] <= 1.0:
        problems.append(f"dockq_best {data['dockq_best']} outside [0, 1]")
    return problems


def collect_dockq(ctx: CollectCtx) -> CollectEach:
    """Per-design dockq collector. Variants (another reference column, another
    residue-match mode) on one table are distinguished by --dir-label."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return
        result_df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row = result_df.iloc[0]
        status = str(row["status"])
        data = _empty()
        try:
            for col in STR_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS:
                if (value := _present(row, col)) is not None:
                    if col in INT_COLUMNS:
                        data[col] = int(value)
                    elif col in FLOAT_COLUMNS:
                        data[col] = float(value)
                    else:
                        data[col] = str(value)
            for col in BOOL_COLUMNS:
                if (value := _present(row, col)) is not None:
                    data[col] = _as_bool(value)
        except ValueError as e:
            yield Collected(status=f"error: unreadable cell ({e})", path="", data=_empty())
            return

        if status != "OK":
            for col in SCIENCE_COLUMNS:
                data[col] = None  # a 0 here would read as a measurement
            yield Collected(status=status, path="", data=data)
            return

        missing = [c for c in REQUIRED_ON_OK if data.get(c) is None]
        if missing:
            yield Collected(
                status=f"error: OK row missing {','.join(missing)}", path="", data=_empty()
            )
            return
        problems = check_ok_row(data)
        if problems:
            yield Collected(status=f"error: {'; '.join(problems)}", path="", data=_empty())
            return

        path = _present(row, "path")
        yield Collected(status=status, path="" if path is None else str(path), data=data)

    return one
