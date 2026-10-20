#!/usr/bin/env python3
"""
Collect hisbond results into the table.

Reads the per-design <name>.tsv files hisbond_worker.py wrote under
<run_dir>/<table>/hisbond[_<label>]/ and merges them back as <prefix>_status,
<prefix>_path (the structure that was measured) and one <prefix>_<field> column per
metric. Every column and how to read it: ``.claude/skills/hisbond/SKILL.md``.

Both spellings are real, in the same tool: the centre-list columns are AMERICAN
(``centers_n``, ``centers_his_ok``, ``centers_arg``, ``centers_spec``) and the
per-centre counts BRITISH (``n_centres_hb_target``, ``centre<k>_idx1``, ...). Copy
names from ``RESULT_COLUMNS``: a wrong one is a KeyError.

EVERY INDEX HERE IS 1-BASED, along the binder chain in FILE ORDER. ``centre<k>_idx1``
is the k-th position the caller supplied with ``--centers``, and the same number is
the prefix of every ``;``-joined cell (``per_centre``, ``per_centre_nd1_target``,
``per_centre_ne2_target``, ``per_centre_nd1_binder``, ``per_centre_ne2_binder``,
``per_centre_role_target``, ``per_centre_role_binder``), so the join key between a
cell entry and a centre column is plainly ``centre<k>_idx1``.

Two column families, trust first:

  trust    echoes of the question (``binder_chain``, ``target_chains``,
           ``centers_arg`` -- the ``--centers`` template verbatim, e.g.
           ``{atomiumph_centre1_idx1},{atomiumph_centre2_idx1}`` --, ``centers_spec``
           -- this design's resolved positions, every one of them --,
           ``centers_lineage_depth`` -- 0 when every island came from the design's
           OWN table, otherwise the depth of the DEEPEST island --, ``env_radius``,
           ``min_dist``/``max_dist``/``max_angle``) and facts about the file
           (``centers_n``, ``centers_his_ok``, ``binder_chain_n_res``,
           ``expect_binder_len``, ``target_chain_n_res``, ``target_n_het_atoms``,
           ``chain_ids_found``, ``structure_has_hydrogens``, ``role_basis``,
           ``plddt_scale``). Read ``role_basis`` before any role column and
           ``plddt_scale`` before ``his_plddt_min``.
  science  contact counts and gates (merged, per nitrogen, and intra-binder), the
           donor/acceptor roles and ``switch_candidate`` gates, per-centre cells,
           ``centre<k>_idx1`` / ``_resnum`` (positions to pin), tautomer and charged
           environment counts, ``his_plddt_min``.

NA (empty) marks a value not measured: every 'error: ...' and 'missing' row has all
science columns NA and keeps only the echo and the trust facts established before the
failure. 0 is a measurement ("this His touches nothing"), so the two must stay apart.

Statuses: 'OK'; 'error: ...' (unresolvable structure, a ``--centers`` island that
resolved to nothing or to a non-integer for this design, a malformed position list,
a missing/unreadable file, an absent chain, a centre position past the end of the
chain, a position whose residue is not a histidine -- which is how an Asp/Glu
(acid_p / all_p) table errors on every row --, an incomplete imidazole, a binder
length disagreeing with --expect-binder-len-col, a crashed worker); 'missing' when the
design has no TSV.
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Keep MAX_CENTRE_COLUMNS in sync by hand with hisbond_worker.py (no shared imports).
# 8, raised from 4: each centre costs five columns, so this is already 40 columns on
# every table hisbond touches; past the cap centers_spec still lists every supplied
# position and centres_reported < centers_n makes the truncation visible.
MAX_CENTRE_COLUMNS = 8

STR_COLUMNS = [
    "binder_chain",
    "target_chains",
    "centers_arg",
    "centers_spec",
    "chain_ids_found",
    "plddt_scale",
    "role_basis",
    "per_centre",
    "per_centre_nd1_target",
    "per_centre_ne2_target",
    "per_centre_nd1_binder",
    "per_centre_ne2_binder",
    "per_centre_role_target",
    "per_centre_role_binder",
    *[f"centre{k}_tautomer" for k in range(1, MAX_CENTRE_COLUMNS + 1)],
]
BOOL_COLUMNS = [
    "structure_has_hydrogens",
    "any_target_hb_proxy",
    "any_target_hb_proxy_nd1",
    "any_target_hb_proxy_ne2",
    "any_switch_candidate",
    "any_switch_candidate_binder",
    "any_switch_candidate_clean_any",
]

INT_COLUMNS = [
    # 0 = every island came from the design's OWN table; otherwise the DEEPEST.
    "centers_lineage_depth",
    "centers_n",
    "centers_his_ok",
    "centres_reported",
    "centres_idx1_eq_resnum",
    *[
        f"centre{k}_{field}"
        for k in range(1, MAX_CENTRE_COLUMNS + 1)
        for field in ("idx1", "resnum", "n_acidic_env", "n_basic_env")
    ],
    "binder_chain_n_res",
    "expect_binder_len",
    "target_chain_n_res",
    "target_n_het_atoms",
    "n_centres_hb_target",
    "n_centres_hb_target_nd1",
    "n_centres_hb_target_ne2",
    "n_centres_acceptor_role_target",
    "n_centres_donor_role_target",
    "n_centres_switch_candidate",
    "n_centres_switch_candidate_clean",
    "n_contacts_incoherent_target",
    "n_contacts_ambiguous_target",
    "n_contacts_unclassified_target",
    "n_centres_no_ring_h",
    "n_centres_hip",
    "n_centres_foreign_h",
    "n_centres_hb_binder",
    # Intra-binder role block, mirroring the target-side counts.
    "n_centres_acceptor_role_binder",
    "n_centres_donor_role_binder",
    "n_centres_switch_candidate_binder",
    "n_centres_switch_candidate_clean_binder",
    "n_contacts_incoherent_binder",
    "n_contacts_ambiguous_binder",
    "n_contacts_unclassified_binder",
    "n_contacts_binder_entries",
    "n_centres_switch_candidate_clean_any",
    "n_centres_polar_dist_only_target",
    "n_centres_too_close_target",
    "n_centres_too_close_binder",
]
FLOAT_COLUMNS = [
    "min_dist",
    "max_dist",
    "max_angle",
    "env_radius",
    "min_dist_target",
    "min_dist_target_nd1",
    "min_dist_target_ne2",
    "min_dist_binder",
    "his_plddt_min",
]
RESULT_COLUMNS = STR_COLUMNS + BOOL_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS

# The measurements, emptied on every error row (NA is never 0); named, never sliced.
SCIENCE_COLUMNS = [
    "per_centre",
    "n_centres_hb_target",
    "any_target_hb_proxy",
    "min_dist_target",
    "n_centres_hb_target_nd1",
    "n_centres_hb_target_ne2",
    "any_target_hb_proxy_nd1",
    "any_target_hb_proxy_ne2",
    "min_dist_target_nd1",
    "min_dist_target_ne2",
    "per_centre_nd1_target",
    "per_centre_ne2_target",
    "per_centre_role_target",
    "n_centres_acceptor_role_target",
    "n_centres_donor_role_target",
    "n_centres_switch_candidate",
    "n_centres_switch_candidate_clean",
    "any_switch_candidate",
    "n_contacts_incoherent_target",
    "n_contacts_ambiguous_target",
    "n_contacts_unclassified_target",
    "n_centres_no_ring_h",
    "n_centres_hip",
    "n_centres_foreign_h",
    # role_basis is a measurement of the file too: emptied on an error row.
    "role_basis",
    "centres_reported",
    "centres_idx1_eq_resnum",
    *[
        f"centre{k}_{field}"
        for k in range(1, MAX_CENTRE_COLUMNS + 1)
        for field in ("idx1", "resnum", "tautomer", "n_acidic_env", "n_basic_env")
    ],
    "n_centres_hb_binder",
    "min_dist_binder",
    "per_centre_nd1_binder",
    "per_centre_ne2_binder",
    "per_centre_role_binder",
    "n_centres_acceptor_role_binder",
    "n_centres_donor_role_binder",
    "n_centres_switch_candidate_binder",
    "n_centres_switch_candidate_clean_binder",
    "any_switch_candidate_binder",
    "n_contacts_incoherent_binder",
    "n_contacts_ambiguous_binder",
    "n_contacts_unclassified_binder",
    "n_contacts_binder_entries",
    "n_centres_switch_candidate_clean_any",
    "any_switch_candidate_clean_any",
    "n_centres_polar_dist_only_target",
    "n_centres_too_close_target",
    "n_centres_too_close_binder",
    "his_plddt_min",
]

# What an OK row MUST carry; the optional ones are NA when nothing was measurable.
OPTIONAL_ON_OK = {
    "expect_binder_len",
    "min_dist_target",
    "min_dist_target_nd1",
    "min_dist_target_ne2",
    "min_dist_binder",
    "his_plddt_min",
    # Centre columns are NA past the row's centres; tautomer is NA without hydrogens.
    *[
        f"centre{k}_{field}"
        for k in range(1, MAX_CENTRE_COLUMNS + 1)
        for field in ("idx1", "resnum", "tautomer", "n_acidic_env", "n_basic_env")
    ],
}

# Leaves run before a column family existed error on collect: re-run with --force.


def _empty() -> dict[str, Any]:
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _as_bool(value: Any) -> bool:
    """A worker-written 'True'/'False' cell as a bool (never ``bool(str)``)."""
    text = str(value).strip().lower()
    if text in ("true", "1"):
        return True
    if text in ("false", "0"):
        return False
    raise ValueError(f"{value!r} is not a boolean")


def collect_hisbond(ctx: CollectCtx) -> CollectEach:
    """Per-design collector; variants on one table are told apart by ``--dir-label``."""

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

        # Echo and trust facts survive on error rows too; no science column is filled.
        for col in STR_COLUMNS + INT_COLUMNS + FLOAT_COLUMNS:
            if (value := _present(row, col)) is not None:
                try:
                    if col in INT_COLUMNS:
                        data[col] = int(value)
                    elif col in FLOAT_COLUMNS:
                        data[col] = float(value)
                    else:
                        data[col] = str(value)
                except ValueError:
                    yield Collected(
                        status=f"error: unreadable {col} cell {value!r}",
                        path="",
                        data=_empty(),
                    )
                    return
        for col in BOOL_COLUMNS:
            if (value := _present(row, col)) is not None:
                try:
                    data[col] = _as_bool(value)
                except ValueError:
                    yield Collected(
                        status=f"error: unreadable {col} cell {value!r}",
                        path="",
                        data=_empty(),
                    )
                    return

        if status != "OK":
            # Every science column stays NA: a 0 here would read as a measurement.
            for col in SCIENCE_COLUMNS:
                data[col] = None
            yield Collected(status=status, path="", data=data)
            return

        missing = [
            c
            for c in RESULT_COLUMNS
            if c not in OPTIONAL_ON_OK and data.get(c) is None
        ]
        if missing:
            # A hole is a contract breach: an error row, not a half-filled OK.
            yield Collected(
                status=f"error: OK row missing {','.join(missing)}",
                path="",
                data=_empty(),
            )
            return
        if data["centers_his_ok"] != data["centers_n"]:
            # Not producible by this worker, so a foreign TSV: never pass silently.
            yield Collected(
                status=(
                    f"error: OK row claims {data['centers_his_ok']} histidine "
                    f"centres of {data['centers_n']}"
                ),
                path="",
                data=_empty(),
            )
            return

        path = _present(row, "path")
        yield Collected(status=status, path="" if path is None else str(path), data=data)

    return one
