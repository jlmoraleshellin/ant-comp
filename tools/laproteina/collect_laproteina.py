#!/usr/bin/env python3
"""
Collect La-Proteina Complexa generation jobs into a child table, one row per binder.

A ``create`` tool, and a **root** create: there is no parent table, so the framework
iterates the design groups the run recorded in the sidecar (``root_designs``) and
this mints one row per backbone each job produced.

What a job leaves behind, after ``laproteina.sh`` has relocated it:

    <out_dir>/jobs/<group>/
        job_<job>_n_<len>_id_<idx>/
            job_<job>_n_<len>_id_<idx>.pdb      <- the binder/target complex
        rewards_<config_name>_<job>.csv          <- per-sample reward scores, if any

The **PDB files are the contract**; the rewards CSV is a bonus. A job that wrote
structures but no CSV still collects cleanly, with the reward columns absent rather
than the rows missing. That ordering is deliberate: the structures are what the next
step needs, and a missing score should not cost a design.

The collected structure is the **complex** (binder plus the target chains LPC was
conditioned on), not the binder alone. Anything that needs the binder by itself --
a self-consistency comparison, an interface metric -- must run ``chainsel`` first.
This mirrors ``bindcraft2``, whose ``path`` is also a complex.
"""

import re
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from prosapia.core import Collected, CollectCtx, CollectEach, DesignCtx

SAMPLES_DIRNAME = "jobs"

# job_<job>_n_<length>_id_<idx>[_<meta>]  -- LPC's own sample directory naming,
# from generate.py. The trailing meta tag is optional and free-form.
_SAMPLE_RE = re.compile(r"^job_(?P<job>\d+)_n_(?P<length>\d+)_id_(?P<idx>\d+)(?P<meta>.*)$")

# Reward columns worth lifting into the table when the CSV is there. Anything else
# in the CSV is ignored; this is a table, not an archive.
_REWARD_HINTS = ("reward", "score", "plddt", "ptm", "iptm", "pae", "rmsd")


def _run_meta(out_dir: Path) -> dict:
    """The sidecar the run wrote. Absent on a hand-made out_dir; treat as empty."""
    meta_path = out_dir / ".meta.json"
    if not meta_path.is_file():
        return {}
    import json

    try:
        return json.loads(meta_path.read_text())
    except (ValueError, OSError):
        return {}


def _samples_root(ctx: CollectCtx, meta: dict) -> Path:
    """Where the tasks parked their output. The run records it; we do not guess."""
    recorded = meta.get("samples_root")
    if recorded:
        p = Path(str(recorded))
        if p.is_dir():
            return p
    # Fall back to the conventional location under this collect's out_dir, so a
    # run_dir moved between filesystems still collects.
    return ctx.out_dir / SAMPLES_DIRNAME


def _read_rewards(group_dir: Path) -> pd.DataFrame | None:
    """The per-sample reward table, if the job wrote one."""
    csvs = sorted(group_dir.glob("rewards_*.csv"))
    if not csvs:
        return None
    frames = []
    for c in csvs:
        try:
            frames.append(pd.read_csv(c))
        except (ValueError, OSError) as e:
            print(f"  (unreadable rewards csv {c.name}: {e})")
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def _reward_row(rewards: pd.DataFrame | None, sample_dir: str) -> dict[str, Any]:
    """Reward columns for one sample, matched by whichever column names it."""
    if rewards is None or rewards.empty:
        return {}
    key_cols = [c for c in rewards.columns if c.lower() in ("sample", "name", "design", "dir", "sample_name")]
    for kc in key_cols:
        hit = rewards[rewards[kc].astype(str) == sample_dir]
        if not hit.empty:
            row = hit.iloc[0]
            out: dict[str, Any] = {}
            for col in rewards.columns:
                if col == kc:
                    continue
                if any(h in col.lower() for h in _REWARD_HINTS):
                    val = row[col]
                    if pd.notna(val):
                        out[str(col)] = val
            return out
    return {}


def collect_laproteina(ctx: CollectCtx) -> CollectEach:
    """Per-job LPC collector. The framework iterates the design groups the run
    recorded and this mints one child row per backbone the job produced. A job with
    no structures yields no rows. The framework stamps status/path/parent_name from
    each Collected."""
    meta = _run_meta(ctx.out_dir)
    samples_root = _samples_root(ctx, meta)
    task_name = meta.get("task_name", "(unknown)")
    print(f"Collecting LPC backbones for target '{task_name}' from {samples_root}")

    if not samples_root.is_dir():
        print(f"  WARNING: {samples_root} does not exist — nothing to collect")

    def one(d: DesignCtx) -> Iterable[Collected]:
        group_dir = samples_root / d.name
        if not group_dir.is_dir():
            print(f"{d.name}: no job dir, skipping")
            return

        sample_dirs = sorted(p for p in group_dir.iterdir() if p.is_dir())
        if not sample_dirs:
            print(f"{d.name}: no sample directories under {group_dir}")
            return

        rewards = _read_rewards(group_dir)
        n = 0
        for sdir in sample_dirs:
            m = _SAMPLE_RE.match(sdir.name)
            if m is None:
                # Not a sample directory (LPC also drops logs and intermediates
                # here). Skipping quietly would hide a naming change upstream, so
                # say it once per directory.
                print(f"{d.name}: {sdir.name} is not a sample dir, skipping")
                continue

            pdbs = sorted(sdir.glob("*.pdb"))
            if not pdbs:
                print(f"{d.name}: {sdir.name} has no pdb, skipping")
                continue
            if len(pdbs) > 1:
                print(f"{d.name}: {sdir.name} has {len(pdbs)} pdbs, taking {pdbs[0].name}")

            data: dict[str, Any] = {
                "job": int(m.group("job")),
                "length": int(m.group("length")),
                "sample": int(m.group("idx")),
                "sample_dir": sdir.name,
                "task_name": task_name,
            }
            data.update(_reward_row(rewards, sdir.name))

            # Strip LPC's redundant job_<id>_ prefix: the design group name already
            # carries the job, and the full dir name makes row keys unreadable.
            suffix = sdir.name[len(f"job_{m.group('job')}_"):]
            yield Collected(
                name=f"{d.name}_{suffix}",
                parent=d.name,
                path=str(pdbs[0]),
                status="OK",
                data=data,
            )
            n += 1

        print(f"{d.name}: OK ({n} backbone(s))")

    return one
