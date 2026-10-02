"""Seed table0 of a run_dir from the inherited BindCraft2 hEGFR complexes.

An input manifest, not a measurement: every prosapia table is normally minted by a
`create` tool, and this chain has no create step at its head (dimer_binder_plan.md
§2.2). Done once, explicitly, and recorded here.

Only the *_hEGFR.cif files at the top level are taken -- not renum/, whose only
difference is that the binder chain continues the target's numbering (B: 194.. rather
than 1..); the target numbering is identical in both, so nothing is gained, and the
per-chain numbering is what chainsel/rpxdock expect downstream.

Columns seeded:
    input_path      the complex on the runs volume (binder chain B + target chain A)
    source_file     the original filename, so the shortened name stays traceable
    ortholog        hEGFR
    binder_len_src  binder length parsed from the filename's l<N> field -- an
                    invariant ifacegeom_binder_len can be checked against, derived
                    from the name rather than from the structure
"""
import re
import sys
from pathlib import Path

import pandas as pd
from prosapia.core import DataManager

run_dir = Path(sys.argv[1])
src_dir = Path("/runs/inputs/bc2_output_for_anthony")

files = sorted(p for p in src_dir.glob("*_hEGFR.cif") if p.is_file())
if not files:
    raise SystemExit(f"no *_hEGFR.cif directly under {src_dir}")

rows = {}
for p in files:
    stem = p.name[: -len(".cif")]
    m = re.search(r"_(v\d+)_multitarget_l(\d+)_([0-9a-f]{8})", stem)
    if not m:
        raise SystemExit(f"cannot parse a name out of {p.name!r}")
    version, length, short_hash = m.groups()
    name = f"{version}_l{length}_{short_hash}"
    if name in rows:
        raise SystemExit(f"name collision on {name!r} -- widen the hash slice")
    rows[name] = {
        "input_path": str(p),
        "source_file": p.name,
        "ortholog": "hEGFR",
        "binder_len_src": int(length),
    }

df = pd.DataFrame.from_dict(rows, orient="index")
df.index.name = "name"

dm = DataManager(run_dir)
table = dm.rm.derive_new_table(None)  # root table -> table0
dm.rm.register_table(table)
dm.write_frame(table.table_name, df)

print(f"seeded {table.table_name} with {len(df)} rows from {src_dir}")
print(f"binder lengths from filenames: {sorted(df['binder_len_src'].unique())}")
print(df.head(3).to_string())
