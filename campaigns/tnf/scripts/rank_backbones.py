"""Rank LPC backbones by how well they engage the intended epitope.

There are no model scores to rank by -- the campaign ran with reward_model=null,
so `sapia collect` had nothing to collect. Everything here is measured from the
structure instead, and it ranks ENGAGEMENT, not affinity. Nothing below is
evidence that a design binds; that needs Boltz.

Score components, in priority order:
  hotspots_hit  how many of the six steered residues the binder actually touches
  composite     does it straddle two protomers (the receptor-blocking geometry)
  contacts      interface size, with a penalty for sprawl -- a binder touching
                70+ target residues is draped over the surface, not docked
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HOTSPOTS = {
    "f1": [("A", 32), ("A", 33), ("A", 145), ("C", 87), ("C", 91), ("C", 92)],
    "f2": [("A", 87), ("A", 88), ("A", 90), ("A", 92), ("A", 131), ("A", 135)],
    "f3": [("A", 113), ("A", 115), ("C", 72), ("C", 75), ("C", 77), ("C", 97)],
    "f4": [("A", 24), ("A", 65), ("A", 67), ("A", 138), ("A", 140), ("A", 141)],
}
CUTOFF = 5.0
IDEAL_CONTACTS = 35  # focused docking; well below E5's unguided ~75


def parse(pdb):
    res, order = {}, []
    with open(pdb) as fh:
        for ln in fh:
            if not ln.startswith(("ATOM  ", "HETATM")):
                continue
            ch = ln[21]
            try:
                num = int(ln[22:26])
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                continue
            if ch not in order:
                order.append(ch)
            res.setdefault((ch, num), []).append(xyz)
    return {k: np.array(v) for k, v in res.items()}, order


def score_table(run_dir, lbl, limit=None):
    df = pd.read_csv(Path(run_dir) / f"table0_{lbl}.tsv", sep="\t")
    df = df[df.laproteina_status == "OK"]
    if limit:
        df = df.head(limit)
    hs = set(HOTSPOTS[lbl])
    out = []
    for _, r in df.iterrows():
        p = Path(str(r.laproteina_path))
        if not p.is_file():
            continue
        res, order = parse(p)
        binder = order[-1]
        B = np.vstack([v for (c, _), v in res.items() if c == binder])
        touched, chains = set(), set()
        for (c, num), P in res.items():
            if c == binder:
                continue
            if np.linalg.norm(B[:, None, :] - P[None, :, :], axis=2).min() < CUTOFF:
                touched.add((c, num)); chains.add(c)
        hit = len(touched & hs)
        comp = len(chains) > 1
        sprawl = abs(len(touched) - IDEAL_CONTACTS) / IDEAL_CONTACTS
        out.append({
            "name": r["name"], "epitope": lbl,
            "hotspots": hit, "contacts": len(touched),
            "composite": "yes" if comp else "no",
            "binder_len": int(r.laproteina_binder_length),
            "score": hit * 10 + (5 if comp else 0) - sprawl * 5,
            "path": str(p),
        })
    return out


if __name__ == "__main__":
    run_dir = sys.argv[1]
    labels = sys.argv[2].split(",") if len(sys.argv) > 2 else ["f1", "f2"]
    rows = []
    for l in labels:
        rows += score_table(run_dir, l)
    d = pd.DataFrame(rows).sort_values("score", ascending=False)
    print(f"scored {len(d)} designs\n")
    print(d.head(12)[["name", "epitope", "hotspots", "contacts", "composite", "binder_len"]].to_string(index=False))
    print("\nPATHS:")
    for p in d.head(6).path:
        print(p)
