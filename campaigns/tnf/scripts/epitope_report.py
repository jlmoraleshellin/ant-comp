"""Per-epitope report over the collected LPC tables.

Answers the question the campaign exists to answer: for each epitope, did the
designs actually land on the hotspots they were steered to?

Runs inside the Modal workstation, so stdlib + numpy/pandas only (no gemmi).
Samples N designs per table rather than all 200 -- the contact calculation is
O(binder_atoms x target_atoms) and a sample of 25 is plenty to compare epitopes.
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
    "f5": [],
}
CUTOFF = 5.0


def parse(pdb):
    """-> {(chain, resnum): Nx3 coords}, chain order of first appearance."""
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


def main(run_dir, n_sample=25):
    run_dir = Path(run_dir)
    rows = []
    for lbl in ("f1", "f2", "f3", "f4", "f5"):
        tsv = run_dir / f"table0_{lbl}.tsv"
        if not tsv.is_file():
            continue
        df = pd.read_csv(tsv, sep="\t")
        df = df[df.laproteina_status == "OK"]
        samp = df.sample(min(n_sample, len(df)), random_state=0)
        hs = set(HOTSPOTS[lbl])

        hit_counts, n_contacts, comps, blens = [], [], 0, []
        for _, r in samp.iterrows():
            p = Path(str(r.laproteina_path))
            if not p.is_file():
                continue
            res, order = parse(p)
            binder = order[-1]
            B = np.vstack([v for (c, _), v in res.items() if c == binder])
            blens.append(sum(1 for (c, _) in res if c == binder))
            touched, chains = set(), set()
            for (c, num), P in res.items():
                if c == binder:
                    continue
                if np.linalg.norm(B[:, None, :] - P[None, :, :], axis=2).min() < CUTOFF:
                    touched.add((c, num))
                    chains.add(c)
            n_contacts.append(len(touched))
            if hs:
                hit_counts.append(len(touched & hs))
            if len(chains) > 1:
                comps += 1

        rows.append(
            {
                "epitope": lbl,
                "n": len(samp),
                "binder_len": f"{int(np.mean(blens))}",
                "contacts": f"{np.mean(n_contacts):.0f}",
                "hotspots_hit": f"{np.mean(hit_counts):.1f}/6" if hit_counts else "n/a",
                "any_hotspot_%": f"{100*np.mean([h > 0 for h in hit_counts]):.0f}%" if hit_counts else "n/a",
                "composite_%": f"{100*comps/max(len(n_contacts),1):.0f}%",
            }
        )
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 25)
