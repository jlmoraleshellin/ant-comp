"""Find binder histidines at the TNF interface, across every collected table.

For the pH switch (bind at 7.4, release at 6.0) a histidine only helps if it is
AT the interface: His is neutral at 7.4 and protonates at 6.0, and the release
comes from that new positive charge being unfavourable. A His buried in the
core, or out on a surface loop, does nothing.

It also matters what the His faces. Protonated His+ next to target Arg/Lys is
repulsive -> binding weakens at low pH, which is what we want. Next to Asp/Glu
it forms a salt bridge -> binding STRENGTHENS at low pH, the opposite. So each
interface His is classified by its nearest charged target residue.

NOTE these are LPC's own all-atom sequences. atomium will redesign them, so this
measures what the generator offers as a starting point, not the final answer.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

CUTOFF = 5.0
NEAR = 8.0


def parse(pdb):
    res, order = {}, []
    with open(pdb) as fh:
        for ln in fh:
            if not ln.startswith(("ATOM  ", "HETATM")):
                continue
            ch, name = ln[21], ln[17:20].strip()
            try:
                num = int(ln[22:26])
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                continue
            if ch not in order:
                order.append(ch)
            res.setdefault((ch, num), {"name": name, "xyz": []})["xyz"].append(xyz)
    for v in res.values():
        v["xyz"] = np.array(v["xyz"])
    return res, order


def analyse(run_dir, labels):
    rows = []
    for lbl in labels:
        tsv = Path(run_dir) / f"table0_{lbl}.tsv"
        if not tsv.is_file():
            continue
        df = pd.read_csv(tsv, sep="\t")
        df = df[df.laproteina_status == "OK"]
        for _, r in df.iterrows():
            p = Path(str(r.laproteina_path))
            if not p.is_file():
                continue
            res, order = parse(p)
            binder = order[-1]
            tgt = {k: v for k, v in res.items() if k[0] != binder}
            T = np.vstack([v["xyz"] for v in tgt.values()])
            his = [(k, v) for k, v in res.items() if k[0] == binder and v["name"] == "HIS"]
            if not his:
                rows.append({"epitope": lbl, "name": r["name"], "his_total": 0,
                             "his_iface": 0, "basic": 0, "acidic": 0})
                continue
            n_if, nb, na = 0, 0, 0
            for k, v in his:
                if np.linalg.norm(v["xyz"][:, None, :] - T[None, :, :], axis=2).min() >= CUTOFF:
                    continue
                n_if += 1
                for tk, tv in tgt.items():
                    if tv["name"] not in ("ARG", "LYS", "ASP", "GLU"):
                        continue
                    if np.linalg.norm(v["xyz"][:, None, :] - tv["xyz"][None, :, :], axis=2).min() < NEAR:
                        if tv["name"] in ("ARG", "LYS"):
                            nb += 1
                        else:
                            na += 1
            rows.append({"epitope": lbl, "name": r["name"],
                         "his_total": len(his), "his_iface": n_if,
                         "basic": nb, "acidic": na})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    d = analyse(sys.argv[1], sys.argv[2].split(","))
    print(f"analysed {len(d)} designs\n")
    print("=== histidines per epitope ===")
    g = d.groupby("epitope").agg(
        designs=("name", "count"),
        mean_his=("his_total", "mean"),
        with_iface_his=("his_iface", lambda s: (s > 0).sum()),
        mean_iface_his=("his_iface", "mean"),
    )
    g["pct_with_iface_his"] = (100 * g.with_iface_his / g.designs).round(0)
    print(g.round(2).to_string())

    cand = d[(d.his_iface > 0) & (d.basic > d.acidic)].copy()
    cand["net_basic"] = cand.basic - cand.acidic
    cand = cand.sort_values(["his_iface", "net_basic"], ascending=False)
    print(f"\n=== pH-SWITCH CANDIDATES: interface His facing a NET BASIC pocket ===")
    print(f"({len(cand)} of {len(d)} designs)\n")
    print(cand.head(15)[["name", "epitope", "his_iface", "basic", "acidic", "net_basic"]].to_string(index=False))
