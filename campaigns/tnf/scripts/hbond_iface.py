"""Binder positions H-bonding to the target: the set worth holding fixed.

Geometric criterion only (no hydrogens in these models): donor/acceptor heavy
atoms (N,O) within 3.5 A across the interface. Cruder than DSSP/HBPLUS but the
right shape of answer, and it needs no extra dependency in the workstation.

Positions are reported 1-indexed within the binder chain, which is what
ProteinMPNN's fixed-positions language expects (every chain renumbered 1..L).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

HB = 3.5
def parse(p):
    res, order = {}, []
    for ln in open(p):
        if not ln.startswith(("ATOM  ", "HETATM")): continue
        ch, an, rn = ln[21], ln[12:16].strip(), ln[17:20].strip()
        try:
            num = int(ln[22:26]); xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError: continue
        if ch not in order: order.append(ch)
        res.setdefault((ch, num), {"n": rn, "a": [], "p": []})
        res[(ch, num)]["a"].append(an); res[(ch, num)]["p"].append(xyz)
    return res, order

def hbond_positions(pdb):
    res, order = parse(pdb)
    b = order[-1]
    bk = sorted([k for k in res if k[0] == b], key=lambda k: k[1])
    idx = {k: i + 1 for i, k in enumerate(bk)}          # 1..L within the chain
    T = [(k, v) for k, v in res.items() if k[0] != b]
    TP = np.vstack([np.array(v["p"])[[i for i, a in enumerate(v["a"]) if a[0] in "NO"]]
                    for _, v in T if any(a[0] in "NO" for a in v["a"])])
    out = []
    for k in bk:
        v = res[k]
        P = np.array(v["p"])[[i for i, a in enumerate(v["a"]) if a[0] in "NO"]]
        if not len(P): continue
        if np.linalg.norm(P[:, None, :] - TP[None, :, :], axis=2).min() < HB:
            out.append((idx[k], v["n"]))
    return out, len(bk)

if __name__ == "__main__":
    run, labels = sys.argv[1], sys.argv[2].split(",")
    rows = []
    for l in labels:
        d = pd.read_csv(Path(run) / f"table0_{l}.tsv", sep="\t")
        d = d[d.laproteina_status == "OK"].head(int(sys.argv[3]) if len(sys.argv) > 3 else 40)
        for _, r in d.iterrows():
            p = Path(str(r.laproteina_path))
            if not p.is_file(): continue
            pos, L = hbond_positions(p)
            rows.append({"epitope": l, "name": r["name"], "binder_len": L,
                         "n_hbond": len(pos), "has_his": any(n == "HIS" for _, n in pos),
                         "positions": ",".join(str(i) for i, _ in pos)})
    df = pd.DataFrame(rows)
    print(f"analysed {len(df)} designs\n")
    print(df.groupby("epitope").agg(designs=("name","count"), mean_hbond=("n_hbond","mean"),
          with_his_hbond=("has_his","sum")).round(2).to_string())
    print("\n=== designs whose H-bond set includes a HIS (pH-switch starting points) ===")
    h = df[df.has_his].sort_values("n_hbond", ascending=False)
    print(f"{len(h)} of {len(df)}\n")
    print(h.head(10)[["name","epitope","binder_len","n_hbond","positions"]].to_string(index=False))
