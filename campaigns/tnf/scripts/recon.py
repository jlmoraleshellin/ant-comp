"""TNF-alpha human/mouse conservation + surface-exposure recon."""
import sys, math, json
import numpy as np, gemmi

PROBE = 1.4
VDW = {"C":1.70,"N":1.55,"O":1.52,"S":1.80,"H":1.20,"SE":1.90}
# Tien et al. 2013 theoretical max SASA
MAXSASA = {"A":129,"R":274,"N":195,"D":193,"C":167,"E":223,"Q":225,"G":104,"H":224,
           "I":197,"L":201,"K":236,"M":224,"F":240,"P":159,"S":155,"T":172,"W":285,
           "Y":263,"V":174}

def load(path):
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    return st

def sphere_pts(n=200):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2*i/n)
    theta = math.pi * (1 + 5**0.5) * i
    return np.stack([np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)], 1)

SP = sphere_pts(200)

def sasa(st, chains):
    """Shrake-Rupley over the given chains together; returns {(chain,resnum): sasa}."""
    pos, rad, key = [], [], []
    for ch in st[0]:
        if ch.name not in chains: continue
        for r in ch:
            for a in r:
                pos.append([a.pos.x, a.pos.y, a.pos.z])
                rad.append(VDW.get(a.element.name.upper(), 1.70))
                key.append((ch.name, r.seqid.num))
    pos = np.array(pos); rad = np.array(rad) + PROBE
    out = {}
    rmax = rad.max()
    for i in range(len(pos)):
        d = np.linalg.norm(pos - pos[i], axis=1)
        nb = np.where((d > 1e-6) & (d < rad[i] + rmax))[0]
        pts = pos[i] + SP * rad[i]
        acc = np.ones(len(pts), bool)
        for j in nb:
            acc &= np.linalg.norm(pts - pos[j], axis=1) >= rad[j]
        a = 4*math.pi*rad[i]**2 * acc.sum()/len(pts)
        out[key[i]] = out.get(key[i], 0.0) + a
    return out

hu, mo = load(sys.argv[1]), load(sys.argv[2])

def seqmap(st):
    return {r.seqid.num: gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
            for r in st[0]["A"].get_polymer()}
H, M = seqmap(hu), seqmap(mo)

print("Computing SASA (trimer context)...", file=sys.stderr)
hs_tri = sasa(hu, {"A","B","C"})
ms_tri = sasa(mo, {"A","B","C"})

rel = {}
for n, aa in H.items():
    s = hs_tri.get(("A", n), 0.0)
    rel[n] = s / MAXSASA.get(aa, 200)

# superposition, CA chain A
def ca(st):
    return {r.seqid.num: r.sole_atom("CA") for r in st[0]["A"].get_polymer() if r.find_atom("CA","*")}
hca, mca = ca(hu), ca(mo)
common = sorted(set(hca) & set(mca))
P = gemmi.Position
sup = gemmi.superpose_positions([hca[n].pos for n in common], [mca[n].pos for n in common])
print(f"\nBackbone superposition (chain A, {len(common)} CA): RMSD = {sup.rmsd:.2f} A")

PATCHES = {"20-25":range(20,26), "31-33":range(31,34), "65-67":range(65,68),
           "72-79":range(72,80), "86-92":range(86,93), "135-140":range(135,141),
           "143-149":range(143,150)}

print("\n" + "="*78)
print("PATCH ANALYSIS  (human numbering; rel.SASA in trimer; * = differs in mouse)")
print("="*78)
summary = []
for label, rng in PATCHES.items():
    rows, nd, nexp, ndexp = [], 0, 0, 0
    for n in rng:
        h = H.get(n); m = M.get(n)
        rs = rel.get(n, float("nan"))
        if h is None: tag, m_s = "absent(H)", "-"
        elif m is None: tag, m_s = "DELETED", "-"
        else: tag, m_s = ("same" if h == m else "DIFF"), m
        exposed = rs >= 0.25
        if tag in ("DIFF","DELETED"): nd += 1
        if exposed: nexp += 1
        if exposed and tag in ("DIFF","DELETED"): ndexp += 1
        rows.append((n, h, m_s, rs, tag, exposed))
    tot = len(list(rng))
    summary.append((label, tot, nd, nexp, ndexp))
    print(f"\n--- patch {label} ---")
    for n, h, m, rs, tag, exp in rows:
        mark = "*" if tag in ("DIFF","DELETED") else " "
        e = "exposed" if exp else "buried "
        print(f"  {mark} {n:3}  human {h}  mouse {m}   relSASA {rs:5.2f} {e}  {tag}")

print("\n" + "="*78)
print("SUMMARY — ranked by cross-reactivity suitability")
print("="*78)
print(f"{'patch':10} {'size':>4} {'diffs':>6} {'exposed':>8} {'exposed+diff':>13}  verdict")
for label, tot, nd, nexp, ndexp in sorted(summary, key=lambda x: (x[4], x[2])):
    v = "CLEAN - fully conserved where it counts" if ndexp == 0 else \
        ("usable - 1 exposed difference" if ndexp == 1 else f"RISKY - {ndexp} exposed differences")
    print(f"{label:10} {tot:4} {nd:6} {nexp:8} {ndexp:13}  {v}")
