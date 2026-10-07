import sys, itertools, numpy as np, gemmi

st = gemmi.read_structure(sys.argv[1]); st.setup_entities()
st.remove_ligands_and_waters(); st.remove_hydrogens(); st.remove_alternative_conformations()
m = st[0]

PATCHES = {"20-25":range(20,26), "31-33":range(31,34), "65-67":range(65,68),
           "72-79":range(72,80), "86-92":range(86,93), "135-140":range(135,141),
           "143-149":range(143,150)}
EXPOSED = {  # rel-SASA >= 0.25 from the previous run
 "20-25":[20,21,23,24,25], "31-33":[31,32,33], "65-67":[65,67],
 "72-79":[72], "86-92":[86,87,88,89,90,91], "135-140":[135,137,138,140],
 "143-149":[144,145,147]}

def atoms(chain, nums):
    out=[]
    for r in m[chain]:
        if r.seqid.num in nums:
            for a in r: out.append([a.pos.x,a.pos.y,a.pos.z])
    return np.array(out)

def mindist(X, Y):
    return min(np.linalg.norm(X[:,None,:]-Y[None,:,:], axis=2).min(), 999) if len(X) and len(Y) else 999

print("="*74)
print("1. Are the EXPOSED parts of these patches spatially contiguous? (chain A)")
print("="*74)
print("   min heavy-atom distance between exposed patch atoms, in Angstrom")
print("   (<8 A = touching/contiguous surface;  <15 A = reachable by one binder)\n")
labels = list(PATCHES)
A = {L: atoms("A", EXPOSED[L]) for L in labels}
hdr = "         " + "".join(f"{L:>9}" for L in labels); print(hdr)
for L1 in labels:
    row = f"{L1:>9}"
    for L2 in labels:
        row += "        -" if L1==L2 else f"{mindist(A[L1],A[L2]):9.1f}"
    print(row)

print("\n" + "="*74)
print("2. Distance from each patch to the NEIGHBOURING protomer (composite epitope?)")
print("="*74)
allB = atoms("B", range(1,200)); allC = atoms("C", range(1,200))
print(f"   {'patch':10} {'->chainB':>10} {'->chainC':>10}   interpretation")
for L in labels:
    dB, dC = mindist(A[L], allB), mindist(A[L], allC)
    near = [n for n,d in (("B",dB),("C",dC)) if d < 8]
    interp = f"at interface with {'/'.join(near)}" if near else "single-protomer face"
    print(f"   {L:10} {dB:10.1f} {dC:10.1f}   {interp}")

print("\n" + "="*74)
print("3. LPC's built-in TNF hotspots: A113 and C73")
print("="*74)
seq = {r.seqid.num: gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
       for r in m["A"].get_polymer()}
print(f"   human 113 = {seq.get(113)}   human 73 = {seq.get(73)}")
print("   mouse 113 = L (human A113 -> mouse L111? see subs list: 111 A>L)")
print("   mouse  73 = ABSENT (deletion)")
d113 = mindist(atoms("A",[113]), allC)
print(f"   A113 to chain C: {d113:.1f} A  -> {'composite epitope' if d113<8 else 'not at that interface'}")
for L in labels:
    d = mindist(A[L], atoms("A",[113]))
    if d < 15: print(f"   patch {L} is {d:.1f} A from A113")
