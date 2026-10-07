import sys, numpy as np, gemmi
st=gemmi.read_structure(sys.argv[1]); st.setup_entities(); st.remove_ligands_and_waters()
m=st[0]
one=lambda r: gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
HS=[("A",32),("A",33),("A",145),("C",87),("C",91),("C",92)]
pts={}
for ch,n in HS:
    for r in m[ch]:
        if r.seqid.num==n:
            pts[(ch,n)]=(np.array([[a.pos.x,a.pos.y,a.pos.z] for a in r]), one(r))
print("chosen hotspots:")
for (ch,n),(P,aa) in pts.items(): print(f"   {ch}{n:<4} {aa}")
ks=list(pts)
print("\npairwise min heavy-atom distance (A):")
print("        "+"".join(f"{ch}{n:>6}" for ch,n in ks))
mx=0
for k1 in ks:
    row=f"{k1[0]}{k1[1]:<5}"
    for k2 in ks:
        if k1==k2: row+="     -"
        else:
            d=np.linalg.norm(pts[k1][0][:,None,:]-pts[k2][0][None,:,:],axis=2).min()
            mx=max(mx,d); row+=f"{d:6.1f}"
    print(row)
allp=np.vstack([P for P,_ in pts.values()])
span=np.linalg.norm(allp[:,None,:]-allp[None,:,:],axis=2).max()
print(f"\ntotal footprint span = {span:.1f} A   (60-110 aa binder covers ~20-30 A)")
cen=allp.mean(0)
print(f"max distance from centroid = {np.linalg.norm(allp-cen,axis=1).max():.1f} A")
