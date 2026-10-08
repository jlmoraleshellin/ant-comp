"""Score candidate hotspot sets on the human TNF trimer."""
import sys, math, numpy as np, gemmi

PROBE=1.4
VDW={"C":1.70,"N":1.55,"O":1.52,"S":1.80,"H":1.20}
MAXSASA={"A":129,"R":274,"N":195,"D":193,"C":167,"E":223,"Q":225,"G":104,"H":224,
         "I":197,"L":201,"K":236,"M":224,"F":240,"P":159,"S":155,"T":172,"W":285,"Y":263,"V":174}
SETS={
 "1 (validated)": [("A",32),("A",33),("A",145),("C",87),("C",91),("C",92)],
 "2":             [("A",79),("A",92),("A",95),("A",134),("A",135),("B",146)],
 "3":             [("A",24),("A",65),("A",67),("A",138),("A",140),("A",141)],
 "4":             [("A",45),("A",47),("A",49),("A",83),("A",85),("A",131)],
}
def sph(n=120):
    i=np.arange(n)+0.5; phi=np.arccos(1-2*i/n); th=math.pi*(1+5**0.5)*i
    return np.stack([np.cos(th)*np.sin(phi),np.sin(th)*np.sin(phi),np.cos(phi)],1)
SP=sph()

hu=gemmi.read_structure(sys.argv[1]); hu.setup_entities(); hu.remove_ligands_and_waters(); hu.remove_hydrogens(); hu.remove_alternative_conformations()
mo=gemmi.read_structure(sys.argv[2]); mo.setup_entities(); mo.remove_ligands_and_waters(); mo.remove_hydrogens(); mo.remove_alternative_conformations()
m=hu[0]
one=lambda r: gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
H={r.seqid.num:one(r) for r in m["A"].get_polymer()}
M={r.seqid.num:one(r) for r in mo[0]["A"].get_polymer()}

pos,rad,key=[],[],[]
for ch in m:
    for r in ch:
        for a in r:
            pos.append([a.pos.x,a.pos.y,a.pos.z]); rad.append(VDW.get(a.element.name.upper(),1.7)+PROBE)
            key.append((ch.name,r.seqid.num))
pos=np.array(pos); rad=np.array(rad); rmax=rad.max()
sasa={}
for i in range(len(pos)):
    d=np.linalg.norm(pos-pos[i],axis=1)
    nb=np.where((d>1e-6)&(d<rad[i]+rmax))[0]
    pts=pos[i]+SP*rad[i]; acc=np.ones(len(pts),bool)
    for j in nb: acc &= np.linalg.norm(pts-pos[j],axis=1)>=rad[j]
    sasa[key[i]]=sasa.get(key[i],0.0)+4*math.pi*rad[i]**2*acc.sum()/len(pts)

def atoms(sel):
    return np.array([[a.pos.x,a.pos.y,a.pos.z] for ch in m for r in ch for a in r if (ch.name,r.seqid.num) in sel])

for label,S in SETS.items():
    print("="*72); print(f"SET {label}")
    exposed=0; diffs=[]; comp={}
    for c,n in S:
        aa=H.get(n,"?"); mm=M.get(n,"-")
        rs=sasa.get((c,n),0)/MAXSASA.get(aa,200)
        e = rs>=0.25
        exposed+=e
        if aa!=mm: diffs.append(f"{c}{n} {aa}>{mm}")
        comp[c]=comp.get(c,0)+1
        print(f"   {c}{n:<4} {aa}  relSASA {rs:4.2f} {'exposed' if e else 'BURIED '}  mouse {mm}{'  *DIFF' if aa!=mm else ''}")
    P=atoms(set(S))
    span=np.linalg.norm(P[:,None,:]-P[None,:,:],axis=2).max()
    # Arg/Lys within 8 A of the patch -> good partners for a pH-switch His
    rk=set()
    for ch in m:
        for r in ch:
            if one(r) not in "RK": continue
            Q=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in r])
            if np.linalg.norm(P[:,None,:]-Q[None,:,:],axis=2).min()<8: rk.add((ch.name,r.seqid.num,one(r)))
    de=set()
    for ch in m:
        for r in ch:
            if one(r) not in "DE": continue
            Q=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in r])
            if np.linalg.norm(P[:,None,:]-Q[None,:,:],axis=2).min()<8: de.add((ch.name,r.seqid.num,one(r)))
    print(f"   -> chains {dict(comp)} {'COMPOSITE' if len(comp)>1 else 'single protomer'}")
    print(f"   -> footprint {span:.1f} A   exposed {exposed}/{len(S)}   mouse diffs: {diffs or 'none'}")
    print(f"   -> Arg/Lys within 8A: {len(rk)}  {sorted(f'{c}{n}{a}' for c,n,a in rk)}")
    print(f"   -> Asp/Glu within 8A: {len(de)}")
    print(f"   -> pH-switch outlook: {'GOOD (basic patch -> His+ repulsion)' if len(rk)>=len(de) else 'POOR (acidic patch -> His+ would STRENGTHEN at low pH)'}")
