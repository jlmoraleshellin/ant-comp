"""Scan the TNF trimer for epitope patches that are exposed, conserved, and basic
(so a binder histidine protonating at pH 6 meets positive charge -> repulsion)."""
import sys, math, numpy as np, gemmi
PROBE=1.4; VDW={"C":1.70,"N":1.55,"O":1.52,"S":1.80}
MAX={"A":129,"R":274,"N":195,"D":193,"C":167,"E":223,"Q":225,"G":104,"H":224,"I":197,
     "L":201,"K":236,"M":224,"F":240,"P":159,"S":155,"T":172,"W":285,"Y":263,"V":174}
def sph(n=120):
    i=np.arange(n)+0.5; phi=np.arccos(1-2*i/n); th=math.pi*(1+5**0.5)*i
    return np.stack([np.cos(th)*np.sin(phi),np.sin(th)*np.sin(phi),np.cos(phi)],1)
SP=sph()
hu=gemmi.read_structure(sys.argv[1]); hu.setup_entities(); hu.remove_ligands_and_waters(); hu.remove_hydrogens(); hu.remove_alternative_conformations()
mo=gemmi.read_structure(sys.argv[2]); mo.setup_entities(); mo.remove_ligands_and_waters(); mo.remove_hydrogens(); mo.remove_alternative_conformations()
m=hu[0]; one=lambda r: gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
H={r.seqid.num:one(r) for r in m["A"].get_polymer()}; M={r.seqid.num:one(r) for r in mo[0]["A"].get_polymer()}
pos,rad,key=[],[],[]
for ch in m:
    for r in ch:
        for a in r:
            pos.append([a.pos.x,a.pos.y,a.pos.z]); rad.append(VDW.get(a.element.name.upper(),1.7)+PROBE); key.append((ch.name,r.seqid.num))
pos=np.array(pos); rad=np.array(rad); rmax=rad.max(); sasa={}
for i in range(len(pos)):
    d=np.linalg.norm(pos-pos[i],axis=1); nb=np.where((d>1e-6)&(d<rad[i]+rmax))[0]
    pts=pos[i]+SP*rad[i]; acc=np.ones(len(pts),bool)
    for j in nb: acc&=np.linalg.norm(pts-pos[j],axis=1)>=rad[j]
    sasa[key[i]]=sasa.get(key[i],0.0)+4*math.pi*rad[i]**2*acc.sum()/len(pts)
cen={}; 
for ch in m:
    for r in ch:
        P=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in r])
        if len(P): cen[(ch.name,r.seqid.num)]=P.mean(0)
# seeds: exposed + conserved
seeds=[k for k,c in cen.items() if sasa.get(k,0)/MAX.get(H.get(k[1],"A"),200)>=0.30 and H.get(k[1])==M.get(k[1])]
out=[]
for s in seeds:
    near=[k for k in cen if np.linalg.norm(cen[k]-cen[s])<11]
    if len(near)<5: continue
    chains={k[0] for k in near}
    cons=sum(1 for k in near if H.get(k[1])==M.get(k[1]))/len(near)
    rk=sum(1 for k in near if H.get(k[1]) in "RK"); de=sum(1 for k in near if H.get(k[1]) in "DE")
    expo=sum(1 for k in near if sasa.get(k,0)/MAX.get(H.get(k[1],"A"),200)>=0.25)
    out.append((rk-de, cons, len(chains), expo, s, sorted(near)))
out.sort(key=lambda x:(-x[0],-x[1],-x[3]))
print("Top conserved+exposed+BASIC patches (rk-de = net basic residues, want high):\n")
seen=set()
for net,cons,nch,expo,s,near in out:
    if any(np.linalg.norm(cen[s]-cen[t])<9 for t in seen): continue
    seen.add(s)
    lbl=", ".join(f"{c}{n}" for c,n in near if sasa.get((c,n),0)/MAX.get(H.get(n,"A"),200)>=0.30)
    print(f"  seed {s[0]}{s[1]} {H.get(s[1])}  net_basic {net:+d}  conserved {cons:.0%}  chains {nch}  exposed {expo}")
    print(f"     exposed members: {lbl}")
    if len(seen)>=6: break
