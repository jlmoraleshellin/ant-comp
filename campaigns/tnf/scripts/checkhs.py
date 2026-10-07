import sys, gemmi
A=[32,33,65,67,113,144,145,147]
Cc=[82,87,91,92,93,94,95,96,123,124,125,126]
one=lambda r: gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
res={}
for tag,p in (("human",sys.argv[1]),("mouse",sys.argv[2])):
    st=gemmi.read_structure(p); st.setup_entities(); st.remove_ligands_and_waters()
    res[tag]={ch.name:{r.seqid.num:one(r) for r in st[0][ch.name].get_polymer()} for ch in st[0]}
print(f"{'hotspot':9} {'human':>6} {'mouse':>6}   status")
bad=0
for ch,nums in (("A",A),("C",Cc)):
    for n in nums:
        h=res["human"][ch].get(n,"ABSENT"); m=res["mouse"][ch].get(n,"ABSENT")
        if h=="ABSENT" or m=="ABSENT": s="!! MISSING"; bad+=1
        elif h!=m: s="!! DIFFERS"; bad+=1
        else: s="ok"
        print(f"{ch}{n:<8} {h:>6} {m:>6}   {s}")
print(f"\n{len(A)+len(Cc)} hotspots, {bad} problems")
