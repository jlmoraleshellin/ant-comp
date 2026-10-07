import sys, gemmi
st = gemmi.read_structure(sys.argv[1])
st.setup_entities()
st.remove_ligands_and_waters()
st.remove_hydrogens()
st.remove_alternative_conformations()
st.remove_empty_chains()
st.write_pdb(sys.argv[2])
m = st[0]
print("wrote", sys.argv[2])
for ch in m:
    p = ch.get_polymer()
    if len(p):
        r = [x.seqid.num for x in p]
        print(f"   chain {ch.name}: {len(p)} res  {min(r)}-{max(r)}")
