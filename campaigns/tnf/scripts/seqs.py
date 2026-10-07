import gemmi, sys

def chain_seq(path, chain="A"):
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    poly = st[0][chain].get_polymer()
    out = []
    for r in poly:
        out.append((r.seqid.num, gemmi.find_tabulated_residue(r.name).one_letter_code.upper()))
    return out

h = chain_seq(sys.argv[1]); m = chain_seq(sys.argv[2])
print("HUMAN 1TNF chain A:", len(h), "residues,", h[0][0], "-", h[-1][0])
print("".join(c for _, c in h))
print()
print("MOUSE 2TNF chain A:", len(m), "residues,", m[0][0], "-", m[-1][0])
print("".join(c for _, c in m))
print()
# gaps in numbering?
for name, s in (("human", h), ("mouse", m)):
    nums = [n for n, _ in s]
    gaps = [(nums[i], nums[i+1]) for i in range(len(nums)-1) if nums[i+1] != nums[i]+1]
    print(f"{name}: numbering gaps -> {gaps if gaps else 'none (contiguous)'}")
