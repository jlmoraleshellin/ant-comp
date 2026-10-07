import gemmi, sys
HU, MO = sys.argv[1], sys.argv[2]

def full_seq(path):
    """SEQRES-level sequence (what's in the construct), plus modelled residues."""
    doc = gemmi.cif.read(path); b = doc.sole_block()
    seqres = [gemmi.find_tabulated_residue(gemmi.cif.as_string(r[0])).one_letter_code.upper()
              for r in b.find("_entity_poly_seq.", ["mon_id"])]
    st = gemmi.read_structure(path); st.setup_entities(); st.remove_ligands_and_waters()
    modelled = {r.seqid.num: gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                for r in st[0]["A"].get_polymer()}
    return seqres, modelled

for name, p in (("HUMAN 1TNF", HU), ("MOUSE 2TNF", MO)):
    seqres, mod = full_seq(p)
    print(f"{name}: SEQRES length = {len(seqres)}, modelled = {len(mod)}")
print()

_, hmod = full_seq(HU)
_, mmod = full_seq(MO)
print("Is mouse 73 in SEQRES but unmodelled, or a true deletion?")
mseqres, _ = full_seq(MO)
print("  mouse SEQRES len:", len(mseqres), "-> if 157 then res73 exists but is disordered")
print("  mouse modelled has 73?", 73 in mmod, " 72:", mmod.get(72), " 74:", mmod.get(74))
print("  human 72/73/74:", hmod.get(72), hmod.get(73), hmod.get(74))
print()

# direct numbering-based comparison over the shared numbering range
shared = sorted(set(hmod) & set(mmod))
ident = [n for n in shared if hmod[n] == mmod[n]]
print(f"Shared modelled positions: {len(shared)}  ({shared[0]}-{shared[-1]})")
print(f"Identical by number      : {len(ident)}  = {100*len(ident)/len(shared):.1f}%")
diffs = [(n, hmod[n], mmod[n]) for n in shared if hmod[n] != mmod[n]]
print(f"Differences              : {len(diffs)}")
print()
print("All substitutions (human -> mouse), by residue number:")
for i in range(0, len(diffs), 6):
    print("   " + "  ".join(f"{n:3}{a}>{b}" for n, a, b in diffs[i:i+6]))
