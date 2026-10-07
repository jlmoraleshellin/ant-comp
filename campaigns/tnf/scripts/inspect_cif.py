import sys, gemmi
for path in sys.argv[1:]:
    doc = gemmi.cif.read(path)
    b = doc.sole_block()
    def g(tag):
        v = b.find_value(tag)
        return gemmi.cif.as_string(v) if v else None
    print("="*70)
    print("FILE      :", path.split("/")[-1])
    print("PDB id    :", g("_entry.id"))
    print("title     :", (g("_struct.title") or "")[:90])
    print("method    :", g("_exptl.method"))
    print("resolution:", g("_refine.ls_d_res_high") or g("_em_3d_reconstruction.resolution"))
    # organism per entity
    print("-- entities --")
    ent = b.find("_entity.", ["id", "type", "pdbx_description"])
    for r in ent:
        print(f"   entity {r[0]}: {r[1]:10} {gemmi.cif.as_string(r[2])[:60]}")
    src = b.find("_entity_src_gen.", ["entity_id", "pdbx_gene_src_scientific_name"])
    for r in src:
        print(f"   SOURCE entity {r[0]}: {gemmi.cif.as_string(r[1])}")
    nat = b.find("_entity_src_nat.", ["entity_id", "pdbx_organism_scientific"])
    for r in nat:
        print(f"   SOURCE(nat) entity {r[0]}: {gemmi.cif.as_string(r[1])}")
    st = gemmi.read_structure(path)
    st.setup_entities()
    m = st[0]
    print("-- chains --")
    for ch in m:
        poly = ch.get_polymer()
        if len(poly):
            rs = [r.seqid.num for r in poly]
            print(f"   {ch.name}: {len(poly):4} res, range {min(rs)}-{max(rs)}")
