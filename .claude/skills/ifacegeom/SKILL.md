---
name: ifacegeom
description: >-
  Which residues of the binder and of the target form the interface, and where do the binder's N- and C-termini sit relative to it? Records the epitope residue list as a column that survives lineage, the target-side histidine count and terminus geometry. Load before composing an ifacegeom run or interpreting its columns. Traps: the vector selection method is the default for a reason, chain selection is easy to get wrong (the default binder is chain A), and every downstream tool parses the residue-label format.
---

# ifacegeom

**Custom tool** (`tools/ifacegeom/`, not bundled with prosapia). **`action: update`**: annotates the table it reads, in place; `-t` is required. CPU only (the manifest builder forces `gpus_per_task = 0`, no `-g 0`); on vib `SAPIA_ACTIVATE_IFACEGEOM` points at the shared `prosapia_python.sh`. **No default input column**: `-i` is effectively required (see Traps).

Measures, for each binder/target **complex**: (1) **which residues form the interface, on both sides** (the binder-side list is the *epitope footprint*; the target-side list carries the histidine count); (2) **where the termini sit relative to that interface**: a signed projection onto the binder-COM → interface-COM axis, plus two distances per terminus.

## Why it exists

When binders are inherited rather than designed here, the epitope is not known and has to be discovered before anything downstream can respect it (a rotamer graft needs the binder residues to transplant, a network design the positions it must not overlap, a tag requirement the C-terminus). It is a **tool** and not a side script because of the word *column*: an epitope list in a sidecar file does not survive the generation boundary into later tables, while columns ride lineage and can be re-filtered at any cutoff without recomputation. `cms` writes per-residue contributions to a **file** (`cms_path`), which cannot be filtered on or carried through lineage; `pyrosetta --interface` gives scalar `if_*` only, with no residue list or terminus geometry.

## Invocation

```bash
sapia run ifacegeom <run_dir> -t table0 -i input_path --binder-chains A --target-chains auto
sapia collect ifacegeom <run_dir> -t table0
```

Every input file must hold **both** binder and target (a monomer has no interface). Run it on whichever column holds the complex: a designed-pose column, or a `boltz_path` to measure a prediction's interface instead (label the second with `-l` so the columns do not collide).

| Flag | Default | What it does |
| --- | --- | --- |
| `--binder-chains` | `A` | Binder chain IDs, comma-joined. **Often wrong for the file at hand** (the binder is frequently the last chain): always pass it |
| `--target-chains` | `auto` | Target chain IDs; `auto` = every protein chain that is not a binder chain (right for one binder and one target of any chain count) |
| `--method` | `vector` | `vector` or `heavy`, see below |
| `--contact-cutoff` | `5.5` | Heavy-atom distance (Å) for the contact test |
| `--cb-dist-cut` | `11.0` | CB–CB distance beyond which a pair is not considered at all |
| `--vector-dist-cut` | `9.0` | CB–CB distance within which the pointing test applies |
| `--vector-angle-cut` | `75.0` | Max angle (°) between CA→CB and the direction to the partner's CB |
| `--designs-per-task` | `200` | Milliseconds per design; this amortises the container cold start |

### The two methods

`vector` (default) reproduces the **idea** of Rosetta's `InterGroupInterfaceByVector`: a residue is at the interface either because it touches the other side or because its side chain *points at* it.

```
nearby   any heavy atom of i within --contact-cutoff of any heavy atom of j  -> both i and j selected
vector   CB(i)-CB(j) within --vector-dist-cut AND the angle between i's CA->CB direction and
         the CB(i)->CB(j) direction is below --vector-angle-cut              -> i selected (and symmetrically j)
```

The pointing test is applied **per residue, not per pair**: each side's list stands on its own. `heavy` keeps the contact test alone. The vector test excludes residues that are merely nearby while facing away, exactly the rotamers that would be wasted work to graft back, so `vector` is the right default. **`vector` is usually, not always, a superset of `heavy`**: it applies the `--cb-dist-cut` prefilter first, so a residue whose CB is farther than that from the target is dropped even if a backbone atom is within the contact cutoff (typical: a brute-force 5.5 Å heavy-atom set of 26 residues vs 25 from `vector`, the missing residue having its backbone O 5.2 Å from the target but a CB–CB distance of 11.2 Å). Raise `--cb-dist-cut` if an edge residue matters. This is a reimplementation of the concept, **not a port of Rosetta's code** (no PyRosetta here; numbers are not bit-identical); the defaults are Rosetta's documented ones.

## Traps

- **`-i` is effectively required.** `default_input_column` is the `"not applicable"` sentinel (as in `cms`, `usalign`, `chainsel`): no honest default complex column exists and a wrong one fails silently, so the builder refuses the run and lists the table's `_path` columns.
- **A chain named but absent is an error for that design, never a smaller selection**: `--binder-chains A,B` against a file holding only `A` fails that row rather than measuring half a binder. After collect check `ifacegeom_binder_chains` and `ifacegeom_target_chains` (with `auto` they show what was used).
- **`warn: no interface residues` almost always means the wrong chains**, not a real non-contact. The row keeps `binder_len`/`target_len`/chain names (which makes the mistake obvious) and every geometry column is NA (no interface, no axis); the `warn:` prefix keeps it out of any `== "OK"` selection.
- **`--contact-cutoff` is an ATOM-to-atom distance.** The ~8–11 Å numbers usually quoted with `InterGroupInterfaceByVector` are its **CB–CB** cutoffs (`--vector-dist-cut`, `--cb-dist-cut`): `--contact-cutoff 8` does not "match Rosetta", it selects most of the binder.
- **The COMs are in each file's own frame.** `binder_com` / `iface_com` are for auditability, not comparison across designs; anything downstream in a different frame (a dock, a superposition) must recompute them.
- **Residue labels are `chain:resnum[icode]`** (`A:12`, not `A12`), the format every downstream consumer of `binder_res` parses. `ringfit --hotspots` uses the *other* convention (`A47`): do not paste one into the other.
- **The numbering is this file's numbering.** The epitope list is meaningful only against a structure with the same residue numbers, and it must survive later renumberings (chainsel's extraction, rpxdock's forced 1..N-per-chain dump, hbdesigner's offset); this tool does not remap, the consumer does. For an inherited file whose target is renumbered from 1, `target_res` is in that local numbering.

## Columns collected (`ifacegeom_` prefix)

`<leaf>_status` is `OK`, `warn: no interface residues`, `error: …` or `missing`. `<leaf>_path` is the per-residue TSV (side, label, chain, resnum, resname, selected, selected_by, min_heavy_dist, min_cb_dist): every residue near the other side, selected or not, the evidence behind the lists. NA (empty) rather than 0 marks a field that did not apply (the design failed or the geometry was undefined).

| Column | Meaning |
| --- | --- |
| `method` | `vector` or `heavy`: which criterion produced the lists |
| `binder_chains`, `target_chains` | the two sides actually used (`auto` resolved) |
| `binder_len`, `target_len` | amino acids per side |
| `binder_res` | **the epitope**: binder interface residues as `A:12,A:15,…` |
| `n_binder_res` | their count |
| `binder_res_seq` | their one-letter codes in the same order (lets a later graft be checked by residue identity) |
| `binder_his`, `n_binder_his` | histidines already in the binder's own epitope |
| `target_res`, `n_target_res` | target interface residues, same format |
| `target_his`, `n_target_his` | **histidines in the target epitope**: the residues whose contacts change when the pH drops |
| `binder_com` | mass-weighted centre of the whole binder, `x,y,z` |
| `binder_iface_com` | centre of the binder-side interface residues = the **epitope COM** |
| `iface_com` | centre of both sides' interface residues = the axis endpoint |
| `axis_len` | `\|iface_com − binder_com\|` |
| `nterm_res`, `cterm_res` | which residues the termini are (`A:1` / `A:115`), so every projection is checkable by hand |
| `nterm_proj`, `cterm_proj` | **signed** projection (Å) of the terminus CA on the `binder_com → iface_com` axis, from `binder_com`: **> 0 = interface side**, **< 0 = far side** |
| `nterm_iface_dist`, `cterm_iface_dist` | terminus CA to `iface_com` (Å); grows with binder size |
| `nterm_iface_min_dist`, `cterm_iface_min_dist` | terminus CA to the **nearest epitope CA** (Å): the number that answers "would a tag sit in the epitope"; means the same on every design; `0.0` = the terminus *is* an interface residue |
| `seconds` | wall time of the design's measurement |

## Reading the numbers

**A C-terminus requirement is two-sided, and neither side is filtered on here.** If a C-terminal tag or fusion partner is needed, the C-terminus must stay clear of the target epitope while remaining reachable from the partner's N-terminus: `cterm_proj` wants to be negative, but not so negative that a linker cannot reach. Both are **recorded**; the thresholds come from the distribution. `cterm_iface_min_dist` is the sharper of the two C-terminus numbers (the centroid distance grows with binder size).

`n_target_his` is a **record, not a gate**: a target epitope with histidines has contacts that change with pH, which a pH-switch strategy may want to avoid, but a pool is finite, so exclude nothing before seeing the distribution.

## Invariants a worker can check

Against numbers the tool did not compute:

- `binder_len` / `target_len` equal an independent residue count (and, for filenames that encode the binder length, that length).
- `--method heavy` reproduces a brute-force heavy-atom contact set at the same cutoff exactly (both sides); `vector` is checked against `heavy` knowing the `--cb-dist-cut` caveat.
- `binder_res_seq` is in register with the structure at known hot-spot residues; known interface histidines appear in `target_his`.
- On a synthetic two-chain case with a known axis, `cterm_proj` and `nterm_proj` take the expected signs.
- Absent chain, the same chain on both sides, and a missing file each give an `error:` status with every metric NA while the task still exits 0.
- Known hotspots of the design intent should be contacted by the binders: the tool's selection recovering the intended epitope is independent evidence it is right.

## What it is blind to

- **Affinity and energetics**: `pyrosetta` (`if_dG`, `if_dSASA`); **shape and area**: `cms`; **hotspot recall and clashes**: `epitope`.
- **Whether the pose is right**: it measures the interface of whatever pose it is given.
- **Renumbering downstream**: the consumer owns the mapping (see Traps).
