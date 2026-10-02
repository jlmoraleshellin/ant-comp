---
name: all-tools
description: The catalog of every tool available in this workspace — what question each one answers, its action (create/update), the column it consumes, the columns it returns, and which skill to load before running it. Also what is NOT available, and the input-column wiring between steps. Load at the start of a campaign, before planning a chain of steps, and before concluding that no tool covers a measurement.
---

# What is available

This page is an **index, not a manual**. It tells you which tool answers which question and what it would put in the table. **Load that tool's own skill before composing the run** — flags, traps and column meanings live there, and only there.

Authoritative check, when this page and reality disagree: ask the orchestrator for `sapia run --help` in the workstation. It lists every registered tool, built-in and custom.

Two things decide where a tool's output lands, and you never name the output table:

- **`create`** mints a **child table** (`gen+1`), one row per new entity, linked to its parent.
- **`update`** annotates the table it reads, **in place**, adding columns. `-t` is required.

Every tool leaf-prefixes its columns and writes `<leaf>_status`; `OK` is the only success signal. **NA is "not applicable", never zero** — do not read a blank as a pass.

## The catalog

### Make entities — `create`, mints a child table

| Tool | Answers | Consumes | Returns (prefix `<tool>_`) | Skill |
| --- | --- | --- | --- | --- |
| **rfdiffusion3** | "Give me backbones." De novo, or conditioned on an input PDB (motif, binder against a target). | `pdb_path`; no `-t` = root run into `table0` | `path`, `iteration`, `rfd3_batch`, `rfd3_model`, `rfd3_ca_rmsd_to_input`, per-chain length | `rfdiffusion3` |
| **proteinmpnn** | "What sequence folds this backbone?" | `rfdiffusion_path` ← **wrong for rfd3, pass `-i rfdiffusion3_path`** | `sequence`, `score` (lower better), `seq_recovery`; rows `<parent>_f1…` | `proteinmpnn` |
| **atomium** | Same question, private noise-conditioned model. Use to diversify against MPNN. | `rfdiffusion3_path` | `sequence`, `sample`, `temperature`, `seq_rec`; rows `<parent>_a1…`. **No score column — you cannot rank the way MPNN allows.** | `atomium` |
| **hbdesigner** | "Put a buried hydrogen-bond network into this backbone." Designs 2–6 polar positions (monomer, or across an interface) and **leaves every other position as glycine** — a network stub, not a sequence. Feeds proteinmpnn with those positions fixed. | `rfdiffusion3_path` | `path`, `rank`, `hb_score_full`, `hb_score_hb`, `avg_burial`, `saturation`, `buried_heavy_unsats`, `buried_unsat_hpol`, `network`, `network_seq`, `n_network_res`, `network_chains`, `fixed_positions`, `fix1…fixK`, `resnum_offset`, `resnum_shift_max`, `grafted`, `graft_identity`, `n_res_total`; rows `<parent>_hb1…`. **Fan-out is ≤ `--top-k`, and a `_hb0` row is a failure record, not a design.** | `hbdesigner` |
| **rpxdock** | "Is there a rigid-body placement of this scaffold that makes a good symmetric interface?" RPXdock docks it into a **one-component** architecture and scores each pose with precomputed residue-pair motif tables. CPU-only. **SCOPE: CYCLIC ONLY (`C2`–`C17`, `CxSTACK`) with a single-chain monomer** — the dihedrals `Dx_y` and the one-component cages are deliberately **out of scope and unverified**; read the scope section of the `rpxdock` skill before touching them. **Pass `--hscore-files afilmv_ehl` on every run in this campaign** — the default `ilv_h` is helix-only. | no default — name the structure column | `path` (**the symmetric assembly**), `score`, `rpx`, `ncontact`, `model`, `rank`, `reslb`/`resub`, `disp`, `architecture`, `nfold`, `hscore`, `scaffold_path`, `result_path`, `n_docks`, `n_res`, `n_chains_in`, `frac_helix`/`frac_sheet`/`frac_loop`, `input_com_dist`, `recentered`; rows `<parent>_d1…`. **Multi-component (T32, I32, AXLE_*, PLUG_*, ASYM, layers) is refused at submit time; so is a multi-chain input — >1 chain for `Cx`, >2 for `Dx_y`/cages. The input must be origin-centred, and a pre-aligned cyclic ASU for cages.** | `rpxdock` |
| **bindcraft2** | "Give me binders against this target" — the *whole* campaign in one step: AF2 hallucination + MPNN + refold + filter, looping until enough are accepted. **`--trajectory-only` stops it at backbones**, to redesign with `atomium`/`proteinmpnn` instead. | `pdb_path` (name the real column); no `-t` = root run from `--target-pdb` / `--shipped-target` | `sequence`, `i_pDAE`, `i_pTM`, `i_pAE`, `pLDDT`, `Interface_Residues`, `Interface_BuriedArea`, `Hotspot_Contact_Fraction`, `Binder_Length`, `rank`, `outcome`, `failed_filters`, `path` (**the complex**). **One task = one campaign; the row count is unknown until collect.** | `bindcraft2` |

### Prepare an input — `update`, cheap, no GPU

| Tool | Answers | Consumes | Returns | Skill |
| --- | --- | --- | --- | --- |
| **mkcomplex** | "Put the target chains back around this binder sequence, so the predictor folds the complex." Pure string work. | `proteinmpnn_sequence` | `sequence` (**the predictor's input column**), `n_chains`, `total_len`, `chain_lens`, `design_chain_index` | `mkcomplex` |
| **chainsel** | "Pull the binder chain out of this complex." Also merges split protomers back into one chain. | no default — you name the column | `path`, `chains`, `n_chains`, `n_res`, `n_atoms` | `chainsel` |

### Predict a structure — `update`, GPU, the expensive step

| Tool | Answers | Consumes | Returns | Skill |
| --- | --- | --- | --- | --- |
| **boltz** | "What does this sequence fold to?" The workhorse. Supports forced templates and per-chain MSA policy. | `proteinmpnn_sequence` (pass `-i mkcomplex_sequence` on a binder table, `-i atomium_sequence` after atomium) | `path`, `confidence_score`, `ptm`, `iptm`, `protein_iptm`, `ligand_iptm`, `complex_plddt`, `complex_iplddt`, `complex_pde`, `complex_ipde` | `boltz` |
| **alphafold3** | Same question, second opinion. | `proteinmpnn_sequence` | `path`, `ranking_score`, `ptm`, `iptm`, `fraction_disordered`, `has_clash` | — none; read the source |
| **colabfold** | Same question, cheaper/older. | `proteinmpnn_sequence` | `path`, `avg_plddt`, `ptm`, `iptm`, `max_pae` | — none; read the source |

### Measure a design — `update`, the gates

| Tool | Answers | Consumes | Returns | Skill |
| --- | --- | --- | --- | --- |
| **usalign** | "Did it fold back to its own backbone?" and "did the target land?" Takes **`--col-a`/`--col-b`**, not `-i`. | two structure columns | `TM1`, `TM2`, `RMSD`, `ID1/ID2/IDali`, `L1/L2/Lali`, `sup_path` | `usalign` |
| **pyrosetta** | "Is it well packed, and what does the interface cost?" ref2015 after an optional FastRelax. | `boltz_path` | `total_score`, `score_per_res`, `score_raw`, `relax_ca_rmsd`, every weighted term, `sasa`, `sasa_hydrophobic`, `packstat`, `buried_unsat`, `dssp`, `if_dG`/`if_dSASA`/`if_hbonds`/`if_delta_unsat` | `pyrosetta` |
| **cms** | "How much real interface is there, and does it fit?" Contact molecular surface + shape complementarity, GPU. | no default — name the structure column | `target`, `binder`, `sc`, `sc_area`, `sc_median_dist`, `n_atoms_binder/_target`, `path` (**per-residue CMS: the epitope map**) | `cms` |
| **ringfit** | "Does the binder straddle two protomers, and would it clash with the rest of the assembly or its lipid belt?" (Very specific, you will almost never need it) | `rfdiffusion3_path` + `--ref-structure` | `align_rmsd`, `seq_match_frac`, `resnum_offset`, `bsa_t1/t2/total`, `bridge_ratio`, `hotspot_recall`, `n_clash`, `min_dist_ring`, `lipid_clash`, `path` | `ringfit` |
| **ifacegeom** | "**Which residues** are the epitope, and where are the termini relative to it?" The interface residue lists on **both** sides, as columns that ride lineage — plus the target-side histidine count and a signed C-/N-terminus projection. CPU-only, milliseconds per design. | no default — name the **complex** structure column | `binder_res` (**the epitope, `A:12,A:15,…`**), `n_binder_res`, `binder_res_seq`, `target_res`, `target_his`/`n_target_his`, `binder_com`, `binder_iface_com`, `iface_com`, `axis_len`, `nterm_res`/`cterm_res`, `cterm_proj`/`nterm_proj` (**signed; >0 = interface side**), `cterm_iface_dist`/`_min_dist`, `binder_len`, `path` (**per-residue evidence**) | `ifacegeom` |

## Question → tool

| You want to know | Ask |
| --- | --- |
| Give me backbones | `rfdiffusion3` |
| Give me binders against this target | `bindcraft2` (the whole campaign in one step), or `rfdiffusion3` binder mode if you want to compose the chain yourself |
| What sequence folds this | `proteinmpnn`, or `atomium` for diversity |
| Dock this scaffold into a C2/C3 and score the new interface | `rpxdock` (`score`, `rpx`, `ncontact`) — **CYCLIC ONLY, single-chain monomer**; cages/dihedrals are out of scope and unverified. Read its premise before trusting a score |
| Put a buried H-bond network in the core (or across an interface) first | `hbdesigner`, then `proteinmpnn` with `--fixed-positions` from its `fix<i>` columns |
| What does this sequence fold to | `boltz` (`alphafold3` / `colabfold` for a second opinion) |
| Did it fold back to its designed backbone | `boltz` → `usalign` (`--col-a` prediction, `--col-b` parent backbone) |
| Did the binder *stay put*, not just fold | `chainsel` the binder out, then `usalign` **in the target frame** — a separate question from fold, see `binder-campaign` |
| Did the target land under a forced template | `usalign` predicted target chains vs. the template. **This gate comes first.** |
| How big and how good is the interface | `cms` (`target`, `sc`), `pyrosetta` (`if_dG`, `if_dSASA`) |
| Which epitope residues does it actually cover | `ifacegeom` (`binder_res` / `target_res` — **as columns**, so they ride lineage into later tables), or `cms` → `cms_path` for a per-residue area table you cannot filter on |
| Does the target epitope contain histidines (pH sensitivity) | `ifacegeom` (`target_his`, `n_target_his`) |
| Where is the C-terminus relative to the binding face | `ifacegeom` (`cterm_proj` — signed, >0 = interface side; `cterm_iface_min_dist`) |
| Is it well packed / any buried unsats | `pyrosetta` (`packstat`, `buried_unsat`, `score_per_res`) |
| Does it bridge two protomers of an oligomer | `ringfit` (`bridge_ratio`, `hotspot_recall`) |
| Would it clash with the rest of the assembly | `ringfit` (`n_clash`, `min_dist_ring`, `lipid_clash`) |
| Fold the complex, not the binder alone | `mkcomplex` first, then `boltz -i mkcomplex_sequence` |

## Input-column wiring

The commonest silent failure in this workspace is a tool reading the wrong column and submitting **nothing**, or scoring the wrong thing. Defaults were written for a rfdiffusion → proteinmpnn → boltz chain; anything else needs `-i`.

| Coming from | Going to | Pass |
| --- | --- | --- |
| rfdiffusion3 | proteinmpnn / atomium | `-i rfdiffusion3_path` — the default is `rfdiffusion_path` (no 3) and matches nothing |
| rfdiffusion3 | hbdesigner | nothing — its default IS `rfdiffusion3_path`, and a missing column raises instead of submitting nothing |
| hbdesigner | proteinmpnn | `-i hbdesigner_path`, plus `--chains-to-design <hbdesigner_network_chains>` and `--fixed-positions '{hbdesigner_fix1},{hbdesigner_fix2}/'`. **Never fold `hbdesigner_path` directly — it is poly-glycine outside the network.** |
| atomium | boltz / af3 | `-i atomium_sequence` |
| rfdiffusion3 / boltz / chainsel | rpxdock | `-i <that>_path` — no default at all. The scaffold must be **origin-centred and single-chain** (a multi-chain input is now refused at submit time: >1 chain for `Cx`, >2 for `Dx_y`/cages — `chainsel` it first); check `rpxdock_n_chains_in` and `rpxdock_input_com_dist` after collecting. Columns are prefixed by the **dir-label**, so a `-l verification` run writes `rpxdock_verification_*` |
| rpxdock | proteinmpnn / atomium / cms / pyrosetta | `-i rpxdock_path` — the symmetric assembly. It is **backbone-only** (N/CA/C/O/CB + a CEN pseudo-atom) unless the run passed `--use-orig-coords` |
| bindcraft2 | boltz / af3 | `-i bindcraft2_sequence` — an **independent** check; bindcraft2's own scores come from the AF2 that designed the binder |
| bindcraft2 `--trajectory-only` | atomium / proteinmpnn | `-i bindcraft2_traj_path` (use `-l traj`, or the leaf collides with a full campaign's). The backbone is the **complex** — pass `--chains-to-design <binder chain>`, and filter on `bindcraft2_traj_completed` |
| bindcraft2 | chainsel / cms / usalign | `-i bindcraft2_path` — it is the **complex**, so `chainsel` the binder out first |
| proteinmpnn, binder campaign | boltz / af3 | `mkcomplex` first, then `-i mkcomplex_sequence` |
| boltz, binder | usalign / cms | `chainsel` first, then the `chainsel_path` it wrote |
| anything | cms / chainsel / ifacegeom | no default at all — you must name the column. `ifacegeom` needs the **complex** (binder + target in one file); a monomer has no interface |
| ifacegeom | chainsel / rpxdock / a later graft | nothing — its columns ride lineage. But `ifacegeom_binder_res` is in **this file's numbering**, and `chainsel`/`rpxdock` both renumber; the consumer does the remapping, not `ifacegeom` |
| anything | usalign | `--col-a` / `--col-b`, **not** `-i` |

Two structure columns exist for every predicted design: the **backbone** it was designed as (`rfdiffusion3_path`, in the parent table) and the **prediction** (`boltz_path`, in this one). Lineage resolves the parent for you; say which you mean.

## Not available for modal

- **`rfdiffusion` (v1), `openfold3`, `make_symmdef`, `align_symm_axis`** are registered but have **no `modal_image.py`**. If runnning on modal, edit them.
- **No tool installs weights.** rfd3's checkpoint Volume was populated by hand. A new tool needing weights is a setup task for the user, not something an agent can do end to end.


## Related skills

| Skill | Load when |
| --- | --- |
| `prosapia` | Composing any `sapia` command — run_dirs, tables, labels, the ready set |
| `binder-campaign` | Starting a binder campaign, and before reading any interface number |
| `authoring-a-tool` / `editing-a-tool` | Commissioning or bending a tool (usually via `tool-creator`) |

## One operational trap

- Custom tools (`atomium`, `bindcraft2`, `chainsel`, `cms`, `hbdesigner`, `mkcomplex`, `ringfit`, `rpxdock`) live in this project's `tools/` and are baked into the workstation image from the **local working directory**. If `sapia modal-shell` is launched from somewhere other than the project root, they are simply absent from `sapia run --help` — the built-ins still work, so it looks like the custom tool was never written rather than like a path problem.
- **No tool names its output table**, and none of them will reorder your campaign. You compose.
