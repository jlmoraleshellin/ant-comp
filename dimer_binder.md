
## Dimer binder design campaign

We want to design a pH-responsive binder for EGFR that only binds EGFR at low pH. The goal is to design a binder that dimerizes at pH 7.4 via a dimerization interface that contains neutral histidines, similar as in https://www.science.org/doi/10.1126/science.aav7897. Upon lowering the pH, the dimerization interface should fall apart, exposing the actual binding interface for the EGFR receptor. Finally, we will link the two binder monomers with a flexible linker to increase the affinity for the intramolecular interface. This strategy allows us to do pH agnostic design, where we ignore protonated histidines during the design because we only need to design a stable interface containing histidine hydrogen bond networks at neutral pH. Additionally, our ideal binder binds to both human and mouse EGFR. 

Other hard requirements:
- The complete sequence (binder monomer - linker - binder monomer) can be no more than 250 amino acids. Let's reserve 20 aa for the linker, that leaves 115 aa per binder monomer.
- The binder interface should not be located too close to the C-terminus - wet lab validation will fuse the C-terminus to a GFP and splitstrep tag, so a binding interface too close to that terminus might be inhibited by the tag.

Here is my intended pipeline:

**Phase 1 covers initial binder design for EGFR.** Crucially, we want to obtain a good binding interface regardless of pH. Ideally the binder targets a hydrophobic patch: suggested hotspots (human EGFR numbering) are L325, P349, F412, V417, I467. We want to avoid binding sites on EGFR that contain histidines, because their protonation at low pH might ruin the interface. 

1. Backbone design: bindcraft 2 with multi-target optimization (human and mouse EGFR). Here we want to filter for binders for which the C-terminus is relatively close to the EGFR target. Ideally, calculate the binder center of mass, the binder-target interface center of mass, calculate the normal plane to that vector. The C-terminus should be located on the interface side of the plane.

2. Redesign of successful binders with Atomium and filtering using alphafold2.  

Alternatively, I want to be able to introduce good binders by other people in our team into my pipeline at this point.

**Phase 2 covers dimerisation and Hbond network design**
1. Input pdbs should have both target and binder. Use Pyrosetta to select residues that are at the binding interface (for example with InterGroupInterfaceByVector selecotr, default 8 angstrom distance). Extract the binder into a single pdb file and retain which binder residues are involved in the interface via an accompanying file. Clean up the pdb file (rename to chain A, renumber sequence if necessary) and we should retain a aa chain:res id mapping for the binder monomer pdb to the binder monomer in the target input. This mapping should be able to requested multiple design steps down the line.
2. Input the binder into RPXdock and run with C2 symmetry.
3. Filter for good docks. Good docks are based on rpxdocks own scoring. We want docked dimers where the dimerization interface is close to the binder-EGFR interface, but ideally does not overlap (completely). This can be achieved by calculating the center of masses of both interfaces and filtering on the distance. Ideally, the C-terminus of one monomer should be placed in reach of the N-terminus of the other binder. A simple distance calculation based on terminal Calphas suffices. These two metrics should be enough to determine correct positioning.
4. Superpose dimerized binder onto monomer binder-EGFR structure (using the binder in the binder-EGFR complex, and one of the monomers of the dimer, but apply the rigid body transformation to the entire dimer) and calculate a clash score for the second binder copy and the EGFR. Both superposition and clash score can be done with pyrosetta (use fa_rep)-> continue with binder dimers that show strong clashes between the second binder copy and EGFR (because we want the dimer to prevent binding). Export a new output file with the binders in the position after the rigid body transformation (without the target) - this will simplify the next step.
5. Identify the binder residues that are part of the binder-EGFR interface (from the pdb-info labels), and transplant them from the binder-EGFR structure to the dimer structure (for both dimers). Do the correct rigid body transformation to place them on the second binder copy. I want the exact same rotamers as in the input file.
5. Use HBdesigner to design Hbond network for a binder dimer  at dimerization interface. Keep binder residues involved with the binder-EGFR interface (which we placed in the previous step) fixed. Prefer HisN --- H donor bonds paired with HisH - carboxylate hbonds for pKa adjustment. HBdesigner can handle symmetric interfaces. Check flaggs for what might help!
6. Filter for designs with histidines in the hbond network
7. Fill in the remainder of the sequence with atomium (keeping the binder-egfr interface residues intact, as well as the hbond network residues)
8. Validate using structure prediction with boltz. Predict the dimer structure, the egfr-monomer structure, and the egfr-dimer structure.
