from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser()
structure = parser.get_structure("1pht_mutated", "1pht_mutated.pdb")
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
