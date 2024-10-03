import sys
import numpy as np
from Bio import PDB

def get_disulfide_bonds(pdb_code):
    # Download the PDB file from the PDB database
    pdb_list = PDB.PDBList()
    pdb_file = pdb_list.retrieve_pdb_file(pdb_code, pdir='.', file_format='pdb')
    
    # Parse the PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('Structure', pdb_file)
    
    # List to store residues that form disulfide bonds
    disulfide_bonds = []

    # Iterate over all models in the structure
    for model in structure:
        # Get all residues
        for chain in model:
            cysteines = []
            for residue in chain:
                # We are only interested in cysteines (CYS)
                if residue.get_resname() == 'CYS':
                    # Check if the residue contains the sulfur atom (SG)
                    if 'SG' in residue:
                        cysteines.append(residue)

            # Check for possible disulfide bond pairs
            for i in range(len(cysteines)):
                for j in range(i + 1, len(cysteines)):
                    # Calculate the distance between SG atoms of cysteines
                    distance = cysteines[i]['SG'] - cysteines[j]['SG']
                    # If the distance is less than 2.5 Å, consider it a disulfide bond
                    if distance < 2.5:
                        disulfide_bonds.append((cysteines[i].get_id()[1], cysteines[j].get_id()[1]))

    # Convert the list of disulfide bonds into a numpy matrix
    bond_matrix = np.array(disulfide_bonds)

    return bond_matrix

def save_matrix_to_file(matrix, filename):
    # Save the matrix to a CSV file
    np.savetxt(filename, matrix, delimiter=',', fmt='%d')

if __name__ == "__main__":
    # Get the PDB code from the command line argument
    pdb_code = sys.argv[1]
    
    # Get the disulfide bond matrix
    matrix = get_disulfide_bonds(pdb_code)
    
    # Save the matrix to a file
    save_matrix_to_file(matrix, 'disulfide_bonds.csv')


