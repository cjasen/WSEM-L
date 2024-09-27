from Bio.PDB import *
import gemmi
import numpy as np

f_in = "henLyzHmod.pdb"  # Archivo PDB a analizar
print('EstructuraPDB =', f_in)
print("pH = 5.6")
N = 130

# Archivos de salida
file = open("1DPX.map", "w")
file2 = open("1DPX_elec.map", "w")
file3 = open("rCalpha.txt", "w")  # Archivo para distancias entre carbonos alfa

st = gemmi.read_structure(f_in)
st.remove_ligands_and_waters()

# Contar átomos y preparar para encontrar la distancia media entre carbonos alfa
alpha_carbons = []
alpha_carbons_residue_indices = []

for model in st:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.name == "CA":  # Carbono alfa
                    alpha_carbons.append(atom)
                    alpha_carbons_residue_indices.append(residue.seqid.num)

# Calcular la distancia media entre carbonos alfa
distances = []
for i in range(len(alpha_carbons) - 1):
    for j in range(i + 1, len(alpha_carbons)):
        d = np.sqrt((alpha_carbons[i].pos.x - alpha_carbons[j].pos.x) ** 2 +
                    (alpha_carbons[i].pos.y - alpha_carbons[j].pos.y) ** 2 +
                    (alpha_carbons[i].pos.z - alpha_carbons[j].pos.z) ** 2)
        distances.append(d)
        print(alpha_carbons_residue_indices[i], alpha_carbons_residue_indices[j], d, file=file3)

mean_distance = np.mean(distances)
print("Distancia media entre carbonos alfa =", mean_distance, "Å")

# Ahora usamos esta distancia media como r_cut
#r_cut = mean_distance
#print("r_cut =", r_cut, "Å")
r_cut=6

# Recuento de átomos
Nat = 0
for model in st:
    for chain in model:
        for residue in chain:
            for atom in residue:
                Nat += 1
Nat += 2
print(Nat)

nombre = [0 for _ in range(Nat)]
residuo = [0 for _ in range(Nat)]
carga = [0 for _ in range(Nat)]
x = [0 for _ in range(Nat)]
y = [0 for _ in range(Nat)]
z = [0 for _ in range(Nat)]

delta = np
