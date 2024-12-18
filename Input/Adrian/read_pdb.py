from Bio.PDB import *
import gemmi
import numpy as np

f_in = "1pht.pdb"  # .pdb es necesario descargarlo del PDB
N = 87  # if .map is an empty file, add two more residues than it has

pdb_code = f_in.split(".")[0].upper()
print('EstructuraPDB =', f_in)
print("pH = 5.6")
r_cut = 6.0
print("N =", N, "r_cut =", r_cut)
file1_name = f"{pdb_code}.map"
file2_name = f"{pdb_code}_elec.map"
file3_name = f"{pdb_code}_rCalpha.txt"

# Abrir los archivos para escritura
file = open(file1_name, "w")
file2 = open(file2_name, "w")
#  not all .pdb start at residue 1, this can cause conflict. I have added some -2 to the calculations, so check that !!!!!!!!!!!!!!!
file3 = open(file3_name, "w") 

j = 0
i = 0
Nat = 0

st = gemmi.read_structure(f_in)
st.remove_ligands_and_waters()
print(st)

for model in st:
    for chain in model:
        for residue in chain:
            for atom in residue:
                Nat = Nat + 1
Nat = Nat + 2

print(Nat)

nombre = [0 for i in range(Nat)]
residuo = [0 for i in range(Nat)]
carga = [0 for i in range(Nat)]
x = [0 for i in range(Nat)]
y = [0 for i in range(Nat)]
z = [0 for i in range(Nat)]

delta = np.zeros((N, N), dtype=int)
dist = np.zeros((N, N))
cp = np.zeros((N, N), dtype=float)
cnp = np.zeros((N, N), dtype=float)

# Bucle en átomos. Asignamos a cada átomo su nombre, residuo y posición
for model in st:
    for chain in model:
        for residue in chain:
            for atom in residue:
                match residue.name:
                    case 'ARG':
                        if(atom.name == 'NE' or atom.name == 'NH1' or atom.name == 'NH2'):
                            carga[atom.serial] = 0.33
                    case 'LYS':
                        if(atom.name == 'NZ'):
                            carga[atom.serial] = 1.0
                    case 'ASP':
                        if(atom.name == 'OD1' or atom.name == 'OD2'):
                            carga[atom.serial] = - 0.5
                    case 'GLU':
                        if(atom.name == 'OE1' or atom.name == 'OE2'):
                            carga[atom.serial] = - 0.5
                    case 'HIS':
                        if(atom.name == 'ND1' or atom.name == 'NE2'):
                            carga[atom.serial] = + 0.5
                nombre[atom.serial] = atom.name
                residuo[atom.serial] = int(residue.seqid.num)
                x[atom.serial] = atom.pos.x
                y[atom.serial] = atom.pos.y
                z[atom.serial] = atom.pos.z

# Crear lista para almacenar posiciones de carbonos alpha
calpha_positions = []
calpha_residues = []

# Identificar posiciones de carbonos alpha
for model in st:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.name == "CA":  # Carbono alpha
                    calpha_positions.append((atom.pos.x, atom.pos.y, atom.pos.z))
                    calpha_residues.append(residue.seqid.num)

# Calcular distancias entre carbonos alpha y escribir en file3
for i in range(len(calpha_positions) - 1):
    for j in range(i + 1, len(calpha_positions)):
        d = np.sqrt(
            (calpha_positions[i][0] - calpha_positions[j][0]) ** 2 +
            (calpha_positions[i][1] - calpha_positions[j][1]) ** 2 +
            (calpha_positions[i][2] - calpha_positions[j][2]) ** 2
        )
        print(calpha_residues[i]-2, calpha_residues[j]-2, d, file=file3)

# Bucle en átomos. Creamos el mapa electrostático considerando todos los contactos (i=j, i=j+1 incluidos)
for i in range(Nat-1):
    for j in range(i+1, Nat):
        d = np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2)
        if(d < r_cut):
            delta[residuo[i], residuo[j]] += 1
        if(carga[i] != 0 and carga[j] != 0):
            print(int(residuo[i]-2), int(residuo[j]-2),
                  carga[i], carga[j], d, file=file2)

# Bucle en residuos. Mapa VdW: contactos i>j+2. Mapa solvatación: todos los contactos
for i in range(1, N-2):
    for j in range(i+2, N):
        if (delta[i][j] != 0):
            print(i-2, j-2, int(delta[i][j]), file=file)

# Cerrar los archivos
file.close()
file2.close()
file3.close()
