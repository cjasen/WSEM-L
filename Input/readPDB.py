from Bio.PDB import*
import gemmi
import numpy as np

f_in="henLyzHmod.pdb" #.pdb es necesario descargarlo del PDB
print('EstructuraPDB =',f_in)
print("pH = 5.6")
N = 130
r_cut = 6.0
print("N =",N,"r_cut =",r_cut)

file = open("1DPX.map", "w")
file2 = open("1DPX_elec.map","w")

j=0
i=0
Nat=0

st = gemmi.read_structure(f_in)
st.remove_ligands_and_waters()
print(st)

for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    Nat=Nat+1 
Nat = Nat +2

print(Nat)

nombre=[0 for i in range(Nat)]
residuo=[0 for i in range(Nat)]
carga = [0 for i in range(Nat)]
x=[0 for i in range(Nat)]
y=[0 for i in range(Nat)]
z=[0 for i in range(Nat)]

delta = np.zeros((N, N), dtype=int)
dist = np.zeros((N, N))
cp = np.zeros((N, N), dtype=float)
cnp = np.zeros((N, N), dtype=float)

#bucle en átomos. asignamos a cada átomo su nombre,residuo,posicióon
for model in st:
    for chain in model:
        for residue in chain:
            for atom in residue:
                match residue.name:
                    case 'ARG':
                        if(atom.name=='NE' or atom.name=='NH1' or atom.name=='NH2'):
                            carga[atom.serial] = 0.33
                    case 'LYS':
                        if(atom.name=='NZ'):
                            carga[atom.serial] = 1.0
                    case 'ASP':
                        if(atom.name=='OD1' or atom.name=='OD2'):
                            carga[atom.serial] = - 0.5
                    case 'GLU':
                        if(atom.name=='OE1' or atom.name=='OE2'):
                            carga[atom.serial]= - 0.5
                    case 'HIS':
                        if(atom.name=='ND1' or atom.name=='NE2'):
                            carga[atom.serial]= + 0.5
                nombre[atom.serial]=atom.name
                residuo[atom.serial]=int(residue.seqid.num)
                x[atom.serial]=atom.pos.x
                y[atom.serial]=atom.pos.y
                z[atom.serial]=atom.pos.z                 

#Bucle en átomos. Creamos el mapa elctrostático considerando todos los cotactos (i=j,i=j+1 incluidos)
for i in range(Nat-1):
    for j in range(i+1,Nat):
        d = np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2)
        if(d<r_cut):
            delta[residuo[i], residuo[j]] += 1
        if(carga[i]!=0 and carga[j]!=0): 
            print(int(residuo[i]),int(residuo[j]),carga[i],carga[j],d,file = file2)




#Bucle en residuos. Mapa VdW: contactos i>j+2. Mapa solvatación: todos los contactos 
for i in range(1,N-2):
    for j in range(i+2,N):
        if (delta[i][j]!=0):
            print(i,j,int(delta[i][j]),file = file)


