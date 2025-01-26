import pandas as pd

# Parámetros dados por el usuario    
a = 1.575E-003  # Valor de a
b = 5.541E-003  # Valor de b 
M = 9.630  # Masa molar (g/mol)

# Leer el archivo de texto original
input_file = "excel_data.txt"
data = pd.read_csv(input_file, delim_whitespace=True, header=None, names=["T_C", "Cp_kcal"])

# Convertir las unidades y calcular Cp'
data["T_K"] = data["T_C"] + 273.15  # De °C a K
data["Cp_kj"] = data["Cp_kcal"] * 4.184 / 1000  # De cal/mol a kJ/mol
data["Cp_prime"] = data["Cp_kj"] +(a + b * data["T_K"]) * M +2.4 # Cp' = Cp + (a + b*T) * M

# Filtrar uno de cada 4 datos
filtered_data = data.iloc[::5]  # Selecciona cada 4 filas

# Seleccionar las columnas necesarias y guardar el nuevo archivo
output_file = "2PHT_expCp.txt"
filtered_data[["T_K", "Cp_prime"]].to_csv(output_file, sep="\t", index=False, header=False)

print(f"Datos convertidos guardados en: {output_file}")
