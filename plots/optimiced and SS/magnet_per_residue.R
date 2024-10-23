# Cargar los datos desde el archivo "magnet.dat"
data <- read.table("magnet.dat", header = FALSE, col.names = c("Temperatura", "Residuo", "m", "s", "m_s"))

# Crear el dataframe
df <- as.data.frame(data)

# Calcular los promedios de "m" y "s" por cada temperatura
df_promedio <- aggregate(cbind(m, s) ~ Temperatura, data = df, FUN = mean)

# Escribir fichero
write.table(df_promedio, file = "ms_average_of_residues.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
