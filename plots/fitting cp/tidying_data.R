# Lee el archivo "profthermo.dat"
data <- read.table("profthermo_full.dat", header = FALSE)

# Define las constantes conocidas
a <- 0#1.575E-003  
b <- 0#5.541E-003  
M <- 14.400  

# Extrae las columnas necesarias
T <- data$V2  # Columna 2 es T
Cp_prime <- data$V8  # Columna 8 es Cp'

# Calcula Cp utilizando la fórmula proporcionada
Cp <- Cp_prime - (a + b * T) * M +5

# Reemplaza la columna Cp' con los nuevos valores de Cp
data$V8 <- Cp

# Guarda el nuevo archivo "profthermo_fixed.dat"
write.table(data, "profthermo_fixed.dat", row.names = FALSE, col.names = FALSE, sep = "\t")

# Haz un gráfico de Cp v. T
library(ggplot2)
ggplot(data, aes(x = T, y = Cp)) +
  geom_line(color = "blue") +
  labs(title = "Gráfico de Cp vs. T", x = "T (Temperatura)", y = "Cp") +
  theme_minimal()
