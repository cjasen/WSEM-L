# Carga las librerías necesarias
library(ggplot2)
library(dplyr)

# Lee el archivo "profthermo.dat"
data <- read.table("profthermo.dat", header = FALSE)

# Extrae las columnas necesarias
T <- data$V2
m <- data$V3
sigma <- data$V4
Cp <- data$V8

# Crea un data frame para los gráficos
df <- data.frame(T, m, sigma, Cp)

# Primer gráfico: m v. T y sigma v. T
ggplot(df, aes(x = T)) +
  geom_line(aes(y = m, color = "m"), size = 1) +
  geom_line(aes(y = sigma, color = "sigma"), size = 1) +
  scale_color_manual(values = c("m" = "purple", "sigma" = "blue")) +
  labs(x = "T", y = "", color = "Legend") +
  theme_minimal()

# Segundo gráfico: Cp v. T y sigma v. T
ggplot(df, aes(x = T)) +
  geom_line(aes(y = Cp/max(Cp)+0.5, color = "Cp"), size = 1) +
  geom_line(aes(y = sigma, color = "sigma"), size = 1) +
  scale_color_manual(values = c("Cp" = "red", "sigma" = "blue")) +
  labs(x = "T",y="", color = "Legend") +
  theme_minimal()


