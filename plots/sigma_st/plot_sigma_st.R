# Cargar librería
library(ggplot2)

# Leer el archivo
data <- read.table("sigma_st.dat", header = FALSE)

# Crear un dataframe con las columnas necesarias
df <- data.frame(
  T = data[[2]], # Temperatura
  loop = data[[3]], # Datos para (10,16)
  a_helix = data[[4]] # Datos para (119,126)
)

# Crear el gráfico con ggplot
ggplot(df, aes(x = T)) +
  geom_line(aes(y = loop, color = "loop (10,16)")) +
  geom_line(aes(y = a_helix, color = "a-helix (119,126)")) +
  labs(
    x = "Temperature (K)",
    y = "<prod m*sigma>",
    color = "Structure",
  ) +
  theme_minimal()
