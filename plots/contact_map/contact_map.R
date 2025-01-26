# Cargar librerías necesarias
library(ggplot2)
library(dplyr)

# Leer el archivo
data <- read.table("2PHT.map", header = FALSE, col.names = c("i", "j", "n"))

# Calcular la longitud del intervalo (j - i)
data <- data %>%
  mutate(interval_length = j - i)

# Agrupar y resumir los datos para contar la frecuencia de n para cada intervalo
heatmap_data <- data %>%
  group_by(interval_length, n) %>%
  summarise(count = n(), .groups = "drop")

# Crear el heatmap
ggplot(heatmap_data, aes(x = interval_length, y = n, fill = count)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "red", "blue"),
    values = scales::rescale(c(0, 1, max(heatmap_data$count))),
    name = "Frecuencia",
    limits = c(0, max(heatmap_data$count))
  ) +
  labs(
    title = "Heatmap de contactos (n) frente a la longitud del intervalo (j-i)",
    x = "Longitud del intervalo (j-i)",
    y = "Número de contactos (n)"
  ) +
  theme_minimal()

############

# Crear un gráfico de densidad 2D
ggplot(data, aes(x = interval_length, y = n)) +
  stat_density_2d(
    aes(fill = ..density..), 
    geom = "raster", 
    contour = FALSE
  ) +
  scale_fill_gradientn(
    colors = c("white", "lightblue", "blue", "red"),
    name = "Densidad"
  ) +
  labs(
    title = "Mapa de densidad 2D de contactos (n) frente a longitud del intervalo (j-i)",
    x = "Longitud del intervalo (j-i)",
    y = "Número de contactos (n)"
  ) +
  theme_minimal()
