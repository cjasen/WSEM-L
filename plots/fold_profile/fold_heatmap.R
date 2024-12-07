library(ggplot2)
library(reshape2)
library(scales) # Para la notación científica en los ejes

# Crear la matriz diagonalizada para que la parte superior izquierda sea NA
matrix_data <- (as.matrix(read.table("fold_profile.txt", header = FALSE)))
matrix_data_upper_white <- matrix_data
matrix_data_upper_white[lower.tri(matrix_data_upper_white, diag = TRUE)] <- NA

# Convertir la matriz modificada a formato largo
matrix_long <- melt(matrix_data_upper_white)
colnames(matrix_long) <- c("Row", "Column", "Value") # Renombrar las columnas

# Convertir columnas a numéricas
matrix_long$Column <- as.numeric(matrix_long$Column)
matrix_long$Row <- as.numeric(matrix_long$Row)

# Crear el mapa de calor con escala logarítmica
heatmap_plot <- ggplot(matrix_long, aes(x = Column, y = Row, fill = Value)) +
  geom_raster(na.rm = FALSE) + # Transición suave sin bordes
  scale_fill_gradientn(
    colors = c("gray", "blue", "red"), # Gradiente de colores
    trans = "log10", # Transformación logarítmica base 10
    na.value = "white", # NA en blanco
    labels = scales::scientific # Etiquetas en notación científica
  ) +
  scale_x_continuous(
    breaks = seq(0, ncol(matrix_data), by = 10), # Mostrar cada 10 columnas
    labels = seq(0, ncol(matrix_data), by = 10) # Etiquetar con los números de las columnas
  ) +
  scale_y_continuous(
    breaks = seq(0, nrow(matrix_data), by = 10), # Mostrar cada 10 filas
    labels = seq(0, nrow(matrix_data), by = 10) # Etiquetar con los números de las filas
  ) +
  labs(
    x = "Residue b",
    y = "Residue a",
    fill = "Value (log10)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

# Mostrar el gráfico
print(heatmap_plot)

