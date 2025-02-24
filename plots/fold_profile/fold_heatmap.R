library(ggplot2)
library(reshape2)
library(scales) # Para la notación científica en los ejes

# Leer todas las matrices desde el archivo y separarlas en una lista
read_matrices <- function(file_path, matrix_size) {
  data <- as.matrix(read.table(file_path, header = FALSE))
  num_matrices <- nrow(data) / matrix_size
  matrices <- lapply(1:num_matrices, function(i) {
    data[((i - 1) * matrix_size + 1):(i * matrix_size), ]
  })
  return(matrices)
}

# Crear heatmap para una matriz con escala global fija
create_heatmap <- function(matrix_data, matrix_index, global_min, global_max) {
  # Crear la matriz diagonalizada para que la parte superior izquierda sea NA
  matrix_data_upper_white <- matrix_data
  matrix_data_upper_white[lower.tri(matrix_data_upper_white, diag = TRUE)] <- NA
  
  # Convertir la matriz modificada a formato largo
  matrix_long <- melt(matrix_data_upper_white)
  colnames(matrix_long) <- c("Row", "Column", "Value") # Renombrar las columnas
  
  # Asegurarse de que las columnas "Row" y "Column" sean numéricas
  matrix_long$Row <- as.numeric(matrix_long$Row)
  matrix_long$Column <- as.numeric(matrix_long$Column)
  
  # Crear el mapa de calor con escala logarítmica
  heatmap_plot <- ggplot(matrix_long, aes(x = Column, y = Row, fill = Value)) +
    geom_raster(na.rm = FALSE) + # Transición suave sin bordes
    scale_fill_gradientn(
      colors = c("gray", "blue", "red"), # Gradiente de colores
      limits = c(global_min, global_max), # Usar límites globales
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
      #fill = "Value (log10)",
      title = paste("T=", sprintf("%02d", matrix_index),"K")
    ) +
    theme_minimal(base_size=22) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10),
      #legend.position = "none" # Eliminar la leyenda
    )
  
  return(heatmap_plot)
}

# Parámetros
file_path <- "sigma_profile cont.txt"
matrix_size <- 83 # Dimensión de las matrices

# Leer las matrices del archivo
matrices <- read_matrices(file_path, matrix_size)

# Calcular el rango global de los valores de todas las matrices
global_min <- 1E-16#min(sapply(matrices, function(x) min(x, na.rm = TRUE)))
global_max <- 1#max(sapply(matrices, function(x) max(x, na.rm = TRUE)))

# Parámetros de temperatura inicial y delta
T_inicial <- 330  
delta <- 40    

# Crear heatmaps para todas las matrices con la misma escala de colores
for (i in seq(length(matrices))) {
  # Calcular la temperatura correspondiente
  temperatura <- T_inicial + (i - 1) * delta
  
  # Crear el heatmap
  heatmap <- create_heatmap(matrices[[i]], temperatura, global_min, global_max)
  
  # Guardar cada heatmap como archivo con nombre basado en la temperatura
  ggsave(
    filename = sprintf("Heatmaps/ss_sigma_cont_%d.jpg", temperatura),
    plot = heatmap,
    width = 8,
    height = 6
  )
}
