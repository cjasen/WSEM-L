# Cargar las librerías necesarias
library(dplyr)
library(ggplot2)

# Función para procesar un archivo y normalizar los datos
process_file <- function(filename) {
  # Leer el archivo
  data <- read.table(filename, header = FALSE, col.names = c("T", "i", "m_i", "s_i", "diff"))
  
  # Definir los intervalos de residuos en estado hélice
  helices <- c(30:40, 46:48, 50:53, 71:73)
  
  # Seleccionar la columna a usar según el archivo
  col_to_use <- ifelse(grepl("magnet2", filename), "diff", "diff")
  
  # Filtrar los residuos en hélice y calcular el promedio por temperatura
  result <- data[data$i %in% helices, ] %>%
    group_by(T) %>%
    summarise(avg_value = mean(.data[[col_to_use]])) %>%
    ungroup()
  
  # Normalizar los valores de avg_value de manera independiente para este archivo
  min_value <- min(result$avg_value)
  max_value <- max(result$avg_value)
  result$avg_value <- (result$avg_value - min_value) / (max_value - min_value)
  
  # Agregar una columna para identificar el archivo con nombres personalizados
  result$file <- case_when(
    grepl("magnet1", filename) ~ "Wild type",
    grepl("magnet2", filename) ~ "Mutant",
    grepl("magnet3", filename) ~ "TCEP"
  )
  return(result)
}

# Función para leer los datos experimentales
process_exp_file <- function(filename) {
  # Leer el archivo experimental
  data <- read.table(filename, header = FALSE, col.names = c("T", "m_i"))
  
  # Normalizar los valores de m_i de manera independiente para este archivo
  min_value <- min(data$m_i)
  max_value <- max(data$m_i)
  data$m_i <- (data$m_i - min_value) / (max_value - min_value)
  
  # Agregar una columna para identificar el archivo con nombres personalizados
  data$file <- case_when(
    grepl("exp_magnet1", filename) ~ "Wild type",
    grepl("exp_magnet2", filename) ~ "Mutant",
    grepl("exp_magnet3", filename) ~ "TCEP"
  )
  return(data)
}

# Listar los archivos simulados y experimentales
sim_files <- c("magnet1.dat", "magnet2.dat", "magnet3.dat")
exp_files <- c("exp_magnet1.txt", "exp_magnet2.txt", "exp_magnet3.txt")

# Procesar los archivos simulados
sim_results <- lapply(sim_files, process_file) %>% bind_rows()

# Procesar los archivos experimentales
exp_results <- lapply(exp_files, process_exp_file) %>% bind_rows()

# Graficar con ggplot2
ggplot() +
  # Agregar las líneas simuladas
  geom_line(data = sim_results, aes(x = T, y = avg_value, color = file, group = file), size = 1) +
  # Agregar los puntos experimentales con relleno de color y diferentes formas
  geom_point(data = exp_results, aes(x = T, y = m_i, fill = file, shape = file), size = 4, color = "black") +
  labs(
    x = "T (K)",
    y = expression("<m(1-" * sigma * ")> for helixes"),
    color = "Protein",
    fill = "Protein",
    shape = "Protein"
  ) +
  theme_minimal(base_size = 18) + # Aumentar el tamaño de todo el texto
  scale_color_manual(
    values = c("Wild type" = "blue", "Mutant" = "red", "TCEP" = "green")
  ) +
  scale_fill_manual(
    values = c("Wild type" = "blue", "Mutant" = "red", "TCEP" = "green")
  ) +
  scale_shape_manual(
    values = c("Wild type" = 22, "Mutant" = 21, "TCEP" = 24) # Asignación de formas: cuadrado, círculo, triángulo
  )