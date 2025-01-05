# Cargar las librerías necesarias
library(ggplot2)
library(dplyr)
library(tidyr)

# Leer los datos del archivo
leer_datos <- function(ruta_archivo) {
  datos <- read.table(ruta_archivo, header = FALSE)
  colnames(datos) <- c("Temperatura", paste0("f", 1:(ncol(datos) - 1)))
  return(datos)
}

# Función para graficar los datos
graficar_todas_las_temperaturas <- function(datos) {
  # Transformar los datos para ggplot
  datos_largos <- datos %>%
    pivot_longer(-Temperatura, names_to = "Indice", values_to = "Valor") %>%
    mutate(Indice = as.numeric(gsub("f", "", Indice)))
  
  # Crear el gráfico
  ggplot(datos_largos, aes(x = Indice, y = Valor, color = factor(Temperatura))) +
    geom_line(size = 1.2) +      # Ajusta el grosor de las líneas
    scale_x_log10() +            # Escala logarítmica en el eje x
    labs(x = "L", y = "f(L)", color = "T (K)") +
    theme_minimal()
}

# Ejemplo de uso
ruta_archivo <- "f_L fold.txt"
datos <- leer_datos(ruta_archivo)
graficar_todas_las_temperaturas(datos)
