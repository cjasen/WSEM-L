# Leer el archivo 2PHT.map
archivo <- "2PHT.map"
datos <- read.table(archivo, header = FALSE, col.names = c("i", "j", "n"))

# Crear un vector para almacenar los contactos totales de cada residuo (1 a 83)
residuos <- 1:83
contactos_por_residuo <- numeric(length(residuos))  # Inicializar con ceros

# Sumar contactos para cada residuo considerando tanto i como j
for (residuo in residuos) {
  contactos_i <- sum(datos$n[datos$i == residuo], na.rm = TRUE)  # Contactos como i
  contactos_j <- sum(datos$n[datos$j == residuo], na.rm = TRUE)  # Contactos como j
  contactos_por_residuo[residuo] <- contactos_i + contactos_j
}

# Crear el data frame con las medias de contactos
contactos_df <- data.frame(Residuo = residuos,
                           MediaContactos = contactos_por_residuo)

# Extraer la información del residuo
contactos_residuo <- contactos_df[contactos_df$Residuo == 57, "MediaContactos"]

# Comparar la media de contactos de cada residuo con el residuo
comparacion <- contactos_df
comparacion$DiferenciaCon57 <- comparacion$MediaContactos - contactos_residuo

# Mostrar resultados
print(comparacion)

library(ggplot2)
ggplot(comparacion, aes(x = Residuo, y = MediaContactos, fill = Residuo == 57)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = contactos_residuo_57, color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("skyblue", "red"), guide = "none") +  # Asigna colores
  labs(title = "Media de contactos por residuo",
       x = "Residuo", y = "Media de contactos") +
  theme_minimal()