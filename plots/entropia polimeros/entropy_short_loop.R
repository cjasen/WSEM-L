library(ggplot2)

# Parámetros
lp <- 6  # Longitud de persistencia en Amstrongs

# Define la función para el primer cálculo de S
calc_S1 <- function(r, delta_ji) {
  1.5 * log(delta_ji + 1) + 1.5 * (r^2 - 3.8^2) / ((delta_ji + 2) * 2 * 20 * 3.8)
}

# Define la función para el segundo cálculo de S y w
calc_S2 <- function(d, delta_ji) {
  lc <- delta_ji * 3.8
  w <- (5.0 * lp / (4.0 * lc)) - (2.0 * d^2 / lc^2) +
    (33.0 * d^4 / (80.0 * lp * lc^3)) + 
    (79.0 * lp^2 / (160.0 * lc^2)) + 
    (329.0 * d^2 * lp / (120.0 * lc^3)) - 
    (6799.0 * d^4 / (1600.0 * lc^4)) + 
    (3441.0 * d^6 / (2800.0 * lp * lc^5)) - 
    (1089.0 * d^8 / (12800.0 * lp^2 * lc^6))
  S <- -1.5 * log(4 * pi * lp * lc / 3) - 3 * d^2 / (4 * lp * lc) + log(1 - w)
  return(S)
}

# Generar los datos para j - i = 0, 1, 2 con los dos términos de S
data <- do.call(rbind, lapply(0:2, function(delta_ji) {
  # Para el primer término de S
  r_values <- seq(0, (delta_ji + 2) * 3.8, length.out = 100)
  S1_values <- sapply(r_values, calc_S1, delta_ji = delta_ji)
  
  # Para el segundo término de S
  d_values <- seq(0, (delta_ji + 2) * 3.8, length.out = 100)
  S2_values <- sapply(d_values, calc_S2, delta_ji = delta_ji)
  
  # Combinar datos
  data.frame(
    r = c(r_values, d_values),
    S = c(S1_values, S2_values),
    delta_ji = factor(delta_ji),
    type = rep(c("j-i=1", "j-i=2"), each = 100)
  )
}))

# Graficar usando ggplot
ggplot(data, aes(x = r, y = S, color = delta_ji, linetype = type)) +
  geom_line(size = 1) +
  labs(x = "L=(j-i+2)a", y = "S",
       color = "Ooka j-i", linetype = "Zhou") +
  theme_minimal() +
  ggtitle("Gráfico de los dos términos de S en función de la distancia para diferentes valores de j - i")

