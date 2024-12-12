library(ggplot2)
library(dplyr)
library(tidyr)

# Parámetros
lp <- 3.04  # Longitud de persistencia en Amstrongs

# Función S1 (Ooka)
calc_S1 <- function(d, delta_ji) {
  lc <- (delta_ji + 2) * 3.8
  d <- d * lc
  S <- 1.5 * log(delta_ji + 2) + 1.5 * (d^2 - 3.8^2) / (lc * 2 * 20)
  return(S)
}

# Función S2 (Zhou)
calc_S2 <- function(d, delta_ji) {
  lc <- (delta_ji + 2) * 3.8
  d <- d * lc
  w <- (5.0 * lp / (4.0 * lc)) - (2.0 * d^2 / lc^2) +
    (33.0 * d^4 / (80.0 * lp * lc^3)) + 
    (79.0 * lp^2 / (160.0 * lc^2)) + 
    (329.0 * d^2 * lp / (120.0 * lc^3)) - 
    (6799.0 * d^4 / (1600.0 * lc^4)) + 
    (3441.0 * d^6 / (2800.0 * lp * lc^5)) - 
    (1089.0 * d^8 / (12800.0 * lp^2 * lc^6))
  S <- 1.5 * log(4 * pi * lp * lc / 3) + 3 * d^2 / (4 * lp * lc) - log(1 - w)
  return(S)
}

# Función S3 (Becker) [S3=-logQ]
calc_Q <- function(d, delta_ji) {
  a <- 14.054
  b <- 0.473
  Lc <- (delta_ji + 2) * 3.8
  kappa <- lp / Lc #lp here is not the same as in other functions
  c <- 1 - (1 + (0.38 * kappa^(-0.95))^-5)^(-1/5)
  d_val <- ifelse(kappa < 1/8, 0, 1 / (0.177 / (kappa - 0.111) + 6.40 * (kappa - 0.111)^0.783))
  
  cij <- matrix(
    c(-3/4, 23/64, -7/64, 
      -1/2, 17/16, -9/16), 
    nrow = 2, byrow = TRUE
  )
  
  J_SYD <- ifelse(kappa <= 1/8, (3/4 * pi * kappa)^(3/2) * (1 - 5 * kappa / 4), 
                  112.04 * kappa^2 * exp(0.246 / kappa - a * kappa))
  
  term1 <- ((1 - c * d^2) / (1 - d^2))^(5/2)
  term2 <- 0
  i_values <- c(-1, 0)
  j_values <- 1:3
  for (i in seq_along(i_values)) {
    for (j in seq_along(j_values)) {
      term2 <- term2 + cij[i, j] * kappa^i_values[i] * d^(2 * j_values[j])
    }
  }
  term2 <- exp(term2 / (1 - d^2))
  term3 <- exp(-d_val * kappa * a * b * (1 + b) * d^2 / (1 - b^2 * d^2))
  term4 <- besselI(-d_val * kappa * a * (1 + b) * d / (1 - b^2 * d^2), nu = 0)
  
  Q_r <- J_SYD * term1 * term2 * term3 * term4
  return(Q_r)
}

# Valores de delta_ji y rango de d
delta_ji <- c(2, 5, 20, 70)
d <- seq(0, 1, 0.01)

# Calcular las tres entropías
data <- expand.grid(d = d, delta_ji = delta_ji) %>%
  mutate(Ooka = mapply(calc_S1, d, delta_ji),
         Zhou = mapply(calc_S2, d, delta_ji),
         Becker = -log(mapply(calc_Q, d, delta_ji)))

# Convertir los datos a formato largo
data_long <- data %>%
  gather(key = "S_type", value = "S_value", Ooka, Zhou, Becker)

# Graficar
ggplot(data_long, aes(x = d, y = S_value, color = factor(delta_ji), linetype = S_type)) +
  geom_line(size = 1.2) + 
  labs(title = "fixed delta_ij", 
       x = "r/lc", 
       y = "Entropía (S)", 
       color = "delta_ij", 
       linetype = "Tipo de S") +
  theme_minimal()


################# FIXED DISTANCE


d <- c(0.1, 0.3 , 0.6,0.85)
delta_ji <- seq(0,120,1)

data2 <- expand.grid(d = d, delta_ji = delta_ji) %>%
  mutate(Ooka = mapply(calc_S1, d, delta_ji),
         Zhou = mapply(calc_S2, d, delta_ji),
         Becker = -log(mapply(calc_Q, d, delta_ji)))

data_long2 <- data2 %>%
  gather(key = "S_type", value = "S_value", Ooka, Zhou, Becker)

ggplot(data_long2, aes(x = delta_ji, y = S_value, color = factor(d), linetype = S_type)) +
  geom_line(size = 1.2) + 
  labs(title = "fixed d (r_ij)", 
       x = "delta_ij", 
       y = "S", 
       color = "d/lc", 
       linetype = "S") +
  theme_minimal()

