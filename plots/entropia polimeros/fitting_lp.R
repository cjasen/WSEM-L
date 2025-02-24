library(ggplot2)
library(dplyr)
library(tidyr)

calc_Q <- function(d, delta_ji,lp) {
  a <- 14.054
  b <- 0.473
  Lc <- (delta_ji + 2) * 3.8
  kappa <- lp / Lc
  d <- d/Lc
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

# Params
L <- 13
lp <- 2.1 

# Valores de r
r_values <- seq(7, 35, by = 1)

# Calcula p(r)
p_values <- 4 * pi * r_values^2 * sapply(r_values, function(r) calc_Q(r, L, lp))

# Leer los datos experimentales
zhou_data <- read.table("zhou_data.txt", header = FALSE, col.names = c("r", "p"))
zhou_data$r <- as.numeric(zhou_data$r)
zhou_data$p <- as.numeric(zhou_data$p)

# Rango de los valores calculados
min_calculated <- min(p_values)
max_calculated <- max(p_values)

# Rango de los valores de zhou_data
min_zhou <- min(zhou_data$p)
max_zhou <- max(zhou_data$p)

# Transformación lineal
zhou_data$p_transformed <- min_calculated + 
  (max_calculated - min_calculated) * (zhou_data$p - min_zhou) / (max_zhou - min_zhou)
# Crea un dataframe para ggplot
df <- data.frame(r = r_values, p = p_values)

ggplot() +
  geom_line(data = df, aes(x = r, y = p), color = "blue", size = 1, linetype = "solid") +
  geom_point(data = zhou_data, aes(x = r, y = p_transformed), color = "red", size = 2) +
  labs(
       x = "r (A)",
       y = "p(r,lp)",
       color = "Datos") +
  theme_minimal(base_size = 18)