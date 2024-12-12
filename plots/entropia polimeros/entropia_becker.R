library(ggplot2)

# Definir la función Q(r)
Q <- function(r, kappa) {
  a <- 14.054
  b <- 0.473
  c <- 1 - (1 + (0.38 * kappa^(-0.95))^-5)^(-1/5)
  d <- ifelse(kappa < 1/8, 0, 1 / (0.177 / (kappa - 0.111) + 6.40 * (kappa - 0.111)^0.783))
  cij <- matrix(
    c(-3/4, 23/64, -7/64, 
      -1/2, 17/16, -9/16), 
    nrow = 2, byrow = TRUE
  )
  
  J_SYD <- ifelse(kappa <= 1/8, (3/4 * pi * kappa)^(3/2) * (1 - 5 * kappa / 4), 112.04 * kappa^2 * exp(0.246 / kappa - a * kappa))
  term1 <- ((1 - c * r^2) / (1 - r^2))^(5/2)
  
  term2 <- 0
  i_values <- c(-1, 0)
  j_values <- 1:3
  for (i in seq_along(i_values)) {
    for (j in seq_along(j_values)) {
      term2 <- term2 + cij[i, j] * kappa^i_values[i] * r^(2 * j_values[j])
    }
  }
  term2 <- exp(term2 / (1 - r^2))
  
  term3 <- exp(-d * kappa * a * b * (1 + b) * r^2 / (1 - b^2 * r^2))
  term4 <- besselI(-d * kappa * a * (1 + b) * r / (1 - b^2 * r^2), nu = 0)
  Q_r <- J_SYD * term1 * term2 * term3 * term4
  return(Q_r)
}

# Valores de delta_ij y sus respectivas kappas
delta_ij <- c(2,5,15,25,60,100)
Lc <- (delta_ij+2) * 3.8 # Longitud de contorno (a = 3.8)
Lp <- 3.04 # Longitud de persistencia
kappa_values <- Lp / Lc

# Generar los valores de d y calcular Q(d) para cada kappa
d_values <- seq(from = 0, to = 1, by = 0.02) # d está normalizado por lc
data_list <- lapply(seq_along(kappa_values), function(i) {
  kappa <- kappa_values[i]
  Q_values <- -log(sapply(d_values, Q, kappa = kappa))
  data.frame(d = d_values, Q = Q_values, delta_ij = delta_ij[i])
})

# Combinar los datos en un único data frame
data <- do.call(rbind, data_list)

# Graficar con ggplot
ggplot(data, aes(x = d, y = Q, color = factor(delta_ij))) +
  geom_line(size = 1) +
  labs(
    title="Fixed delta_ij",
    x = "r/lc",
    y = "log Q(d)",
    color = "Delta_ij"
  ) +
  theme_minimal()
