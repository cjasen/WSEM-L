library(fields)

# Definir parámetros
a <- 3.8
A <- 20
k_B <- 1 #1.38e-23  # Constante de Boltzmann en Joules/Kelvin

# Leer el archivo rCalpha.txt que contiene i, j y r_ij
data <- read.table("rCalpha.txt", header=FALSE)
colnames(data) <- c("i", "j", "r_ij")

# Crear una matriz vacía para S_ij
n <- 120
S_matrix <- matrix(NA, n, n)  # Usamos NA en vez de 0 para identificar valores no calculados

# Función para calcular la entropía S'(L)
entropy_Ooka <- function(L, r, a, A, k_B) {
  term1 <- log(L)
  term2 <- (r^2 - a^2) / (2 * A * a * L)
  return(-3/2 * k_B * (term1 + term2))
}

# Llenar la matriz S_ij
for (row in 1:nrow(data)) {
  i <- data[row, "i"]
  j <- data[row, "j"]
  r_ij <- data[row, "r_ij"]
  
  # Verificar que i y j estén dentro de los límites y que i < j
  if (i < j && i >= 1 && j <= n) {
    L <- j - i
    S_matrix[i, j] <- entropy_Ooka(L, r_ij, a, A, k_B)
  }
}

# Reemplazar los NA por ceros (opcional, si prefieres)
S_matrix[is.na(S_matrix)] <- 0

# Graficar la matriz S_ij con leyenda
image.plot(1:n, 1:n, t(S_matrix), col=topo.colors(100), xlab="i", ylab="j", 
           main="S_ij (i<j) with Kb=1")
grid()

