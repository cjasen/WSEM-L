# Configura el terminal y salida
set terminal pngcairo size 1000,800
set output 'graficos_compartidos.png'

# Activa el modo multiplot con layout
set multiplot layout 2,1 title "Gráficos de m, sigma y Cp vs T"

# Configura el eje X solo una vez y activa la cuadrícula
set xlabel 'T (k)'
set grid

# Gráfico 1: m v. T y sigma v. T
set key outside
set ylabel 'Valores de m y sigma'
plot 'profthermo.dat' using 2:3 with lines title 'm' linecolor rgb 'purple' lw 2, \
     'profthermo.dat' using 2:4 with lines title 'sigma' linecolor rgb 'blue' lw 2

# Gráfico 2: Cp v. T y sigma v. T
# Calcula el valor máximo de Cp (columna 8)
stats 'profthermo.dat' using 8 nooutput
max_Cp = STATS_max

# Segundo gráfico con el mismo eje X
set ylabel 'Valores de Cp y sigma'
plot 'profthermo.dat' using 2:($8/max_Cp + 0.5) with lines title 'Cp' linecolor rgb 'red' lw 2, \
     'profthermo.dat' using 2:4 with lines title 'sigma' linecolor rgb 'blue' lw 2

# Desactiva el modo multiplot
unset multiplot
