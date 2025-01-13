# Configura el terminal y salida
set terminal pngcairo size 1000,800
set output 'graficos.png'

# Activa el modo multiplot
set multiplot layout 2,1

# Gráfico 1: m v. T y sigma v. T
set xlabel 'T (K)'
set ylabel ''
set grid
set key outside
plot 'profthermo.dat' using 2:3 with lines title '<m>' linecolor rgb 'purple' lw 2, \
     'profthermo.dat' using 2:4 with lines title '<{/Symbol s}>' linecolor rgb 'blue' lw 2

# Gráfico 2: Cp v. T y sigma v. T
# Calcula el valor máximo de Cp (columna 8)
stats 'profthermo.dat' using 8 nooutput
max_Cp = STATS_max

# Escala Cp y ajusta la posición sumando 0.5
set xlabel 'T (K)'
set ylabel 'KJ/(molK)'
plot 'profthermo.dat' using 2:8 with lines title 'Cp' linecolor rgb 'red' lw 2

# Desactiva el modo multiplot
unset multiplot