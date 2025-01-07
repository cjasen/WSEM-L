# Configura el terminal y salida
set terminal pngcairo size 1000,800
set output 'sim_vs_exp.png'


# Gráfico 1: m v. T y sigma v. T
set xlabel 'T'
set ylabel 'Cp'
set grid
set key outside
plot 'profthermo.dat' using 2:8 with lines title 'Sim' linecolor rgb 'red' lw 2, "2PHT_expCp_noSS.txt" lc rgb "blue" title "Exp"
