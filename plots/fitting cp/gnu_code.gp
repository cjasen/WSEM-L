# Configuración de Gnuplot
set terminal pngcairo size 800,600
set output 'fit_cp_curve.png'
set xlabel "T (K)"
set ylabel "Cp (KJ/mol)"
set grid

# Definición de constantes
R = 8.314/1000  # Constante de los gases en KJ/(mol·K)

# Definición de la función para Cp(T)
Cp(T) = a*exp( -(T-Tm)**2/(2*s**2) ) + DeltaCp*(1/(1+exp(-(T-Tm)))) + b

# Valores iniciales de los parámetros para el ajuste
a=20.0
b=3.0
DeltaCp = 1.0
s = 10.0  # En KJ/mol
Tm = 330.0     # En K

# Realiza el ajuste a los datos
fit Cp(x) 'profthermo.dat' using 2:8 via DeltaCp, s, Tm,a,b

# Gráfico de Cp
plot 'profthermo.dat' using 2:8 with points pt 7 title 'Data', \
     Cp(x) with lines lw 2 lt rgb "red" title 'Fit'

#p DeltaH*(1-x/Tm)+DeltaCp*(x-Tm-x*log(x/Tm)) lt rgb "blue" title "DeltaG", "profthermo.dat" u 2:8 lt rgb "red" title "Cp", 0 lt rgb "black" title ""
