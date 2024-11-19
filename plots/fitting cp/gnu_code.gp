# Configuración de Gnuplot
set terminal pngcairo size 800,600
set output 'fit_cp_curve.png'
set xlabel "Temperatura (K)"
set ylabel "Cp (KJ/mol)"
set grid

# Definición de constantes
R = 8.314/1000  # Constante de los gases en KJ/(mol·K)

# Definición de la función para Cp(T)
#Cp(T) = (DeltaH * exp((DeltaH * (1 - T / Tm)) / (R * T))) / (R * Tm * (1 + exp((DeltaH * (1 - T / Tm)) / (R * T)))**2) +  DeltaCp*(1/(1+exp(-(T-Tm))))
Cp(T) = a*exp( -(T-Tm)**2/(2*DeltaH**2) ) +  DeltaCp*(1/(1+exp(-(T-Tm))))

# Valores iniciales de los parámetros para el ajuste
DeltaCp = 1.0
DeltaH = 10.0  # En KJ/mol
Tm = 360.0     # En K

# Realiza el ajuste a los datos
fit Cp(x) 'profthermo_fixed.dat' using 2:8 via DeltaCp, DeltaH, Tm,a

# Gráfico de Cp
plot 'profthermo_fixed.dat' using 2:8 with points pt 7 title 'Data', \
     Cp(x) with lines lw 2 lt rgb "red" title 'Fit'

###  p DH*(1-x/Tm)+Cp*(x-Tm-x*log(x/Tm)) lt rgb "blue" title "DeltaG", "profthermo.dat" u 2:8 lt rgb "red" title "Cp", 0 lt rgb "black" title ""
