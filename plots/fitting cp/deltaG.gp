# Configuración de Gnuplot
set terminal pngcairo size 800,600
set output 'deltaG.png'
set xlabel "T (K)"
set ylabel "Cp and DeltaG (KJ/mol)"
set xrange [170:365]
set yrange [-15:40]
set grid


# Valores sacados del fit de Cp y la integración
DeltaH=290 # calculado integrando a mano Cp
DeltaCp=3.7
Tm=331.5

p DeltaH*(1-x/Tm)+DeltaCp*(x-Tm-x*log(x/Tm)) lt rgb "blue" title "DeltaG", "profthermo_fixed.dat" u 2:8 lt rgb "red" title "Cp" w lines, 0 lt rgb "black" title ""
