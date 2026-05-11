set terminal pngcairo size 900,600
set output "gamma.png"

set title "Gamma function"
set xlabel "x"
set ylabel "Gamma(x)"

set xrange [0:6]
set yrange [0:25]

set zeroaxis
set grid
set key top right

plot "gamma.data" using 1:2 with lines title "Stirling approximation with recurrence", \
     "gamma_tab.data" using 1:2 with points pointtype 7 pointsize 1.4 title "factorials"