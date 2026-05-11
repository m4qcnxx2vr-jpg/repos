set terminal pngcairo size 900,600
set output "lngamma.png"

set title "Logarithm of the Gamma function"
set xlabel "x"
set ylabel "ln Gamma(x)"

set xrange [0:10]
set yrange [-1:15]

set zeroaxis
set grid
set key top left

plot "lngamma.data" using 1:2 with lines title "lngamma approximation", \
     "lngamma_tab.data" using 1:2 with points pointtype 7 pointsize 1.4 title "ln((n-1)!)"