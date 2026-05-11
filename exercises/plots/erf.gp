set terminal pngcairo size 900,600
set output "erf.png"

set title "Error function"
set xlabel "x"
set ylabel "erf(x)"

set grid
set key top left

plot "erf.data" using 1:2 with lines title "erf approximation", \
     "erf_tab.data" using 1:2 with points pointtype 7 pointsize 1.4 title "tabulated values"