set terminal pngcairo size 900,600
set output "running_time.png"

set title "Harmonic sum running time"
set xlabel "Number of threads"
set ylabel "Time in seconds"

set grid
set key top right

plot "out.times" using 1:2 with linespoints title "real time", \
     "out.times" using 1:3 with linespoints title "user CPU time"