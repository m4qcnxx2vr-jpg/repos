CXXFLAGS += -std=c++23 -Wall -O
LDLIBS += -lstdc++ -lm -lpthread
SHELL := /bin/bash
TIME = time --format "$$nthreads %e %U" --append --output $@

times.svg: out.times
	gnuplot <<'GP' | tee log.times
	set terminal svg
	set output "$@"
	set xlabel "number of threads"
	set ylabel "real time (seconds)"
	a=1; b=10 ;\
	f(x) = a + b/x
	fit f(x) "$<" using 1:2 via a,b ;\
	plot \
	 	"$<" using 1:2 with linespoints title "times", \
		f(x) title "fit with 1/n";\
	' | tee log.times | gnuplot

N=1e9

out.times: main Makefile
	: > $@
	: > log.main
	for nthreads in 1 2 3 4 5 6 7 8 9 10; do \
		$(TIME) ./main -terms $(N) -threads $$nthreads >> log.main ; \
	done
