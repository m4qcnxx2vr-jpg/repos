CXXFLAGS += -std=c++23 -Wall -O
LDLIBS += -lstdc++ -lm -lpthread
SHELL := /bin/bash
TIME = gtime --format "$$nthreads %e %U" --append --output $@


times.svg: out.times Makefile
	echo '\
	set terminal svg ;\
	set output "times.svg" ;\
	set xlabel "Threads" ;\
	set ylabel "Real time (seconds)" ;\
	a=1; b=10 ;\
	f(x)=a+b/x ;\
	fit f(x) "$<" using 1:2 via a,b ;\
	plot "$<" using 1:2 with linespoints title "times", \
	     f(x) title "fit with 1/n" ;\
	' | tee log.times | gnuplot


N=1e9

out.times: main Makefile
	: > $@
	: > log.main
	for nthreads in $$(seq 8); do \
		$(TIME) ./main -terms $(N) -threads $$nthreads >> log.main ; \
	done

all: out1 out2 out3 out4 
out1: main ; ./main -threads 1 -terms 1e8 > $@ 
out2: main ; ./main -threads 1 -terms 2e8 > $@ 
out3: main ; ./main -threads 1 -terms 3e8 > $@ 
out4: main ; ./main -threads 1 -terms 4e8 > $@ 

test: main 
	./main ; ./main -threads 1 -terms 1e8 > out5 &
	./main ; ./main -threads 1 -terms 2e8 > out6 & 
	./main ; ./main -threads 1 -terms 3e8 > out7 & 
	./main ; ./main -threads 1 -terms 4e8 > out8 & 

clean:
	$(RM) log* [Oo]ut* main *.svg fit.log