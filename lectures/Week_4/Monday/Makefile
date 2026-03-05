CXXFLAGS += -Wall -std=c++23
LDLIBS   += -lstdc++ -lm

all: fig1.svg fig2.png fig3.svg

test.txt: erf.data Makefile
	>$@
	cat erf.data | while read a b c; do echo $$a $$c >> $@; done

fig3.svg: test.txt Makefile
	graph \
		-X "x" \
		-Y "y" \
		-T svg \
		-m 1 \
		$< > $@


fig2.png: erf.data Makefile
	echo '\
	set terminal png ;\
	set output "$@" ;\
	set tics out ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	plot [][-2:2] \
		"$<" index 1 using 1:2 with linespoints \
		,"$<" index 1 using 1:3 with lines ;\
	' | tee log.fig2.gpi | gnuplot

fig1.svg: fig1.gpi erf.data
	gnuplot fig1.gpi

erf.data: main Makefile
	./$< -xmin -2 -xmax 2 -dx 0.125 > $@ 2> log.erf
	echo "" >> $@
	echo "" >> $@
	./$< -xmin 0 -xmax 10 -dx 0.125 >> $@ 2>> log.erf