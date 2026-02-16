CXXFLAGs += -Wall -std=c++23
LDLIBS += -lstdc++ -lm

fig1.svg: fig1.gpi erf.data
	gnuplot fig1.gpi

erf.data: main Makefile
	./$< -xmin -2 -xmax 2 dx 0.25> $@ 2>log.erf
	