CXX = g++
CXXFLAGS = -O2 -std=c++23 -Wall -Wextra -pedantic

OBJS = main.o matrix.o jacobi.o

all: main Out.txt convergence_dr.svg convergence_rmax.svg wavefunctions.svg

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o main $(OBJS)

Out.txt convergence_dr.dat convergence_rmax.dat state_0.dat state_1.dat state_2.dat state_3.dat: main
	./main -rmax 10 -dr 0.3 -nstates 4 > Out.txt

convergence_dr.svg: convergence_dr.dat
	gnuplot plot_convergence_dr.gpi

convergence_rmax.svg: convergence_rmax.dat
	gnuplot plot_convergence_rmax.gpi

wavefunctions.svg: state_0.dat state_1.dat state_2.dat state_3.dat
	gnuplot plot_wavefunctions.gpi

main.o: main.cc matrix.h jacobi.h
	$(CXX) $(CXXFLAGS) -c main.cc

matrix.o: matrix.cc matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

jacobi.o: jacobi.cc jacobi.h matrix.h
	$(CXX) $(CXXFLAGS) -c jacobi.cc

clean:
	rm -f *.o main Out.txt *.dat *.svg