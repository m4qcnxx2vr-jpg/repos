CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra

all: main orbits.svg figure8.svg

main: main.o solver.o stepper.o matrix.o
	$(CXX) $(CXXFLAGS) $^ -o $@

1.txt 2.txt 3.txt 8.txt Out.txt: main
	./main

orbits.svg: 1.txt 2.txt 3.txt plot.gpi
	gnuplot plot.gpi

figure8.svg: 8.txt plot_8.gpi
	gnuplot plot_8.gpi

main.o: main.cc solver.h stepper.h matrix.h
	$(CXX) $(CXXFLAGS) -c main.cc

solver.o: solver.cc solver.h
	$(CXX) $(CXXFLAGS) -c solver.cc

stepper.o: stepper.cc stepper.h
	$(CXX) $(CXXFLAGS) -c stepper.cc

matrix.o: matrix.cc matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

clean:
	rm -f *.o main *.txt orbits.svg