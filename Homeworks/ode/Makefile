CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra

all: main pendulum.svg orbits.svg figure8.svg

main: main.o solver.o stepper.o matrix.o
	$(CXX) $(CXXFLAGS) $^ -o $@

pendulum.txt 1.txt 2.txt 3.txt 8.txt Out.txt: main
	./main > Out.txt

pendulum.svg: pendulum.txt plotA.gpi
	gnuplot plotA.gpi

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
	rm -f *.o main *.txt *.svg Out.txt