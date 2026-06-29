CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lm

DATA = erf_data.dat erf_accuracy.dat
PLOTS = erf_comparison.svg erf_diff.svg erf_accuracy.svg

all: Out.txt PartB.txt $(PLOTS)

$(PLOTS): $(DATA) plot.gpi
	gnuplot plot.gpi

Out.txt $(DATA): main
	./main > Out.txt

PartB.txt: test_partB
	./test_partB > PartB.txt

main: main.o integrate.o
	$(CXX) main.o integrate.o -o main $(LDLIBS)

test_partB: test_partB.o integrate.o
	$(CXX) test_partB.o integrate.o -o test_partB $(LDLIBS)

main.o: main.cc integrate.h
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

test_partB.o: test_partB.cc integrate.h
	$(CXX) $(CXXFLAGS) -c test_partB.cc -o test_partB.o

integrate.o: integrate.cc integrate.h
	$(CXX) $(CXXFLAGS) -c integrate.cc -o integrate.o

.PHONY: all clean test

clean:
	$(RM) main test_partB Out.txt PartB.txt *.o $(DATA) $(PLOTS)

test:
	@echo $(CXX)
	@echo $(RM)