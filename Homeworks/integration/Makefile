CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lm

DATA = erf_data.dat erf_accuracy.dat
PLOTS = erf_comparison.svg erf_diff.svg erf_accuracy.svg

all: $(PLOTS)

$(PLOTS): $(DATA) plot.gpi
	gnuplot plot.gpi

Out.txt $(DATA): main
	./main > Out.txt

main: main.o integrate.o
	$(CXX) main.o integrate.o -o main $(LDLIBS)

main.o: main.cc integrate.h
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

integrate.o: integrate.cc integrate.h
	$(CXX) $(CXXFLAGS) -c integrate.cc -o integrate.o

.PHONY: all clean test

clean:
	$(RM) main Out.txt *.o $(DATA) $(PLOTS)

test:
	@echo $(CXX)
	@echo $(RM)