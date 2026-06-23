CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
GNUPLOT = gnuplot

COMMON_SRCS = ann.cc qr.cc matrix.cc
COMMON_OBJS = $(COMMON_SRCS:.cc=.o)


all: build run-main run-test plot

build: main test_ann

main: main.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

test_ann: test_ann.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

run-main: main
	./main

run-test: test_ann
	./test_ann > out.txt
	@echo "wrote out.txt"

plot: curve.txt plot.gpi
	@$(GNUPLOT) plot.gpi || \
		echo "(skipping plot: gnuplot not available -- install it and run 'make plot')"

curve.txt: main
	./main

clean:
	rm -f *.o main test_ann params.txt curve.txt out.txt plot.png

.PHONY: all build run-main run-test plot clean
