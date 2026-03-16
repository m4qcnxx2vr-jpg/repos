CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

OBJS = main.o matrix.o qr.o lsfit.o

all: plot

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o main $(OBJS)

main.o: main.cc matrix.h qr.h lsfit.h
	$(CXX) $(CXXFLAGS) -c main.cc

matrix.o: matrix.cc matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

qr.o: qr.cc qr.h matrix.h
	$(CXX) $(CXXFLAGS) -c qr.cc

lsfit.o: lsfit.cc lsfit.h qr.h matrix.h
	$(CXX) $(CXXFLAGS) -c lsfit.cc

run: main
	./main > Out.txt

plot: run decay_plot.gpi
	gnuplot decay_plot.gpi

clean:
	rm -f *.o main Out.txt decay_data.txt decay.svg fit_params.gpi