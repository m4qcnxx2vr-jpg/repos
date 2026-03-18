CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

OBJS = main.o linterp.o qspline.o fqinterp.o

all: plots

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o main $(OBJS)

main.o: main.cc linterp.h qspline.h fqinterp.h
	$(CXX) $(CXXFLAGS) -c main.cc

linterp.o: linterp.cc linterp.h
	$(CXX) $(CXXFLAGS) -c linterp.cc

qspline.o: qspline.cc qspline.h
	$(CXX) $(CXXFLAGS) -c qspline.cc

fqinterp.o: fqinterp.cc fqinterp.h
	$(CXX) $(CXXFLAGS) -c fqinterp.cc

run: main
	./main > Out.txt
	@echo "----- Output -----"
	@cat Out.txt

plots: run
	gnuplot taskA_plot.gpi
	gnuplot taskB_plot.gpi
	gnuplot taskC_plot.gpi

clean:
	rm -f *.o main Out.txt \
	linspline_data.txt linspline_curve.txt taskA.svg \
	qspline_data.txt qspline_curve.txt taskB.svg \
	fqinterp_curve.txt taskC.svg