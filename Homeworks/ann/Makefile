CXX = g++
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra

EXE = main
OBJ = main.o ann.o minimize.o qr.o matrix.o

all: Out.txt plots

Out.txt: $(EXE)
	./$(EXE) > Out.txt

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXE)

main.o: main.cc ann.h matrix.h
	$(CXX) $(CXXFLAGS) -c main.cc

ann.o: ann.cc ann.h minimize.h matrix.h
	$(CXX) $(CXXFLAGS) -c ann.cc

minimize.o: minimize.cc minimize.h qr.h matrix.h
	$(CXX) $(CXXFLAGS) -c minimize.cc

qr.o: qr.cc qr.h matrix.h
	$(CXX) $(CXXFLAGS) -c qr.cc

matrix.o: matrix.cc matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

plots: ann.tsv plot.gpi
	gnuplot plot.gpi

ann.tsv: Out.txt

clean:
	rm -f *.o $(EXE) Out.txt ann.tsv *.sv