CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

OBJS = main.o jacobi.o matrix.o

all: run

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o main $(OBJS)

main.o: main.cc matrix.h jacobi.h
	$(CXX) $(CXXFLAGS) -c main.cc

jacobi.o: jacobi.cc jacobi.h matrix.h
	$(CXX) $(CXXFLAGS) -c jacobi.cc

matrix.o: matrix.cc matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cc

run: main
	./main > Out.txt

clean:
	rm -f *.o main Out.txt