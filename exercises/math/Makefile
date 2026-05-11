CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2

all: Out.txt

Out.txt: main
	./main > Out.txt

main: main.o sfuns.o
	$(CXX) $(CXXFLAGS) -o main main.o sfuns.o

main.o: main.cc sfuns.h
	$(CXX) $(CXXFLAGS) -o main.o -c main.cc

sfuns.o: sfuns.cc sfuns.h
	$(CXX) $(CXXFLAGS) -o sfuns.o -c sfuns.cc

clean:
	rm -f main *.o Out.txt