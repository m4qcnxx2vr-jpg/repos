CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23

Out.txt: main
	./main > Out.txt


main: main.o approx.o
	$(CXX) $(CXXFLAGS) -o main main.o approx.o

main.o: main.cc approx.h
	$(CXX) $(CXXFLAGS) -c main.cc

approx.o: approx.cc approx.h
	$(CXX) $(CXXFLAGS) -c approx.cc

clean:
	rm -f *.o main Out.txt