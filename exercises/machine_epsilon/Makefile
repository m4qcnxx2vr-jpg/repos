CXX = c++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2

PROGRAM = main
OBJECTS = main.o epsilon.o compare.o

all: run

$(PROGRAM): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(PROGRAM)

main.o: main.cc epsilon.h compare.h
	$(CXX) $(CXXFLAGS) -c main.cc

epsilon.o: epsilon.cc epsilon.h
	$(CXX) $(CXXFLAGS) -c epsilon.cc

compare.o: compare.cc compare.h
	$(CXX) $(CXXFLAGS) -c compare.cc

run: $(PROGRAM)
	./$(PROGRAM) > Out.txt

clean:
	rm -f $(PROGRAM) $(OBJECTS) Out.txt