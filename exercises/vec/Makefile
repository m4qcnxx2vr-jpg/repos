CXX = c++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2

PROGRAM = main
OBJECTS = main.o vec.o

all: Out.txt

Out.txt: $(PROGRAM)
	./$(PROGRAM) > Out.txt

$(PROGRAM): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(PROGRAM)

main.o: main.cc vec.h
	$(CXX) $(CXXFLAGS) -c main.cc

vec.o: vec.cc vec.h
	$(CXX) $(CXXFLAGS) -c vec.cc

clean:
	rm -f $(PROGRAM) $(OBJECTS) Out.txt