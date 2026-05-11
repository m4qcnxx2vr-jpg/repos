CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -pedantic

PROGRAM = main
SOURCES = main.cc

all: $(PROGRAM)
	./$(PROGRAM) < Input.txt > Out.txt

$(PROGRAM): $(SOURCES) genlist.h
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(PROGRAM)

clean:
	rm -f $(PROGRAM) Out.txt