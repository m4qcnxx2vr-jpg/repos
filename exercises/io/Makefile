CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -pedantic

PROGRAM = main
SOURCES = main.cc table.cc

all: $(PROGRAM)
	./$(PROGRAM) -n 1 -n 2 -n 3 -n 4 -n 5 > Out.txt

$(PROGRAM): $(SOURCES) table.h
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(PROGRAM)

stdin: $(PROGRAM)
	./$(PROGRAM) < Input.txt > Out_stdin.txt

file: $(PROGRAM)
	./$(PROGRAM) --input FileInput.txt --output FileOutput.txt

clean:
	rm -f $(PROGRAM) Out.txt Out_stdin.txt FileOutput.txt