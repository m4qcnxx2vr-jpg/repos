CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

TARGET = main
SRC = main.cc integrate.cc

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

run: $(TARGET)
	./$(TARGET) > Out.txt
	@echo "Output written to Out.txt"
	@cat Out.txt

plot: Out.txt
	gnuplot plot.gpi

clean:
	rm -f $(TARGET) Out.txt erf_data.dat erf_accuracy.dat \
	      erf_comparison.svg erf_diff.svg erf_accuracy.svg

.PHONY: all run plot clean