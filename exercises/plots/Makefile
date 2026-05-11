CXX = c++
CXXFLAGS = -std=c++23 -O3 -Wall -Wextra -pedantic
TARGET = main
OBJS = main.o SpecialFunctions.o

all: Out.txt erf.png gamma.png lngamma.png

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

main.o: main.cc SpecialFunctions.h
	$(CXX) $(CXXFLAGS) -c main.cc

SpecialFunctions.o: SpecialFunctions.cc SpecialFunctions.h
	$(CXX) $(CXXFLAGS) -c SpecialFunctions.cc

Out.txt: $(TARGET)
	./$(TARGET) > Out.txt

erf.png: Out.txt erf.gp
	gnuplot erf.gp

gamma.png: Out.txt gamma.gp
	gnuplot gamma.gp

lngamma.png: Out.txt lngamma.gp
	gnuplot lngamma.gp

clean:
	rm -f $(TARGET) $(OBJS) Out.txt
	rm -f erf.data erf_tab.data gamma.data gamma_tab.data
	rm -f lngamma.data lngamma_tab.data
	rm -f erf.png gamma.png lngamma.png

.PHONY: all clean