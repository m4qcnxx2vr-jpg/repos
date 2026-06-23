CXX = c++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2

TARGET = main
OBJECTS = main.o mc.o

all: $(TARGET)
	./$(TARGET) > Out.txt
	gnuplot plotA2.gpi

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

main.o: main.cc mc.h
	$(CXX) $(CXXFLAGS) -c main.cc

mc.o: mc.cc mc.h
	$(CXX) $(CXXFLAGS) -c mc.cc

clean:
	rm -f $(TARGET) $(OBJECTS) Out.txt circle_errors.txt *.png