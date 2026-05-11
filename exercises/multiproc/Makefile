CXX = c++
CXXFLAGS = -std=c++23 -O3 -Wall -Wextra -pedantic
TARGET = main
OBJS = main.o Harmonic.o

NTERMS = 1000000000
THREADS = 1 2 3 4 5 6 7 8
REPEATS = 3

all: Out.txt out.times running_time.png

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

main.o: main.cc Harmonic.h
	$(CXX) $(CXXFLAGS) -c main.cc

Harmonic.o: Harmonic.cc Harmonic.h
	$(CXX) $(CXXFLAGS) -c Harmonic.cc

Out.txt: $(TARGET)
	./$(TARGET) -threads 8 -terms $(NTERMS) > Out.txt

out.times: $(TARGET)
	> out.times
	> log.main
	for nthreads in $(THREADS); do \
		for repeat in $$(seq 1 $(REPEATS)); do \
			gtime --format "$$nthreads %e %U" --append --output out.times \
			./$(TARGET) -threads $$nthreads -terms $(NTERMS) >> log.main ; \
		done \
	done

running_time.png: out.times plot_times.gp
	gnuplot plot_times.gp

clean:
	rm -f $(TARGET) $(OBJS) Out.txt out.times log.main running_time.png

.PHONY: all clean