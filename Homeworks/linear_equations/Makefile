CXXFLAGS += -std=c++23 -Wall -Wextra
LDLIBS   += -lm

# Default target
all: Out.txt qr_scaling.svg

main: main.o matrix.o qr.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.cc matrix.h qr.h
matrix.o: matrix.cc matrix.h
qr.o: qr.cc qr.h matrix.h

# Run program and save output
Out.txt: main
	./$< > $@

# Timing measurements
out.times.data: main
	>$@
	for N in $$(seq 100 50 600); do  \
		gtime --format "$$N %e" --output $@ --append \
		./$< -size $$N 1>out 2>err ;\
	done

# Plot runtime scaling
qr_scaling.svg: out.times.data plot.gpi
	gnuplot plot.gpi


.PHONY: clean

clean:
	$(RM) [Oo]ut* main *.o *.svg err fit.log