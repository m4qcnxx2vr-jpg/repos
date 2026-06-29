CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O2
LDLIBS = -lstdc++ -lm

TIME := $(shell command -v gtime >/dev/null 2>&1 && echo gtime || echo /usr/bin/time)

all: Out.txt

times: out.times.svg

comparison: svd_vs_jacobi_timing.svg

Out.txt: main
	./main > Out.txt

out.times.data: main
	>$@
	for N in $$(seq 5 5 200); do \
		$(TIME) --format "$$N %e" --output $@ --append \
		./$< -size $$N 1>out 2>err ; \
	done

out.times.svg: out.times.data timing.gpi
	gnuplot timing.gpi

svd.times.data: main
	>$@
	for N in $$(seq 5 5 200); do \
		$(TIME) --format "$$N %e" --output $@ --append \
		./$< -size $$N -method svd 1>out 2>err ; \
	done

jacobi.times.data: main
	>$@
	for N in $$(seq 5 5 200); do \
		$(TIME) --format "$$N %e" --output $@ --append \
		./$< -size $$N -method jacobi 1>out 2>err ; \
	done

svd_vs_jacobi_timing.svg: svd.times.data jacobi.times.data svd_vs_jacobi.gpi
	gnuplot svd_vs_jacobi.gpi

main: main.o matrix.o svd.o jacobi.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.cc matrix.h svd.h jacobi.h
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

matrix.o: matrix.cc matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cc -o matrix.o

svd.o: svd.cc svd.h matrix.h
	$(CXX) $(CXXFLAGS) -c svd.cc -o svd.o

jacobi.o: jacobi.cc jacobi.h matrix.h
	$(CXX) $(CXXFLAGS) -c jacobi.cc -o jacobi.o

.PHONY: all clean times comparison test

clean:
	$(RM) main *.o Out.txt \
	      out.times.data out.times.svg \
	      svd.times.data jacobi.times.data svd_vs_jacobi_timing.svg \
	      fit.log out err