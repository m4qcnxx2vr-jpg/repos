CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm

TIME := $(shell command -v gtime >/dev/null 2>&1 && echo gtime || echo /usr/bin/time)

all : Out.txt out.times.svg

Out.txt : main
	./main > Out.txt

out.times.data : main
	>$@
	for N in $$(seq 5 5 200); do \
		$(TIME) --format "$$N %e" --output $@ --append \
		./$< -size $$N 1>out 2>err ;\
	done

out.times.svg : out.times.data timing.gpi
	gnuplot timing.gpi

main : main.o matrix.o svd.o
	$(CXX) main.o matrix.o svd.o -o main $(LDLIBS)

main.o : matrix.h svd.h main.cc
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

matrix.o : matrix.h matrix.cc
	$(CXX) $(CXXFLAGS) -c matrix.cc -o matrix.o

svd.o : matrix.h svd.h svd.cc
	$(CXX) $(CXXFLAGS) -c svd.cc -o svd.o

.PHONY: all clean test

clean:
	$(RM) main main.o matrix.o svd.o Out.txt out.times.data out.times.svg fit.log out err *.o