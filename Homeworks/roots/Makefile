CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm

GPI_FILES = wavefunction.gpi conv_rmax.gpi conv_rmin.gpi conv_acc.gpi conv_eps.gpi
PNG_FILES = wavefunction.png conv_rmax.png conv_rmin.png conv_acc.png conv_eps.png
DATA_FILES = wavefunction.txt conv_rmax.txt conv_rmin.txt conv_acc.txt conv_eps.txt

all : Out.txt $(PNG_FILES)

Out.txt : main
	./main > Out.txt

main : main.o matrix.o qr.o roots.o
	$(CXX) main.o matrix.o qr.o roots.o -o main $(LDLIBS)

main.o : roots.h matrix.h main.cc
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

matrix.o : matrix.h matrix.cc
	$(CXX) $(CXXFLAGS) -c matrix.cc -o matrix.o

qr.o : qr.h matrix.h qr.cc
	$(CXX) $(CXXFLAGS) -c qr.cc -o qr.o

roots.o : roots.h qr.h matrix.h roots.cc
	$(CXX) $(CXXFLAGS) -c roots.cc -o roots.o

$(DATA_FILES) : Out.txt

$(PNG_FILES) : %.png : %.gpi %.txt
	gnuplot $<

.PHONY: clean all
clean:
	$(RM) main Out.txt *.o $(DATA_FILES) $(PNG_FILES)
