CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm

GPI_FILES = higgs.gpi
PNG_FILES = higgs.png
DATA_FILES = higgs_fit.txt higgs_curve.txt

all : Out.txt $(PNG_FILES)

Out.txt : main higgs.data.txt
	./main < higgs.data.txt > Out.txt

main : main.o matrix.o qr.o minimize.o
	$(CXX) main.o matrix.o qr.o minimize.o -o main $(LDLIBS)

main.o : minimize.h matrix.h main.cc
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

matrix.o : matrix.h matrix.cc
	$(CXX) $(CXXFLAGS) -c matrix.cc -o matrix.o

qr.o : qr.h matrix.h qr.cc
	$(CXX) $(CXXFLAGS) -c qr.cc -o qr.o

minimize.o : minimize.h qr.h matrix.h minimize.cc
	$(CXX) $(CXXFLAGS) -c minimize.cc -o minimize.o

$(DATA_FILES) : Out.txt

$(PNG_FILES) : %.png : %.gpi $(DATA_FILES)
	gnuplot $<

.PHONY: clean all
clean:
	$(RM) main Out.txt *.o $(DATA_FILES) $(PNG_FILES)
