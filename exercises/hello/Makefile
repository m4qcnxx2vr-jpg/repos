CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -pedantic


Out.txt : hello
	./hello > Out.txt
hello: hello.cc
	c++ -o hello hello.cc

clean:
	rm -f hello Out.txt