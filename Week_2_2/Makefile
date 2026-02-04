#CXX = g++ -std=c++23
CXX = g++ -std=c++23
CXXLAGS += -Wall -Wextra #-Werror -o
LDLIBS += -lstdc++ -lm

all: Out.txt
	cat $<
Out.txt: main
#	./main > Out.txt
	./$< > $@

main: main.o hello.o
main.o = hello.h


