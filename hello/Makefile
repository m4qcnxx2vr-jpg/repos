CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O 
LDLIBS = -lstdc++  -lm

Out.txt : main
	./main > Out.txt

main : main.o hello.o
	$(CXX) main.o hello.o -o main $(LDLIBS)

main.0 : hello.h main.cc
	$(CXX) $(CXXGLAFS) -c main.cc -o main.o 

hello.o: hello.h hello.cc
	$(CXX) $(CXXFLAGS) -c hello.cc -o hello.o

.PHONEY:clean
clean:
	$(RM) main main.o Out.txt *.o


test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
