CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O 
LDLIBS += -lstdc++  -lm
LINK.o = $(CXX) $(LDFLAGS)

Out.txt : main
	./main > Out.txt

main : main.o hello.o sfuns.o
#	$(CXX) $(LDFLAGS) main.o hello.o sfuns.o $(LDLIBS) -o main 
#	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o main

main.o : hello.h sfuns.h main.cc
#	$(CXX) $(CXXGLAFS) -c main.cc -o main.o 

hello.o: hello.h hello.cc
#	$(CXX) $(CXXFLAGS) -c hello.cc -o hello.o
sfuns.o: sfuns.h sfuns.cc
#	$(CXX) $(CXXFLAGS) -c sfuns.cc -o sfuns.o


.PHONEY:clean
clean:
	$(RM) main main.o Out.txt *.o


test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
	echo $(SHELL)
