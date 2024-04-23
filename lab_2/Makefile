all: integral simple global

integral: integral_with_collective.cpp 
	clang++ integral_with_collective.cpp -lm -pthread -O2 -Wall -std=c++14 -o integral_with_collective.out 

simple: integral_simple.cpp
	clang++ integral_simple.cpp -lm -pthread -O2 -Wall -std=c++14 -o integral_simple.out 

global: integral_global_stack.cpp
	clang++ integral_global_stack.cpp -lm -pthread -O2 -Wall -std=c++14 -o integral_global_stack.out 


clean:
	rm -f *.out
