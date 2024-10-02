.PHONY: iterative
iterative:
	clang -Wall -DITERATIVE_ALGO  -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp -O0 Strassen_Algorithm.c main.c Matrix.c -o run.out 

.PHONY: recursive
recursive:
	clang -Wall -DRECURSIVE_ALGO -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp -O0 Strassen_Algorithm.c main.c Matrix.c -o run.out 

.PHONY: parallel
parallel:
	clang -DPARALLEL_ALGO -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp -O0 Strassen_Algorithm.c main.c Matrix.c -o run.out 
 