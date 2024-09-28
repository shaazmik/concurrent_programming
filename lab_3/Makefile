simple:
	clang -Wall -fsanitize=address Strassen_Algorithm.c simple.c -o simple.out

parallel:
	clang -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp parallel.c -o parallel 
 