unparallel:
	clang lab1_simple.cpp -O3 -Wall -o unparallel.out

run_unparallel:
	time ./unparallel.out > output.txt

parallel:
	mpicc lab1_parallel.cpp -O3 -o parallel.out

run_parallel:
	time mpirun -n 8 parallel.out

run_kernel:
	time mpirun -n 1 parallel.out

clean:
	rm output.txt *out
