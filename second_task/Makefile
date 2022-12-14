task1: task1.out
	mpirun -n 16 --oversubscribe --allow-run-as-root ./task1.out
task1.out: task1.c
	mpicc task1.c -o task1.out

task2: task2.out
	mpirun -n 5 --map-by :OVERSUBSCRIBE --allow-run-as-root ./task2.out
task2.out: task2.c
	mpicc task2.c -o task2.out

clean:
	rm -rf *.out
