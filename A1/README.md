# Assignment-1

### Command for sequential code
	./compile.sh s
	./a.out matrix_size check file_path write

### Command for parallel code
	./compile.sh x
	./a.out matrix_size number_of_threads check file_path write

### Notations
1. x = o for openmp and p for pthreads
2. check = 1 if you want to print the L2,1 norm of PA-LU matrix
3. write = 1 if you want the matrices written to file