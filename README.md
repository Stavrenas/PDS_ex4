# PDS_ex4
This project is part of the "Parallel and Distributed Systems" at Aristotle University of Thessaloniki. The aim of the code is fast boolean matrix multiplication using both MPI and OpenMP . The code production was unsupervised

In order to compile and run the executables you must have installed the following libraries:
* [OpenMP](https://www.openmp.org/)
* [Open MPI](https://www.open-mpi.org/)
* [make](https://www.gnu.org/software/make/)

## Compiling
* Type `make` to compile all versions; v0,v1,v2 and v3.
* Type `make v0` to compile v0 only, `make v1` to compile v1 etc..
* Type `make clean` to remove all unnecessary files.

## Note
Every version implements F.(A * B) where F == A == B.

## Running
* `v0` is a simple, single threaded BMM algorithm. The only argument is the name of the matrix , without the .mtx extension.
Example: `./v0 test` to run the algorithm with the matrix in the file test.mtx
* `v1` is a multithreaded BMM algorithm. As above, there is a single command line argument, the name of the matrix.
* `v2` is a single threaded BMM, using a block algorithm. In this case, the argument is only the name of the matrix. The blocksize is calculated automatically. For exmaple: `./v2 test`
* `v3` is an algorithm with 2 levels of parallelization. The arguments are the same as v2, but we execute as follows: `mpirun -np N ./v3 test BLOCKSIZE` where N is the number of MPI nodes. During our testing, even by setting manually the number of OMP Threads allocated to each MPI node, the system simply does not utilize them.

The code production was unsupervised and was done in colaboration with [AntoniosOurdas](https://github.com/AntoniosOurdas) .
