CC=gcc
MPICC=mpicc
CFLAGS=-O3 -lm -g
MP=-fopenmp
N=2

default: all

all: v0 v1 v2 v3 

v0: v0.c utilities.c read.c mmio.c
	$(CC) $(CFLAGS) -o $@ $^
v1: v1.c utilities.c omp_utilities.c read.c mmio.c
	$(CC) $(CFLAGS) $(MP) -o $@ $^
v2: v2.c utilities.c omp_utilities.c read.c mmio.c
	$(CC) $(CFLAGS) $(MP) -o $@ $^
v3: v3.c utilities.c omp_utilities.c mpi_utilities.c read.c mmio.c
	$(MPICC) $(CFLAGS) $(MP) -o $@ $^

test: test.c utilities.c omp_utilities.c read.c mmio.c
	$(MPICC) $(CFLAGS) $(MP) -o $@ $^
run:
	mpirun -np $(N) ./v3

clean:
	rm -f v0 v1 v2 v3 test *.txt
