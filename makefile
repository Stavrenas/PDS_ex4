CC=gcc
MPICC=mpicc
CFLAGS=-O3 -lm
MP=-fopenmp

default: all

all: v0 v1 v2 v3 

v0: v0.c utilities.c read.c mmio.c
	$(CC) $(CFLAGS) -o $@ $^
v1: v1.c utilities.c omp_utilities.c read.c mmio.c
	$(CC) $(CFLAGS) $(MP) -o $@ $^
v2: v2.c utilities.c omp_utilities.c read.c mmio.c
	$(CC) $(CFLAGS) $(MP) -o $@ $^
v3: v3.c utilities.c omp_utilities.c mpi_utilities.c read.c mmio.c
	$(MPICC) $(CFLAGS) $(MP) -lrt -o $@ $^ 

clean:
	rm -f v0 v1 v2 v3 *.txt
