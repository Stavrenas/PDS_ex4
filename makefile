CC=gcc
MPICC=mpicc
CFLAGS=-O3 -lpthread -lm
MP=-fopenmp 
N=2

default: all

all: v0

v0: v0.c utilities.c read.c controller.c mmio.c
	$(MPICC) $(CFLAGS) -o $@ $^
test: test.c utilities.c read.c controller.c mmio.c
	$(MPICC) $(CFLAGS) -o $@ $^
run:
	mpirun -np $(np) ./v0

clean:
	rm -f v0 v1 v2



