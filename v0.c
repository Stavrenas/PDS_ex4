#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h> // sqrt
#include <mpi.h>
#include <string.h>
#include <omp.h>
#include "utilities.h"
#include "controller.h"
#include "read.h"
#include "mmio.h"

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank, name_len;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);

    Matrix *A = malloc(sizeof(Matrix));
    Matrix *B = malloc(sizeof(Matrix));
    Matrix *res = malloc(sizeof(Matrix));
    BlockedMatrix *blockA = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    char filename[] = "dblp.mtx";
    readMatrix(filename, A);
    readMatrix(filename, B);

    struct timeval start = tic();
    blockMatrix(A, 100, blockA);
    printf("Time for block: %f\n", toc(start));

    //printBlockedMatrix(blockA);

    //printMatrix(A);

    // struct timeval start = tic();

    // cscBMM2(A, B, res);

    // printf("Time for mult: %f\n", toc(start));

    // printMatrix(res);
    //saveMatrix(res, "mycielskianPARALLEL.txt");

    // Blocking algorithms: BCSR or CSB

    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);
    MPI_Finalize();

    // free memory

    clearMatrix(A);
    clearMatrix(B);
    clearBlockedMatrix(blockA);
}


