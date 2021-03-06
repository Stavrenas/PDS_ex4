#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "utilities.h"
#include "omp_utilities.h"
#include "mpi_utilities.h"
#include "read.h"
#include "mmio.h"

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    Matrix *A = (Matrix *)malloc(sizeof(Matrix));
    Matrix *B = (Matrix *)malloc(sizeof(Matrix));
    Matrix *C = (Matrix *)malloc(sizeof(Matrix));

    BlockedMatrix *blockA = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockB = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *mask = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    // BlockedMatrix *blockResult = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    char *matrix = (char *)malloc(40 * sizeof(char));
    int blocksize = 1;

    if (argc == 1)
    {
        sprintf(matrix, "%s", "mycielskian");
        blocksize = 250;
    }
    else if (argc == 3)
    {
        sprintf(matrix, "%s", argv[1]);
        blocksize = atoi(argv[2]);
    }

    if (world_rank == 0)
    {
        if (argc != 1 && argc != 3)
            printf("Usage: ./v3 matrix_name blocksize %d\n", argc);
    }

    char *filenameA = (char *)malloc(40 * sizeof(char));
    char *filenameB = (char *)malloc(40 * sizeof(char));
    char *name = (char *)malloc(40 * sizeof(char));
    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    struct timeval start = tic();

    //temp, only for sustoixia
    blocksize = A->size / world_size;
    //temp, only for sustoixia

    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);

    MPI_Barrier(MPI_COMM_WORLD); //sync MPI threads

    mask = blockA;
    C = MPI_MultMasked(blockA, blockB, mask);

    if (world_rank == 0)
    {
        sprintf(name, "%s_blockedMPI_%d.txt", matrix, world_size);
        saveMatrix(C, name);
        printf("MPI mult time: %f\n", toc(start));
    }

    MPI_Finalize();
}
