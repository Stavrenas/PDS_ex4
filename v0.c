#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "utilities.h"
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
    Matrix *C = malloc(sizeof(Matrix));
    BlockedMatrix *blockA = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockB = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockC = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockResult = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    char matrix[] = "50";
    char *filenameA = (char *)malloc(25 * sizeof(char));
    char *filenameB = (char *)malloc(25 * sizeof(char));
    char *name = (char *)malloc(25 * sizeof(char));
    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);
    struct timeval start = tic();

    // multMatrix2(A, B, C);
    // sprintf(name, "%s_serial.txt", matrix);
    // printf("Time for serial mult: %f\n", toc(start));
    // saveMatrix(C, name);
    // printMatrix(C);

    // start = tic();
    // multMatrixParallel(A, B, C);
    // sprintf(name, "%s_parallel.txt", matrix);
    // printf("Time for parallel mult: %f\n", toc(start));
    // //saveMatrix(C, name);
    // printMatrix(C);

    printf("\n");

    start = tic();
    int blocksize = 5;

     blockMatrix(A, blocksize, blockA);
     printMatrix(A);
     unblockMatrix(blockA, C);
    // blockMatrix(B, blocksize, blockB);
    // multBlockedMatrix(blockA, blockB, blockC);
    // unblockMatrix(blockC, C);
    // sprintf(name, "%s_blocked.txt", matrix);
    // printf("Time for blocked mult: %f\n", toc(start));
     printMatrix(C);
    // saveMatrix(C, name);

    MPI_Finalize();

    // free memory

    // clearMatrix(A);
    // clearMatrix(B);
    // clearMatrix(C);
    // clearBlockedMatrix(blockA);
}
