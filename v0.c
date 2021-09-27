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

    char matrix[] = "mycielskian";
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

    start = tic();
    multMatrixParallel(A, B, C);
    printf("Time for parallel mult: %f\n", toc(start));

    sprintf(name, "%s_parallel.txt", matrix);
    saveMatrix(C, name);
    //printMatrix(C);

    printf("\n");

    start = tic();
    int blocksize = 50;

    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);

    multBlockedMatrix(blockA, blockB, blockC);

    unblockMatrix(blockC, C);
    printf("Time for blocked mult: %f\n", toc(start));
    //printMatrix(C);

    sprintf(name, "%s_blocked.txt", matrix);
    saveMatrix(C, name);

    MPI_Finalize();

    // free memory

    // clearMatrix(A);
    // clearMatrix(B);
    // clearMatrix(C);
    // clearBlockedMatrix(blockA);
}

void MPI_Mult(BlockedMatrix* A, BlockedMatrix *B, BlockedMatrix *C)
{
    // Given n MPI proccesses
    // Each proccess must calculate (maxBlocks*maxBlocks)/n blocks
    // using multMatrixParallel(A, B, C)
    // and then send data to proccess 0 to merge into final result
    // Example n = maxBlocks then each proccess p will calculate
    // blocks C(p,1) to C(p,maxBlocks)
    // Then send these blocks to proccess 0 using MPI_Send

    // Get number of tasks and rank
    int rank = 0, num_tasks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    // Proccess matrix chunk of blocks
    Matrix *Cp = malloc(sizeof(Matrix));

    // Matrix multiplication should have limits (p1, p2, q1, q2)
    // Here we say that multBlockedMatrix should calulate only one 'block-row'
    int p1 = rank + 1;
    int p2 = rank + 1;
    int q1 = rank*(A->size) + 1;
    int q2 = rank*(A->size) + A->maxBlocks;

    multBlockedMatrix(A, B, Cp, p1, p2, q1, q2);

    // All proccesses but 0 send data to 0
    if(rank != 0)
    {
        MPI_Recv();
        // Merge matrices Cp into C
    } else
    {
        MPI_Send();
        // maybe some memory deallocation
        return;
    }
}
