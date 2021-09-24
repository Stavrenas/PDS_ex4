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

    char filenameA[] = "12.mtx";
    char filenameB[] = "12.mtx";

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    // addMatrix(A, B, C);
    // printMatrix(C);

    blockMatrix(A, 4, blockA);
    blockMatrix(B, 4, blockB);

    multMatrix(A, B, C);
    blockMatrix(C, 4, blockResult);
    printf("====BMM Result====\n\n");
    printBlockedMatrix(blockResult);

    multBlockedMatrix(blockA, blockB, blockC);
    printf("====BlockBMM Result====\n");
    printBlockedMatrix(blockC);

    //saveMatrix(C, "mycielskianPARALLEL.txt");

    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);
    MPI_Finalize();

    // free memory

    // clearMatrix(A);
    // clearMatrix(B);
    // clearMatrix(C);
    // clearBlockedMatrix(blockA);
}


