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

void MPI_Mult(BlockedMatrix *A, BlockedMatrix *B, Matrix *C);

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank, name_len;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);

    Matrix *A = (Matrix*)malloc(sizeof(Matrix));
    Matrix *B = (Matrix*)malloc(sizeof(Matrix));
    Matrix *C = NULL;
    if (world_rank == 0)
    {
        C = (Matrix*)malloc(sizeof(Matrix));
    }

    BlockedMatrix *blockA = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockB = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockC = NULL;
    if (world_rank == 0)
    {
        blockC = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));;
    }
    // BlockedMatrix *blockResult = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    char matrix[] = "12";
    char *filenameA = (char *)malloc(25 * sizeof(char));
    char *filenameB = (char *)malloc(25 * sizeof(char));
    char *name = (char *)malloc(25 * sizeof(char));
    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);
    // struct timeval start = tic();
    //
    //multMatrix2(A, B, C);
    //sprintf(name, "%s_serial.txt", matrix);
    // // printf("Time for serial mult: %f\n", toc(start));
    //saveMatrix(C, name);
    // // printMatrix(C);
    //
    // start = tic();
    // multMatrixParallel(A, B, C);
    // printf("Time for parallel mult: %f\n", toc(start));
    //
    // sprintf(name, "%s_parallel.txt", matrix);
    // saveMatrix(C, name);
    // //printMatrix(C);
    //
    // printf("\n");
    //
    // start = tic();
    int blocksize = 4;
    //
    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);
    //
    //multBlockedMatrix(blockA, blockB, blockC);
    //printBlockedMatrix(blockC);
    //
    //unblockMatrix(blockC, C);
    // printf("Time for blocked mult: %f\n", toc(start));
    // //printMatrix(C);
    //
    //sprintf(name, "%s_blocked.txt", matrix);
    //saveMatrix(C, name);
    //
    MPI_Mult(blockA, blockB, C);
    // printf("blockResult = %p\n", blockResult->list);
    // printMatrix(C);
    // sprintf(name, "%s_blockedMPI.txt", matrix);
    // saveMatrix(C, name);

    // // free memory
    // clearMatrix(A);
    // clearMatrix(B);
    // if (world_rank == 0)
    // {
    //     clearMatrix(C);
    // }
    // clearBlockedMatrix(blockA);
    // clearBlockedMatrix(blockB);
    // if (world_rank == 0)
    // {
    //     clearBlockedMatrix(blockC);
    // }
    // printf("Done here\n");

    MPI_Finalize();
}

void MPI_Mult(BlockedMatrix *A, BlockedMatrix *B, Matrix *C)
{
    // Given n MPI proccesses
    // Each proccess will calculate (A->size / n) rows
    // The remaining m rows (if any) are assigned to the first m processes

    // Get rank and number of tasks
    int rank = 0, num_tasks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    // Define result matrix and blocked matrix
    Matrix *Cp = (Matrix*)malloc(sizeof(Matrix));
    BlockedMatrix *Cp_blocked = malloc(sizeof(BlockedMatrix));

    int maxBlocks = floor(A->size / A->blockSize) + 1;
    if((A->size) % (A->blockSize) == 0)
        maxBlocks--;
    

    if(rank == 0)
    {
        printf("Size = %d\n", A->size);
        printf("Block size = %d\n", A->blockSize);
        printf("Max blocks = %d\n", maxBlocks);
    }

    // Number of rows per proccess
    int size = floor(maxBlocks / num_tasks);
    // Remaining rows
    int rest = (maxBlocks) % num_tasks;
    // Rows' size
    int rows_size = size;
    // Current proccess rows
    uint32_t *rows = (uint32_t *)malloc(rows_size * sizeof(uint32_t));

    // Assign rows to proccesses
    for (int i = 0; i < size; ++i)
        rows[i] = rank + i * num_tasks;

    // Assign remaining rows on proccesses 0 to (rest - 1)
    if (rank < rest)
    {
        rows = realloc(rows, (size + 1) * sizeof(uint32_t));
        rows[size] = size * num_tasks + rank;
        // Increase size
        rows_size = size + 1;
    }
    printf("Thread %d: ", rank);
    for (int i = 0; i < rows_size; ++i)
        printf("%d ", rows[i]);

    printf("\n");

    // Perform multiplication
    multBlockedMatrixMPI(A, B, Cp_blocked, rows, rows_size);
    printf("Mult %d\n", rank);
    printBlockedMatrix(Cp_blocked);
    // Unblock matrix Cp
    printf("++Before unblock %d\n", rank);

    unblockMatrix(Cp_blocked, Cp);

    printf("--Unblocked %d\n", rank);

    MPI_Status status;
    MPI_Request request[num_tasks - 1];

    // All proccesses except 0 send data to 0
    int tag = 99;
    if (rank != 0)
    {
        // Send csc_idx size to proccess 0
        MPI_Send(&Cp->csc_elem[Cp->size], 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
        if (Cp->csc_elem[Cp->size] != 0)
        {
            MPI_Send(Cp->csc_idx, Cp->csc_elem[Cp->size], MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
            MPI_Send(Cp->csc_elem, Cp->size, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
        }
    }
    else
    {
        // Set list of matrices to be received from other proccesses
        Matrix **C_recv = (Matrix **)malloc(0 * sizeof(Matrix *));

        uint32_t elements = 0;

        // Size of each proccess' matrix
        uint32_t *idx_size = (uint32_t *)malloc(num_tasks * sizeof(uint32_t));
        idx_size[0] = Cp->csc_elem[Cp->size];

        // Receive size of each matrix
        for (int i = 1; i < num_tasks; ++i)
            MPI_Irecv(&idx_size[i], 1, MPI_UINT32_T, i, 99, MPI_COMM_WORLD, &request[i - 1]);

        // Wait to receive all sizes
        for (int i = 1; i < num_tasks; i++)
            MPI_Wait(&request[i - 1], &status);

        // Allocate memory for non-empty matrices to be received
        C_recv = malloc(num_tasks * sizeof(Matrix *));  //better malloc than realloc
        for (int i = 1; i < num_tasks; ++i)
        {
            printf("rank: %d, size: %d\n", rank, idx_size[i]);
            // Check if matrix contains non-zero blocks
            if (idx_size[i] != 0)
            {
                // Allocate memory for new matrix
                elements++;
                C_recv[elements - 1] = (Matrix *)malloc(1 * sizeof(Matrix));
                C_recv[elements - 1]->size = A->size;
                C_recv[elements - 1]->csc_idx = (uint32_t *)malloc(idx_size[i] * sizeof(uint32_t));
                C_recv[elements - 1]->csc_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
            }
        }

        // Start receiving csc_elem data asychronously
        elements = 0;
        for (int i = 1; i < num_tasks; ++i)
        {
            if (idx_size[i] != 0)
            {
                printf("receiving %d idx_elements\n", idx_size[i]);
                elements++;
                MPI_Irecv(C_recv[elements - 1]->csc_idx, idx_size[i], MPI_UINT32_T, i, 99, MPI_COMM_WORLD, &request[i - 1]);
            }
        }

        // Wait to receive all csc_idx data
        for (int i = 1; i < num_tasks; ++i)

            MPI_Wait(&request[i - 1], &status);

        // Start receiving csc_idx data asychronously
        elements = 0;
        for (int i = 1; i < num_tasks; ++i)
        {
            if (idx_size[i] != 0)
            {
                printf("receiving %d elem_elements\n", idx_size[i]);
                elements++;
                MPI_Irecv(C_recv[elements - 1]->csc_elem, C_recv[elements - 1]->size, MPI_UINT32_T, i, 99, MPI_COMM_WORLD, &request[i - 1]);
            }
        }

        // Wait to receive all csc_elem data
        for (int i = 1; i < num_tasks; ++i)
            MPI_Wait(&request[i - 1], &status);

        // Merge matrices by adding them
        printf("Thread %d: adding\n", rank);
        for (int i = 0; i < elements; ++i)
            addMatrix(C_recv[i], Cp, Cp);

        // return;
    }
    clearMatrix(Cp);
    clearBlockedMatrix(Cp_blocked);
}
