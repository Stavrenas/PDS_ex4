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

Matrix *MPI_Mult(BlockedMatrix *A, BlockedMatrix *B);

Matrix *MPI_Mult2(BlockedMatrix *A, BlockedMatrix *B);

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank, name_len;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);

    Matrix *A = (Matrix *)malloc(sizeof(Matrix));
    Matrix *B = (Matrix *)malloc(sizeof(Matrix));
    Matrix *C = NULL;
    if (world_rank == 0)
        C = (Matrix *)malloc(sizeof(Matrix));

    BlockedMatrix *blockA = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockB = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockC = NULL;
    if (world_rank == 0)
        blockC = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    // BlockedMatrix *blockResult = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    char matrix[] = "mycielskian";
    char *filenameA = (char *)malloc(40 * sizeof(char));
    char *filenameB = (char *)malloc(40 * sizeof(char));
    char *name = (char *)malloc(40 * sizeof(char));
    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    int blocksize = 4;

    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);

    C = MPI_Mult(blockA, blockB);
    printf("exited (%d)\n", world_rank);

    if (world_rank == 0)
    {
        sprintf(name, "%s_blockedMPI_%d.txt", matrix, world_size);
        saveMatrix(C, name);
    }

    // printf("Done here\n");

    MPI_Finalize();
}

Matrix *MPI_Mult2(BlockedMatrix *A, BlockedMatrix *B)
{
    // Given n MPI proccesses
    // Each proccess will calculate (A->size / n) rows
    // The remaining m rows (if any) are assigned to the first m processes

    // Get rank and number of tasks
    int rank, num_tasks;
    MPI_Status mpistat;
    MPI_Request mpireq; //initialize MPI environment
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    // Define result matrix and blocked matrix
    Matrix *Cp = (Matrix *)malloc(sizeof(Matrix));
    BlockedMatrix *Cp_blocked = malloc(sizeof(BlockedMatrix));

    int maxBlocks = floor(A->size / A->blockSize) + 1;
    if ((A->size) % (A->blockSize) == 0)
        maxBlocks--;

    if (rank == 0)
    {
        printf("Size = %d\n", A->size);
        printf("Block size = %d\n", A->blockSize);
        printf("Max blocks = %d\n", maxBlocks);
    }

    // Number of rows per proccess
    int size = floor(maxBlocks / num_tasks);
    // Remaining rows
    int rest = (maxBlocks) % num_tasks;
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
        rows_size = size + 1;
    }

    // printf("Thread %d: ", rank);
    // for (int i = 0; i < rows_size; ++i)
    //     printf("%d ", rows[i]);
    // printf("\n");

    multBlockedMatrixMPI(A, B, Cp_blocked, rows, rows_size);
    // printf("Mult %d\n", rank);
    // printBlockedMatrix(Cp_blocked);

    unblockMatrix(Cp_blocked, Cp);

    // printf("Original Cp for thread %d is\n", rank);
    // printMatrix(Cp);

    // All proccesses except 0 send data to 0
    int tag = 99;
    if (rank != 0)
    {
        // Send csc_idx size to proccess 0
        MPI_Send(&Cp->csc_elem[Cp->size], 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
        MPI_Send(Cp->csc_idx, Cp->csc_elem[Cp->size], MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);

        printf("Sent %d idx elements(%d)\n", Cp->csc_elem[Cp->size], rank);

        MPI_Send(Cp->csc_elem, Cp->size + 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);

        printf("Sent %d elem elements (%d)\n", Cp->size + 1, rank);
    }
    else
    {
        Matrix **C_recv = (Matrix **)malloc(num_tasks * sizeof(Matrix *));
        uint32_t *idx_size = (uint32_t *)malloc(num_tasks * sizeof(uint32_t));
        idx_size[0] = Cp->csc_elem[Cp->size];

        // Receive size of each matrix
        for (int i = 1; i < num_tasks; ++i)
        {
            MPI_Recv(&idx_size[i], 1, MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);
            printf("received: idx_size = %d (%d)\n", idx_size[i], i);
        }

        // // Wait to receive all sizes
        // for (int i = 1; i < num_tasks; i++)
        //     MPI_Wait(&request[i - 1], &status);

        uint32_t elements = 0;
        for (int i = 1; i < num_tasks; i++)
        {
            //printf("rank: %d, size: %d\n", rank, idx_size[i]);
            // Check if matrix contains non-zero blocks
            // Allocate memory for new matrix

            C_recv[i - 1] = (Matrix *)malloc(sizeof(Matrix));
            C_recv[i - 1]->size = A->size;
            C_recv[i - 1]->csc_idx = (uint32_t *)malloc(idx_size[i] * sizeof(uint32_t));
            C_recv[i - 1]->csc_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
        }

        // Start receiving csc_elem data asychronously

        for (int i = 1; i < num_tasks; ++i)
        {
            MPI_Recv(C_recv[i - 1]->csc_idx, idx_size[i], MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);
            printf("Received %d idx elements from %d\n", idx_size[i], i);
        }
        // // Wait to receive all csc_idx data
        // for (int i = 1; i < num_tasks; ++i)
        //     MPI_Wait(&request[i - 1], &status);

        // Start receiving csc_idx data asychronously
        for (int i = 1; i < num_tasks; i++)
        {

            MPI_Recv(C_recv[i - 1]->csc_elem, C_recv[i - 1]->size + 1, MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);
            printf("received %d elem_elements from %d\n", C_recv[i - 1]->size + 1, i);
        }

        // Merge matrices by adding them
        printf("Thread %d: adding\n", rank);

        Matrix *temp = (Matrix *)malloc(sizeof(Matrix));
        for (int i = 0; i < num_tasks - 1; i++)
        {
            // printMatrix(C_recv[i]);
            // printMatrix(Cp);
            addMatrix(C_recv[i], Cp, temp);
            Cp = temp;
        }
        printf("Thread %d: added\n", rank);
        //return Cp;
    }
}

Matrix *MPI_Mult(BlockedMatrix *A, BlockedMatrix *B)
{
    int rank, num_tasks;
    MPI_Status mpistat;
    MPI_Request mpireq; //initialize MPI environment
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    // Define result matrix and blocked matrix
    Matrix *Cp = (Matrix *)malloc(sizeof(Matrix));
    BlockedMatrix *Cp_blocked = malloc(sizeof(BlockedMatrix));

    int maxBlocks = floor(A->size / A->blockSize) + 1;
    if ((A->size) % (A->blockSize) == 0)
        maxBlocks--;

    // if (rank == 0)
    // {
    //     printf("Size = %d\n", A->size);
    //     printf("Block size = %d\n", A->blockSize);
    //     printf("Max blocks = %d\n", maxBlocks);
    // }

    // printf("Thread %d: ", rank);
    // for (int i = 0; i < rows_size; ++i)
    //     printf("%d ", rows[i]);
    // printf("\n");

    multBlockedMatrixMPI(A, B, Cp_blocked);
    // printf("Mult %d\n", rank);
    // printBlockedMatrix(Cp_blocked);

    unblockMatrix(Cp_blocked, Cp);

    // printf("Original Cp for thread %d is\n", rank);
    // printMatrix(Cp);

    // All proccesses except 0 send data to 0
    int tag = 99;
    if (rank != 0)
    {
        // Send csc_idx size to proccess 0
        MPI_Send(&Cp->csc_elem[Cp->size], 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
        MPI_Send(Cp->csc_idx, Cp->csc_elem[Cp->size], MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);

        printf("Sent %d idx elements(%d)\n", Cp->csc_elem[Cp->size], rank);

        MPI_Send(Cp->csc_elem, Cp->size + 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);

        printf("Sent %d elem elements (%d)\n", Cp->size + 1, rank);
    }
    else
    {
        Matrix **C_recv = (Matrix **)malloc(num_tasks * sizeof(Matrix *));
        uint32_t *idx_size = (uint32_t *)malloc(num_tasks * sizeof(uint32_t));
        idx_size[0] = Cp->csc_elem[Cp->size];

        // Receive size of each matrix
        for (int i = 1; i < num_tasks; ++i)
        {
            MPI_Recv(&idx_size[i], 1, MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);
            printf("received: idx_size = %d (%d)\n", idx_size[i], i);
        }

        // // Wait to receive all sizes
        // for (int i = 1; i < num_tasks; i++)
        //     MPI_Wait(&request[i - 1], &status);

        uint32_t elements = 0;
        for (int i = 1; i < num_tasks; i++)
        {
            //printf("rank: %d, size: %d\n", rank, idx_size[i]);
            // Check if matrix contains non-zero blocks
            // Allocate memory for new matrix

            C_recv[i - 1] = (Matrix *)malloc(sizeof(Matrix));
            C_recv[i - 1]->size = A->size;
            C_recv[i - 1]->csc_idx = (uint32_t *)malloc(idx_size[i] * sizeof(uint32_t));
            C_recv[i - 1]->csc_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
        }

        // Start receiving csc_elem data asychronously

        for (int i = 1; i < num_tasks; ++i)
        {
            MPI_Recv(C_recv[i - 1]->csc_idx, idx_size[i], MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);
            printf("Received %d idx elements from %d\n", idx_size[i], i);
        }
        // // Wait to receive all csc_idx data
        // for (int i = 1; i < num_tasks; ++i)
        //     MPI_Wait(&request[i - 1], &status);

        // Start receiving csc_idx data asychronously
        for (int i = 1; i < num_tasks; i++)
        {

            MPI_Recv(C_recv[i - 1]->csc_elem, C_recv[i - 1]->size + 1, MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);
            printf("received %d elem_elements from %d\n", C_recv[i - 1]->size + 1, i);
        }

        // Merge matrices by adding them
        printf("Thread %d: adding\n", rank);

        Matrix *temp = (Matrix *)malloc(sizeof(Matrix));
        for (int i = 0; i < num_tasks - 1; i++)
        {
            //  printMatrix(C_recv[i]);
            //  printMatrix(Cp);
            addMatrix(C_recv[i], Cp, temp);
            Cp = temp;
        }
        printf("Thread %d: added\n", rank);
        return Cp;
    }
    return Cp;
}