#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#include "utilities.h"
#include "omp_utilities.h"
#include "mpi_utilities.h"
#include "mmio.h"
#include "read.h"

void multBlockedMatrixMPI(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C)
{
    uint32_t size, totalBlocks, maxBlocks, blockSize;

    int rank, num_tasks;
    MPI_Status mpistat;
    MPI_Request mpireq; //initialize MPI environment
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    blockSize = A->list[0]->size;
    maxBlocks = floor(A->size / blockSize) + 1;
    if (A->size % blockSize == 0)
        maxBlocks--;
    size = 1;
    totalBlocks = 0;

    //initialize result matrix
    C->list = (Matrix **)malloc(size * sizeof(Matrix *));
    C->offsets = (uint32_t *)malloc(size * sizeof(uint32_t));
    C->row_ptr = (int *)calloc(maxBlocks, sizeof(int));

    for (uint32_t blockY = 1 + rank; blockY <= maxBlocks; blockY += num_tasks)
    {
        C->row_ptr[blockY - 1] = totalBlocks;

        for (uint32_t blockX = 1; blockX <= maxBlocks; blockX++)
        {
            //Create block: Cp,q (p = BlockY, q = BlockX)
            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            Matrix *result = (Matrix *)malloc(sizeof(Matrix)); //used for mult

            //initialize block
            block->size = blockSize;
            block->csc_elem = (uint32_t *)malloc((blockSize + 1) * sizeof(uint32_t));
            block->csc_idx = (uint32_t *)malloc((0) * sizeof(uint32_t));

            for (int i = 0; i <= blockSize; i++)
                block->csc_elem[i] = 0;

            //find indexes of Ap1 and B1q
            uint32_t indexA, indexB;
            for (int i = 1; i <= maxBlocks; i++)
            {
                indexA = findIndex(A, (blockY - 1) * maxBlocks + i);
                if (indexA != -1)
                    break; //stop when we find the first nonzero block in row p
            }

            for (int i = 1; i <= maxBlocks; i++)
            {
                indexB = findIndex(B, maxBlocks * (i - 1) + blockX);
                if (indexB != -1)
                    break; //stop when we find the first nonzero block in column q
            }

            //maxBlocks is the maximum number of mults for Cp,q. Variable s is not used
            for (int s = 1; s <= maxBlocks; s++)
            {

                //if either block does not exist
                if (indexA == -1 || indexB == -1)
                    break;

                uint32_t offsetA = A->offsets[indexA];
                uint32_t offsetB = B->offsets[indexB];

                //break if blocks are out of the desired row/col
                if (offsetA > blockY * maxBlocks || offsetB % maxBlocks > blockX)
                    break;
                //or if we run out of blocks
                if (indexB > B->totalBlocks || indexA > A->totalBlocks)
                    break;

                //check if the blocks match
                uint32_t sA = (offsetA - 1) % maxBlocks;
                uint32_t sB = floor((offsetB - 1) / maxBlocks);

                if (sA == sB)
                {
                    if (blockSize <= 40)
                        multMatrix2(A->list[indexA], B->list[indexB], result);
                    else
                        multMatrixParallel(A->list[indexA], B->list[indexB], result);

                    addMatrix(result, block, block);

                    //find block Bsq
                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks * i);
                        if (indexB != -1)
                            break;
                    }

                    indexA++; //go to the next block in the same line of A
                }

                else if (sA > sB)
                {
                    //find block Bsq
                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks * i);
                        if (indexB != -1)
                            break;
                    }
                }

                else if (sA < sB)
                    indexA++; //go to the next block in the same line of A
            }

            // if the mult results in a nonzero block, add it to the result matrix
            if (block->csc_elem[blockSize] != 0)
            {
                C->list[totalBlocks] = block;
                C->offsets[totalBlocks] = (blockY - 1) * maxBlocks + blockX;
                totalBlocks++;
                if (size == totalBlocks)
                {
                    size++;
                    C->list = realloc(C->list, size * sizeof(Matrix *));
                    C->offsets = realloc(C->offsets, size * sizeof(uint32_t *));
                }
            }
            free(result); //result matrix will not be needed in the future, counter to "block" matrix
        }
    }

    //fill the missing offsets
    for (int i = 0; i < totalBlocks; i++)
    {
        uint32_t row = (C->offsets[i] - 1) / maxBlocks;
        C->row_ptr[row + 1] = i + 1;
    }

    for (int i = 1; i < maxBlocks; i++)
    {
        if (C->row_ptr[i] == 0)
            C->row_ptr[i] = C->row_ptr[i - 1];
    }

    C->size = A->size;
    C->blockSize = A->blockSize;
    C->totalBlocks = totalBlocks;
}

void multBlockedMatrixMPIMasked(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C, BlockedMatrix *mask)
{
    uint32_t size, totalBlocks, maxBlocks, blockSize;

    int rank, num_tasks;
    MPI_Status mpistat;
    MPI_Request mpireq; //initialize MPI environment
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    blockSize = A->list[0]->size;
    maxBlocks = floor(A->size / blockSize) + 1;
    if (A->size % blockSize == 0)
        maxBlocks--;
    size = 1;
    totalBlocks = 0;

    //initialize result matrix
    C->list = (Matrix **)malloc(size * sizeof(Matrix *));
    C->offsets = (uint32_t *)malloc(size * sizeof(uint32_t));
    C->row_ptr = (int *)calloc(maxBlocks, sizeof(int));

    for (uint32_t index = 0; index < mask->totalBlocks; index++)
    {
        uint32_t offset = mask->offsets[index];
        uint32_t blockY = (offset - 1) / maxBlocks + 1;

        if ((blockY - 1) % num_tasks == rank) //chech if the row corresponds to the correct thread
        {
            uint32_t blockX = (offset - 1) % maxBlocks + 1;
            //Create block: Cp,q (p = BlockY, q = BlockX)
            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            Matrix *result = (Matrix *)malloc(sizeof(Matrix)); //used for mult

            //initialize block and result
            block->size = blockSize;
            block->csc_elem = (uint32_t *)calloc((blockSize + 1), sizeof(uint32_t));
            block->csc_idx = (uint32_t *)malloc((0) * sizeof(uint32_t));

            //find indexes of Ap1 and B1q
            uint32_t indexA, indexB;
            for (int i = 1; i <= maxBlocks; i++)
            {
                indexA = findIndex(A, (blockY - 1) * maxBlocks + i);
                if (indexA != -1)
                    break; //stop when we find the first nonzero block in row p
            }

            for (int i = 1; i <= maxBlocks; i++)
            {
                indexB = findIndex(B, maxBlocks * (i - 1) + blockX);
                if (indexB != -1)
                    break; //stop when we find the first nonzero block in column q
            }

            //maxBlocks is the maximum number of mults for Cp,q. Variable s is not used
            for (int s = 1; s <= maxBlocks; s++)
            {
                //if either block does not exist
                if (indexA == -1 || indexB == -1)
                    break;

                uint32_t offsetA = A->offsets[indexA];
                uint32_t offsetB = B->offsets[indexB];

                //break if blocks are out of the desired row/col
                if (offsetA > blockY * maxBlocks || offsetB % maxBlocks > blockX)
                    break;
                //or if we run out of blocks
                if (indexB > B->totalBlocks || indexA > A->totalBlocks)
                    break;

                //check if the blocks match
                uint32_t sA = (offsetA - 1) % maxBlocks;
                uint32_t sB = (offsetB - 1) / maxBlocks;

                if (sA == sB)
                {
                    //choose best algorithm for speeeeed
                    if (blockSize <= 40)
                        multMatrixMasked(A->list[indexA], B->list[indexB], result, mask->list[index]);
                    else
                        multMatrixParallelMasked(A->list[indexA], B->list[indexB], result, mask->list[index]);

                    addMatrix(result, block, block);

                    //find block Bsq
                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks * i);
                        if (indexB != -1)
                            break;
                    }

                    indexA++; //go to the next block in the same line of A
                }

                else if (sA > sB)
                {
                    //find block Bsq
                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks * i);
                        if (indexB != -1)
                            break;
                    }
                }

                else if (sA < sB)
                    indexA++; //go to the next block in the same line of A
            }

            // if the mult results in a nonzero block, add it to the result matrix
            if (block->csc_elem[blockSize] != 0)
            {
                C->list[totalBlocks] = block;
                C->offsets[totalBlocks] = offset;
                totalBlocks++;
                if (size == totalBlocks)
                {
                    size++;
                    C->list = realloc(C->list, size * sizeof(Matrix *));
                    C->offsets = realloc(C->offsets, size * sizeof(uint32_t *));
                }
            }

            C->row_ptr[blockY] = totalBlocks;
            free(result);
        }
    }

    for (int i = 0; i < totalBlocks; i++)
    {
        uint32_t row = (C->offsets[i] - 1) / maxBlocks;
        C->row_ptr[row + 1] = i + 1;
    }

    for (int i = 1; i < maxBlocks; i++)
    {
        if (C->row_ptr[i] == 0)
            C->row_ptr[i] = C->row_ptr[i - 1];
    }

    C->size = A->size;
    C->blockSize = A->blockSize;
    C->totalBlocks = totalBlocks;
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

    struct timeval start = tic();
    multBlockedMatrixMPI(A, B, Cp_blocked);
    //printf("MultBlocked time(%d) : %f\n", rank, toc(start));
    //start = tic();

    unblockMatrix(Cp_blocked, Cp);
    //printf("Unblock time(%d) : %f\n", rank, toc(start));

    // All proccesses except 0 send data to 0
    int tag = 99;
    if (rank != 0)
    {
        // Send csc_idx size to proccess 0
        MPI_Send(&Cp->csc_elem[Cp->size], 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
        MPI_Send(Cp->csc_idx, Cp->csc_elem[Cp->size], MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
        MPI_Send(Cp->csc_elem, Cp->size + 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
    }
    else
    {
        Matrix **C_recv = (Matrix **)malloc(num_tasks * sizeof(Matrix *));
        uint32_t *idx_size = (uint32_t *)malloc(num_tasks * sizeof(uint32_t));
        idx_size[0] = Cp->csc_elem[Cp->size];

        // Receive size of each matrix
        for (int i = 1; i < num_tasks; ++i)
            MPI_Recv(&idx_size[i], 1, MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);

        uint32_t elements = 0;
        for (int i = 1; i < num_tasks; i++)
        {
            // Check if matrix contains non-zero blocks
            // Allocate memory for new matrix

            C_recv[i - 1] = (Matrix *)malloc(sizeof(Matrix));
            C_recv[i - 1]->size = A->size;
            C_recv[i - 1]->csc_idx = (uint32_t *)malloc(idx_size[i] * sizeof(uint32_t));
            C_recv[i - 1]->csc_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
        }

        // Start receiving csc_elem data asychronously

        for (int i = 1; i < num_tasks; ++i)
            MPI_Recv(C_recv[i - 1]->csc_idx, idx_size[i], MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);

        // Start receiving csc_idx data asychronously
        for (int i = 1; i < num_tasks; i++)
            MPI_Recv(C_recv[i - 1]->csc_elem, C_recv[i - 1]->size + 1, MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);

        // Merge matrices by adding them
        //start = tic();
        Matrix *temp = (Matrix *)malloc(sizeof(Matrix));
        for (int i = 0; i < num_tasks - 1; i++)
        {
            addMatrix(C_recv[i], Cp, temp);
            Cp = temp;
        }
        //printf("Add time : %f\n", toc(start));
        return Cp;
    }
    return Cp;
}

Matrix *MPI_MultMasked(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *mask)
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

    struct timeval start = tic();
    multBlockedMatrixMPIMasked(A, B, Cp_blocked, mask);
    //printf("MultBlocked time(%d) : %f\n", rank, toc(start));
    //start = tic();

    unblockMatrix(Cp_blocked, Cp);
    //printf("Unblock time(%d) : %f\n", rank, toc(start));

    // All proccesses except 0 send data to 0
    int tag = 99;
    if (rank != 0)
    {
        // Send csc_idx size to proccess 0
        MPI_Send(&Cp->csc_elem[Cp->size], 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
        MPI_Send(Cp->csc_idx, Cp->csc_elem[Cp->size], MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
        MPI_Send(Cp->csc_elem, Cp->size + 1, MPI_UINT32_T, 0, tag, MPI_COMM_WORLD);
    }
    else
    {
        Matrix **C_recv = (Matrix **)malloc(num_tasks * sizeof(Matrix *));
        uint32_t *idx_size = (uint32_t *)malloc(num_tasks * sizeof(uint32_t));
        idx_size[0] = Cp->csc_elem[Cp->size];

        // Receive size of each matrix
        for (int i = 1; i < num_tasks; ++i)
            MPI_Recv(&idx_size[i], 1, MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);

        uint32_t elements = 0;
        for (int i = 1; i < num_tasks; i++)
        {
            // Check if matrix contains non-zero blocks
            // Allocate memory for new matrix

            C_recv[i - 1] = (Matrix *)malloc(sizeof(Matrix));
            C_recv[i - 1]->size = A->size;
            C_recv[i - 1]->csc_idx = (uint32_t *)malloc(idx_size[i] * sizeof(uint32_t));
            C_recv[i - 1]->csc_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
        }

        // Start receiving csc_elem data asychronously

        for (int i = 1; i < num_tasks; ++i)
            MPI_Recv(C_recv[i - 1]->csc_idx, idx_size[i], MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);

        // Start receiving csc_idx data asychronously
        for (int i = 1; i < num_tasks; i++)
            MPI_Recv(C_recv[i - 1]->csc_elem, C_recv[i - 1]->size + 1, MPI_UINT32_T, i, tag, MPI_COMM_WORLD, &mpistat);

        // Merge matrices by adding them
        //start = tic();
        Matrix *temp = (Matrix *)malloc(sizeof(Matrix));
        for (int i = 0; i < num_tasks - 1; i++)
        {
            addMatrix(C_recv[i], Cp, temp);
            Cp = temp;
        }
        //printf("Add time : %f\n", toc(start));
        return Cp;
    }
    return Cp;
}