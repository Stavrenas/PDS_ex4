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

    //printf("Rows size is %d\n", rows_size);
    for (uint32_t blockY = 1+rank; blockY <= maxBlocks; blockY+=num_tasks)
    {
        C->row_ptr[blockY - 1] = totalBlocks;

        //printf("blockY is %d\n", blockY);
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
                //printf("sA %d, sB %d, offsetA %d, offsetB %d\n",sA, sB, offsetA,offsetB);

                if (sA == sB)
                {
                    //printf("Mult %d with %d\n",offsetA,offsetB);

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
                //printf("Found offset %d\n",(blockY - 1) * maxBlocks + blockX);

                if (size == totalBlocks)
                {
                    size++;
                    C->list = realloc(C->list, size * sizeof(Matrix *));
                    C->offsets = realloc(C->offsets, size * sizeof(uint32_t *));
                }
            }
            //free(result); //result matrix will not be needed in the future, counter to "block" matrix
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