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

    //free memory
}

void blockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix)
{
    uint32_t maxBlocks = ceil(mtr->size / blockSize);
    uint32_t totalBblocks = 0;
    uint32_t listSize = 1;  //also equals to offset size

    blockedMatrix->list = (Matrix **)malloc(1 * sizeof(Matrix *));
    blockedMatrix->offsets = (uint32_t *)malloc(1 * sizeof(uint32_t)); //maximum size of blocks

    for (uint32_t blockY = 1; blockY <= maxBlocks; blockY++)
    {
        for (uint32_t blockX = 1; blockX <= maxBlocks; blockX++)
        {

            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            int *block_idx, *block_elem, elements, idx_size;

            block_idx = (int *)malloc(sizeof(int));
            block_elem = (int *)malloc((blockSize + 1) * sizeof(int));

            idx_size = 1;
            block_elem[0] = 0;
            elements = 0;

            for (int row = (blockY - 1) * blockSize + 1; row < blockY * blockSize + 1; row++)
            {
                // Check if row exceeds matrix size
                if (row <= mtr->size)
                {
                    uint32_t start = mtr->csc_elem[row - 1];
                    uint32_t end = mtr->csc_elem[row];

                    for (int j = start; j < end; j++)
                    {

                        if (mtr->csc_idx[j] > blockX * blockSize) //check if it is worth it
                            break;

                        else if (mtr->csc_idx[j] <= blockX * blockSize && mtr->csc_idx[j] > (blockX - 1) * blockSize)
                        {

                            block_idx[elements] = mtr->csc_idx[j] - (blockX - 1) * blockSize;
                            elements++;

                            if (elements == idx_size)
                            {
                                idx_size ++;
                                block_idx = realloc(block_idx, idx_size * sizeof(int));
                            }
                        }
                    }

                    block_elem[row - (blockY - 1) * blockSize] = elements;
                }

                else // zero padding vertically
                    block_elem[row - (blockY - 1) * blockSize] = block_elem[row - (blockY - 1) * blockSize - 1];
            }

            if (elements != 0)
            {
                block->size = blockSize;
                block->csc_idx = block_idx;
                block->csc_elem = block_elem;
                blockedMatrix->list[totalBblocks] = block;

                blockedMatrix->offsets[totalBblocks] = (blockY - 1) * maxBlocks + blockX;
                totalBblocks++;

                if (listSize == totalBblocks)
                {
                    listSize ++;
                    blockedMatrix->list = realloc(blockedMatrix->list, listSize * sizeof(Matrix *));
                    blockedMatrix->offsets = realloc(blockedMatrix->offsets, listSize * sizeof(uint32_t *));

                }
            }
        }
    }

    blockedMatrix->size = totalBblocks;

    printf("Max blocks are %d and current blocks: %d. Non zero blocks: %f \n",maxBlocks * maxBlocks,totalBblocks, (float)(totalBblocks)/(maxBlocks * maxBlocks));
}
