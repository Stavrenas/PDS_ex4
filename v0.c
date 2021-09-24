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
    char filenameB[] = "12a.mtx";

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    addMatrix(A, B, C);
    printMatrix(C);

    // blockMatrix(A, 4, blockA);
    // blockMatrix(B, 4, blockB);

    // multMatrix(A, B, C);
    // blockMatrix(C, 4, blockResult);
    // printf("====BMM Result====\n\n");
    // printBlockedMatrix(blockResult);

    // multBlockedMatrix(blockA, blockB, blockC);
    // printf("====BlockBMM Result====\n");
    //printBlockedMatrix(blockC);

    //saveMatrix(C, "mycielskianPARALLEL.txt");

    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);
    MPI_Finalize();

    // free memory

    // clearMatrix(A);
    // clearMatrix(B);
    // clearMatrix(C);
    // clearBlockedMatrix(blockA);
}

void multBlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C)
{
    uint32_t size, totalBlocks;

    uint32_t blockSize = A->list[0]->size;
    uint32_t maxBlocks = ceil(A->size / blockSize);

    size = 1;
    C->list = (Matrix **)malloc(size * sizeof(Matrix *));
    C->offsets = (uint32_t *)malloc(size * sizeof(uint32_t));

    for (uint32_t blockY = 1; blockY <= maxBlocks; blockY++)
    {
        for (uint32_t blockX = 1; blockX <= maxBlocks; blockX++)
        {
            printf("Row %d Col %d\n", blockY, blockX);
            //Create block: Cp,q (p = BlockY, q = BlockX)

            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            Matrix *result = (Matrix *)malloc(sizeof(Matrix)); //used for mult

            block->csc_elem = (uint32_t *)malloc((block->size + 1) * sizeof(uint32_t));
            block->csc_idx = (uint32_t *)malloc((0) * sizeof(uint32_t));
            block->size = A->list[0]->size; //get blocksize

            //initialize block
            for (int i = 0; i <= block->size; i++)
                block->csc_elem[i] = 0;

            uint32_t indexA, indexB; //find indexes of Ap1 and B1q

            for (int i = 1; i <= maxBlocks; i++)
            {
                indexA = findIndex(A, (blockY - 1) * maxBlocks + i);
                if (indexA != -1)
                    break; //stop when we find the first nonzero block in row p
            }
            //printf("After find index\n");
            for (int i = 1; i <= maxBlocks; i++)
            {
                indexB = findIndex(B, maxBlocks * (i - 1) + blockX);
                if (indexB != -1)
                    break; //stop when we find the first nonzero block in column q
            }

            printf("indexA is %d and indexB is %d\n", indexA, indexB);
            uint32_t blocksAdded = 0;

            for (int s = 1; s <= maxBlocks; s++)
            {
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

                printf("OffsetA is %d and OffsetB is %d\n", offsetA, offsetB);

                //check if the blocks match
                if (offsetA % maxBlocks == floor((offsetB - 1) / maxBlocks) + 1)
                {
                    // printf("s=%d\n", s);

                    printf("****Multipling A: %d with B: %d ****\n", offsetA, offsetB);
                    // printf("A\n");
                    // printMatrix(A->list[indexA]);
                    // printf("B\n");
                    // printMatrix(B->list[indexB]);

                    multMatrix2(A->list[indexA], B->list[indexB], result);

                    // printf("result\n");
                    // printMatrix(result);
                    // printf("Block\n");
                    // printMatrix(block);

                    addMatrix(result, block, block);

                    // printf("Block\n");
                    // printMatrix(block);

                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks);
                        if (indexB != -1)
                            break; //stop when we find the first nonzero block in column q
                    }
                    indexA++;
                    blocksAdded = 1;
                }
                else if (offsetA % maxBlocks > floor(offsetB / maxBlocks) + 1)
                {
                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks);
                        if (indexB != -1)
                            break; //stop when we find the first nonzero block in column q
                    }
                }

                else if (offsetA % maxBlocks < floor(offsetB / maxBlocks) + 1)
                    indexA++; //go to the next block in the same line
            }

            if (blocksAdded == 1)
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
            free(result);
        }
    }
    C->size = A->size;
    C->totalBlocks = totalBlocks;
}
