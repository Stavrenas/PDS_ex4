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
    BlockedMatrix *blockC = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    char filenameA[] = "12.mtx";
    char filenameB[] = "12a.mtx";

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    struct timeval start = tic();

    blockMatrix(A, 3, blockA);

    blockBMM(blockA, blockA, blockC);

    printf("Time for BlockBMM: %f\n", toc(start));

    printBlockedMatrix(blockC);

    //printMatrix(C);

    // cscBMM2(A, B, C);

    // printf("Time for mult: %f\n", toc(start));

    // printMatrix(C);
    //saveMatrix(C, "mycielskianPARALLEL.txt");

    // Blocking algorithms: BCSR or CSB

    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);
    MPI_Finalize();

    // free memory

    clearMatrix(A);
    clearMatrix(B);
    clearMatrix(C);
    // clearBlockedMatrix(blockA);
}

void blockBMM(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C)
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
            uint32_t offset = (blockY - 1) * maxBlocks + blockX;

            //Create block: Cp,q (p = BlockY, q = BlockX)

            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            Matrix *result = (Matrix *)malloc(sizeof(Matrix));

            uint32_t indexA, indexB; //find indexes of Ap1 and B1q

            for (int i = 0; i < maxBlocks; i++)
            {
                indexA = findIndex(A, (blockY - 1) * maxBlocks + i);
                if (indexA != -1)
                    break; //stop when we find the first nonzero block in row p
            }

            for (int i = 0; i < maxBlocks; i++)
            {
                indexB = findIndex(B, maxBlocks * i + blockX);
                if (indexB != -1)
                    break; //stop when we find the first nonzero block in column q
            }
            uint32_t blocksAdded = 0;
            for (int s = 1; s <= maxBlocks; s++)
            {
                //check if the blocks match
                uint32_t offsetA = A->offsets[indexA];
                uint32_t offsetB = B->offsets[indexB];

                if (offsetA % maxBlocks == offsetB / maxBlocks)
                {
                    cscBMM2(A->list[indexA], B->list[indexB], result);
                    addMatrix(result, block, block);
                    indexB++;
                    indexA++;
                    blocksAdded = 1;
                }
                else if (offsetA % maxBlocks > offsetB / maxBlocks)
                    indexB++;

                else if (offsetA % maxBlocks < offsetB / maxBlocks)
                    indexA++;

                if (offsetA > blockY * maxBlocks || offsetB % maxBlocks > blockX)
                    break;
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
        }
    }

    C->size = A->size;
    C->totalBlocks = totalBlocks;
}

