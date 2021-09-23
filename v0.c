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

    char filenameA[] = "12.mtx";
    char filenameB[] = "12a.mtx";

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    struct timeval start = tic();
    for (int i = 0; i < 10000000; i++)
        addMatrix(A, B, C);

    //blockMatrix(A, 100, blockA);
    printf("Time for add: %f\n", toc(start));

    //printBlockedMatrix(blockA);

    printMatrix(C);

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

void addMatrix(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, idx_size;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));

    idx_size = 1;
    c_elem[0] = 0;
    elements = 0;

    uint32_t size, start_a, end_a, start_b, end_b, indexA, indexB;

    size = A->size;

    for (uint32_t row = 1; row <= size; row++)
    { //go to each row of mtr A
        //printf("row %d\n",row);

        start_a = A->csc_elem[row - 1];
        end_a = A->csc_elem[row];

        start_b = B->csc_elem[row - 1];
        end_b = B->csc_elem[row];

        for (uint32_t a = start_a, b = start_b;;)
        { //go to each element in row of mtr A and mtr B

            if (a == end_a && b == end_b)
                break;
            else
            {
                indexA = A->csc_idx[a];
                indexB = B->csc_idx[b];

                //printf("a is %d and b is %d, indexA is %d and indexB is %d\n", a, b, indexA, indexB);

                if (indexA < indexB)
                {
                    if (a != end_a)
                    {
                        c_idx[elements] = indexA;
                        elements++;
                        a++;
                    }
                    else
                    {
                        c_idx[elements] = indexB;
                        elements++;
                        b++;
                    }
                }

                else if (indexA > indexB)
                {
                    if (b != end_b)
                    {
                        c_idx[elements] = indexB;
                        elements++;
                        b++;
                    }
                    else
                    {
                        c_idx[elements] = indexA;
                        elements++;
                        a++;
                    }
                }

                else
                {
                    c_idx[elements] = indexB;
                    elements++;
                    if (a != end_a)
                        a++;
                    if (b != end_b)
                        b++;
                }

                if (elements == idx_size)
                {
                    idx_size *= 2;
                    c_idx = realloc(c_idx, idx_size * sizeof(uint32_t));
                }
            }
        }

        c_elem[row] = elements;
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
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

int findIndex(BlockedMatrix *mtr, uint32_t indx)
{
    return binarySeach(mtr->offsets, 0, mtr->totalBlocks - 1, indx);
}

uint32_t binarySeach(uint32_t *list, uint32_t left, uint32_t right, uint32_t index)
{
    if (right >= left)
    {
        int mid = left + (right - left) / 2;

        // If the element is present at the middle
        // itself
        if (list[mid] == index)
            return mid;

        // If element is smaller than mid, then
        // it can only be present in left subarray
        if (list[mid] > index)
            return binarySearch(list, left, mid - 1, index);

        // Else the element can only be present
        // in right subarray
        return binarySearch(list, mid + 1, right, index);
    }

    // We reach here when element is not
    // present in array
    return -1;
}