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

void blockBMM(BlockedMatrix *A, BlockedMatrix *B, Matrix *C)
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

        c_elem[row] = elements;
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}
