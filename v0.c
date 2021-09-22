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

    // free memory

    clearMatrix(A);
    clearMatrix(B);
    clearBlockedMatrix(blockA);
}

void addCSC(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, *temp, indexB, indexA, last;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    temp = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    uint32_t idx_size = 1;
    uint32_t end = 0;

    c_elem[0] = 0;
    last = -1;
    elements = 0;

    uint32_t sizeA, sizeB, start_a, end_a, start_b, end_b, indexA, indexB;

    sizeA = A->size;
    sizeB = B->size;

    for (uint32_t row = 1; row <= sizeA; row++)
    { //go to each row of mtr A

        start_a = A->csc_elem[row - 1];
        end_a = A->csc_elem[row];

        start_b = B->csc_elem[row - 1];
        end_b = B->csc_elem[row - 1];

        for (uint32_t a = start_a, b = start_b; a < end_a, b < end_b;)
        { //go to each element in row of mtr A and mtr B

            indexA = A->csc_idx[a];
            indexB = B->csc_idx[b];

            if (indexA < indexB)
            {
                c_idx[elements] = indexA;
                a++;
            }

            else if (indexB < indexA)
            {
                c_idx[elements] = indexB;
                b++;
            }

            else
            {
                c_idx[elements] = indexB;
                b++;
                a++;
            }
            
            elements++;
            if (elements == idx_size)
            {
                idx_size *= 2;
                c_idx = realloc(c_idx, idx_size * sizeof(uint32_t));
            }
        }

        c_elem[row] = elements;
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}
