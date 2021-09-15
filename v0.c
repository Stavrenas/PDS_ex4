#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h> // sqrt
#include <mpi.h>
#include <string.h>
#include "utilities.h"
#include "controller.h"
#include "read.h"
#include "mmio.h"

int main(int argc, char **argv)
{
  // this is me
    MPI_Init(NULL, NULL);

    int world_size, world_rank, name_len;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);

    char filename[] = "12.mtx";
    Matrix *mtr = malloc(sizeof(Matrix));
    Matrix *res = malloc(sizeof(Matrix));
    readMatrix(filename, mtr);

    struct timeval start = tic();
    for (int i = 0; i < 10000000; i++)
        cscSymmetricBMM(mtr, mtr, res);

    printf("time for mult: %f\n", toc(start));

    printMatrix(res);

    // Blocking algorithms: BCSR or CSB

    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);

    MPI_Finalize();
}

void cscSymmetricBMM(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, idx_size, new_size, *temp;
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    idx_size = 0 * sizeof(uint32_t);

    bool added;

    temp = malloc(A->size * sizeof(uint32_t));

    c_elem[0] = 0;

    for (int row = 1; row <= A->size; row++)
    { //go to each row of mtr A

        elements = 0;
        c_elem[row] = c_elem[row - 1] + elements;
        for (int col = 1; col <= B->size; col++)
        { //go to each col of mtr Β

            for (int a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
            { //go to each element in row of mtr A

                //printf("\n(%d ,%d):", A->csc_idx[a], row);
                int index = A->csc_idx[a];

                for (int b = B->csc_elem[index - 1]; b < B->csc_elem[index]; b++)
                {
                    added = false;

                    if (B->csc_idx[b] > col)
                        break;

                    else if (B->csc_idx[b] == col)
                    {
                        if (elements > 1 && temp[elements - 1] == col)
                            break;

                        temp[elements] = col;
                        elements++;
                        added = true;
                        break;
                    }
                }
            }
        }
        c_elem[row] = c_elem[row - 1] + elements;

        if (elements != 0)
        {
            idx_size += elements * sizeof(uint32_t); //in bytes
            c_idx = realloc(c_idx, idx_size);
            uint32_t end = idx_size / sizeof(uint32_t);

            for (int i = end - elements; i < end; i++)
                c_idx[i] = temp[i - end + elements];
        }
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}
