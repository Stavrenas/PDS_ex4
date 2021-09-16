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
    MPI_Init(NULL, NULL);

    int world_size, world_rank, name_len;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);

    char filename[] = "mycielskian.mtx";

    //dblp 1200 sec
    //mycielskian 100 sec

    Matrix *A = malloc(sizeof(Matrix));
    Matrix *B = malloc(sizeof(Matrix));
    Matrix *res = malloc(sizeof(Matrix));

    readMatrix(filename, A);
    readMatrix(filename, B);

    struct timeval start = tic();

    for (int i = 0; i < 1; i++)
        cscBMM(A, B, res);

    printf("Time for mult: %f\n", toc(start));

    //printMatrix(res);

    // Blocking algorithms: BCSR or CSB

    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);

    MPI_Finalize();
}

void cscBMM(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, *temp, indexB, last;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    temp = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    uint32_t idx_size = 0;
    uint32_t end = 0;

    c_elem[0] = 0;
    last = -1; //the first element is set to -1 and is unused. Thats why you see elements+1 on temp
    int perc = 0;
    for (int row = 1; row <= A->size; row++)
    { //go to each row of mtr A

        // if (row * 100 / A->size != perc)
        // {
        //     printf("%d%% \n", perc);
        //     perc = row * 100 / A->size;
        // }
        elements = 0;
        c_elem[row] = c_elem[row - 1] + elements;

        for (int col = 1; col <= B->size; col++)
        { //go to each col of mtr Β

            for (int a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
            { //go to each element in row "row" of mtr A

                int indexA = A->csc_idx[a];

                for (int b = B->csc_elem[indexA - 1]; b < B->csc_elem[indexA]; b++)
                { //check if there is a match in col "col"of mtr B
                    indexB = B->csc_idx[b];

                    if (indexB > col)
                        break;

                    else if (indexB == col)
                    {
                        if (last == col) //check if the element is already added
                        {
                            continue;
                            //do not add it
                        }

                        else
                        {
                            temp[elements] = col;
                            last = col;
                            elements++;
                        }
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
            end += elements;

            for (int i = end - elements; i < end; i++)
                c_idx[i] = temp[i - end + elements];
        }
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}
