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

    //dblp 1200 sec
    //mycielskian 116 sec

    Matrix *A = malloc(sizeof(Matrix));
    Matrix *B = malloc(sizeof(Matrix));
    Matrix *res = malloc(sizeof(Matrix));

    char filename[] = "12.mtx";
    readMatrix(filename, A);
    readMatrix(filename, B);

    struct timeval start = tic();

    for (int i = 0; i < 1; i++)
        cscBMMparallel(A, B, res);

    printf("Time for mult: %f\n", toc(start));

    printMatrix(res);

    // Blocking algorithms: BCSR or CSB

    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);

    MPI_Finalize();
}

void cscBMM2(Matrix *A, Matrix *B, Matrix *C)
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

    int perc = 0;
    int sizeA, sizeB, start_a, end_a, start_b, end_b;

    sizeA = A->size;
    sizeB = B->size;
    for (int row = 1; row <= sizeA; row++)
    { //go to each row of mtr A

        // if (row * 100 / A->size != perc)
        // {
        //     printf("%d%% \n", perc);
        //     perc = row * 100 / A->size;
        // }
        //printf("Row %d of %d\n",row,A->size);

        for (int col = 1; col <= sizeB; col++)
        { //go to each col of mtr Β

            start_a = A->csc_elem[row - 1];
            end_a = A->csc_elem[row];
            for (int a = start_a; a < end_a; a++)
            { //go to each element in row "row" of mtr A

                indexA = A->csc_idx[a];
                start_b = B->csc_elem[indexA - 1];
                end_b = B->csc_elem[indexA];

                for (int b = start_b; b < end_b; b++)
                { //check if there is a match in col "col"of mtr B
                    indexB = B->csc_idx[b];

                    if (indexB > col)
                        break;

                    else if (indexB == col)
                    {
                        if (last == col) //check if the element is already added
                        {
                            //do not add it
                        }

                        else
                        {
                            c_idx[elements] = col;
                            elements++;
                            if (elements == idx_size)
                            {
                                idx_size *= 2;
                                c_idx = realloc(c_idx, idx_size * sizeof(uint32_t));
                            }

                            last = col;
                        }
                        break;
                    }
                }
            }
        }

        c_elem[row] = elements;
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
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

void cscBMMparallel(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *elements, *c_idx;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    elements = (uint32_t *)malloc((A->size) * sizeof(uint32_t));
    uint32_t idx_size = 0;
    uint32_t end = 0;

    c_elem[0] = 0;

    int perc = 0;
    int sizeA, sizeB, start_a, end_a, start_b, end_b;
    sizeA = A->size;

    //printf("thread %d, of %d\n",omp_get_thread_num(),omp_get_num_threads());
#pragma omp parallel
    {
        uint32_t *temp, indexB, indexA, last, localElements, tempSize;
        temp = (uint32_t *)malloc(sizeof(uint32_t));
        tempSize = 1;
        last = -1;
        localElements = 0;

        int nthreads = omp_get_num_threads();
        int id = omp_get_thread_num();
        //printf("Hello from %d \n", id);

        for (int row = 1 + id; row <= A->size; row += nthreads)
        { //go to each row of mtr A

            // if (row * 100 / A->size != perc)
            // {
            //     printf("%d%% \n", perc);
            //     perc = row * 100 / A->size;
            // }
            // printf("Row %d of %d\n", row, A->size);

            localElements = 0;

            for (int col = 1; col <= B->size; col++)
            { //go to each col of mtr Β

                for (int a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
                { //go to each element in row "row" of mtr A

                    indexA = A->csc_idx[a];

                    for (int b = B->csc_elem[indexA - 1]; b < B->csc_elem[indexA]; b++)
                    { //check if there is a match in col "col"of mtr B
                        indexB = B->csc_idx[b];

                        if (indexB > col)
                            break;

                        else if (indexB == col)
                        {
                            if (last == col) //check if the element is already added
                            {
                                //do not add it
                            }

                            else
                            {
                                temp[localElements] = col;
                                localElements++;
                                if (localElements == tempSize)
                                {
                                    tempSize *= 2;
                                    temp = realloc(temp, tempSize * sizeof(uint32_t));
                                }

                                last = col;
                            }
                            break;
                        }
                    }
                }
            }

            elements[row] = localElements;
        }

#pragma omp barrier //sync
        if (id == 0)
        {
            for (int row = 1; row <= A->size; row++)
                c_elem[row] = c_elem[row - 1] + elements[row];
            c_idx = realloc(c_idx, c_elem[A->size] * sizeof(uint32_t));
        }
#pragma omp barrier //sync

        int index = 0;
        for (int row = 1 + id; row <= A->size; row += nthreads)
        {

            for (int j = c_elem[row - 1]; j < c_elem[row], index < localElements; j++, index++)
                c_idx[j] = temp[index];
        }
    }
    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}
