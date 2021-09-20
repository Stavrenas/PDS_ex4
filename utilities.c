#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h> // sqrt
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include "utilities.h"
#include "mmio.h"
#include "read.h"

void swap(double *n1, double *n2)
{
    double temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void swapInts(int *n1, int *n2)
{
    int temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void dividePoints(int n, int tasks, int *array)
{
    int points = n / tasks;
    for (int i = 0; i < n; i++)
    {
        array[i] = points;
    }

    int pointsLeft = n % tasks;
    for (int i = 0; pointsLeft > 0; i++)
    {
        array[i]++;
        pointsLeft--;
    }
}

int findDestination(int id, int NumTasks)
{
    if (id == NumTasks - 1)
        return 0;
    else
        return (id + 1);
}

int findSender(int id, int NumTasks)
{
    if (id == 0)
        return (NumTasks - 1);
    else
        return (id - 1);
}

struct timeval tic()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv;
}

double toc(struct timeval begin)
{
    struct timeval end;
    gettimeofday(&end, NULL);
    double stime = ((double)(end.tv_sec - begin.tv_sec) * 1000) +
                   ((double)(end.tv_usec - begin.tv_usec) / 1000);
    stime = stime / 1000;
    return (stime);
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
    uint32_t sizeA, sizeB, start_a, end_a, start_b, end_b;

    sizeA = A->size;
    sizeB = B->size;
    for (uint32_t row = 1; row <= sizeA; row++)
    { //go to each row of mtr A

        // if (row * 100 / A->size != perc)
        // {
        //     printf("%d%% \n", perc);
        //     perc = row * 100 / A->size;
        // }
        //printf("Row %d of %d\n",row,A->size);

        for (uint32_t col = 1; col <= sizeB; col++)
        { //go to each col of mtr Β

            start_a = A->csc_elem[row - 1];
            end_a = A->csc_elem[row];
            for (uint32_t a = start_a; a < end_a; a++)
            { //go to each element in row "row" of mtr A

                indexA = A->csc_idx[a];
                start_b = B->csc_elem[indexA - 1];
                end_b = B->csc_elem[indexA];

                for (uint32_t b = start_b; b < end_b; b++)
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
    last = -1; //the first element is set to -1 and is used to avoid adding the same element twice

    for (uint32_t row = 1; row <= A->size; row++)
    { //go to each row of mtr A

        elements = 0;
        c_elem[row] = c_elem[row - 1] + elements;

        for (uint32_t col = 1; col <= B->size; col++)
        { //go to each col of mtr Β

            for (uint32_t a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
            { //go to each element in row "row" of mtr A

                uint32_t indexA = A->csc_idx[a];

                for (uint32_t b = B->csc_elem[indexA - 1]; b < B->csc_elem[indexA]; b++)
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

            for (uint32_t i = end - elements; i < end; i++)
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

    c_elem[0] = 0;
    int perc = 0;

#pragma omp parallel
    {
        uint32_t *temp, indexB, indexA, last, localElements, totalElements, tempSize;
        temp = (uint32_t *)malloc(sizeof(uint32_t));
        tempSize = 1;
        last = -1;

        localElements = 0;
        totalElements = 0;

        int nthreads = omp_get_num_threads();
        int id = omp_get_thread_num();

        for (uint32_t row = 1 + id; row <= A->size; row += nthreads)
        { //go to each row of mtr A

            localElements = 0;

            for (uint32_t col = 1; col <= B->size; col++)
            { //go to each col of mtr Β

                for (uint32_t a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
                { //go to each element in row "row" of mtr A

                    indexA = A->csc_idx[a];

                    for (uint32_t b = B->csc_elem[indexA - 1]; b < B->csc_elem[indexA]; b++)
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
                                temp[totalElements] = col; //replaced local elements with total elements
                                totalElements++;
                                if (totalElements == tempSize)
                                {
                                    tempSize *= 2;
                                    temp = realloc(temp, tempSize * sizeof(uint32_t));
                                }
                                localElements++;
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
            for (uint32_t row = 1; row <= A->size; row++)
                c_elem[row] = c_elem[row - 1] + elements[row];
            c_idx = realloc(c_idx, c_elem[A->size] * sizeof(uint32_t));
        }

#pragma omp barrier //sync

        uint32_t index = 0;
        uint32_t start, end;
        for (uint32_t row = 1 + id; row <= A->size; row += nthreads)
        {

            start = c_elem[row - 1];
            end = c_elem[row];

            for (uint32_t j = start; j < end; j++)
                c_idx[j] = temp[j - start];
        }
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}




