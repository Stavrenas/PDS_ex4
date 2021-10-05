#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include "utilities.h"
#include "mmio.h"
#include "read.h"

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

void multMatrix(Matrix *A, Matrix *B, Matrix *C)
{
    //alocate memory for result matrix
    uint32_t *c_elem, *c_idx, elements, *temp, indexB, last;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    temp = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    uint32_t idx_size = 0;
    uint32_t end = 0;

    c_elem[0] = 0;

    for (uint32_t row = 1; row <= A->size; row++)
    { //go to each row of mtr A

        elements = 0;
        last = -1; //the first element is set to -1 and is used to avoid adding the same element twice
        c_elem[row] = c_elem[row - 1] + elements;

        for (uint32_t col = 1; col <= B->size; col++)
        { //go to each col of mtr Β

            for (uint32_t a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
            { //go to each element in "row"-th row  of mtr A

                uint32_t indexA = A->csc_idx[a];

                for (uint32_t b = B->csc_elem[indexA - 1]; b < B->csc_elem[indexA]; b++)
                { //check if there is a match in col "col" of mtr B
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
        { //add the indexes and realloc

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

void multMatrix2(Matrix *A, Matrix *B, Matrix *C)
{
    //allocate memory for result mult
    uint32_t *c_elem, *c_idx, elements, *temp, indexB, indexA, last;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    temp = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    uint32_t idx_size = 1;
    uint32_t end = 0;

    c_elem[0] = 0;

    elements = 0;

    //variables used in for loops
    uint32_t sizeA, sizeB, start_a, end_a, start_b, end_b;

    sizeA = A->size;
    sizeB = B->size;

    for (uint32_t row = 1; row <= sizeA; row++)
    {              //go to each row of mtr A
        last = -1; //last element added in each row
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
                { //check if there is a match in col "col" of mtr B
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
                            { //faster realloc than multMatrix
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

void addMatrix(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, idx_size;
    //initialize result mult

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    idx_size = 1; //used in realloc
    c_elem[0] = 0;
    elements = 0;

    //variables used in for loops
    uint32_t size, start_a, end_a, start_b, end_b, indexA, indexB;
    size = A->size;

    if (B->csc_elem[B->size] == 0 && A->csc_elem[A->size] == 0)
    {
        for (int i = 0; i <= A->size; i++)
            c_elem[i] = 0;
    }
    else
    {

        for (uint32_t row = 1; row <= size; row++)
        { //go to each row of mtr C

            start_a = A->csc_elem[row - 1];
            end_a = A->csc_elem[row];

            start_b = B->csc_elem[row - 1];
            end_b = B->csc_elem[row];

            for (uint32_t a = start_a, b = start_b;;)
            { //  merge  current row of A and  B

                if (a == end_a && b == end_b)
                    break;
                else
                {
                    if (a != end_a)
                        indexA = A->csc_idx[a];
                    if (b != end_b)
                        indexB = B->csc_idx[b];

                    if (indexA < indexB)
                    {
                        //check if row has elements
                        if (a != end_a)
                        {
                            c_idx[elements] = indexA;
                            elements++;
                            a++;
                        }
                        //check if row has elements
                        else if (b != end_b)
                        {
                            c_idx[elements] = indexB;
                            elements++;
                            b++;
                        }
                    }

                    else if (indexA > indexB)
                    {
                        //check if row has elements
                        if (b != end_b)
                        {
                            c_idx[elements] = indexB;
                            elements++;
                            b++;
                        }
                        //check if row has elements
                        else if (a != end_a)
                        {
                            c_idx[elements] = indexA;
                            elements++;
                            a++;
                        }
                    }

                    else
                    { //if indexA = indexB
                        c_idx[elements] = indexB;
                        elements++;
                        if (a != end_a)
                            a++;
                        if (b != end_b)
                            b++;
                    }

                    if (elements == idx_size)
                    {
                        idx_size++;
                        c_idx = realloc(c_idx, idx_size * sizeof(uint32_t));
                    }
                }
            }

            c_elem[row] = elements;
        }
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}


int findIndex(BlockedMatrix *mtr, uint32_t indx)
{
    //check if a block with offset = indx exists
    return binarySearch(mtr->offsets, 0, mtr->totalBlocks - 1, indx);
}

int binarySearch(uint32_t *list, uint32_t left, uint32_t right, uint32_t index)
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
