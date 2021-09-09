#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>  // sqrt
#include <float.h>
#include "utilities.h"
#include "mmio.h"


void swap(double *n1, double *n2) {
    double temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void swapInts(int *n1, int *n2) {
    int temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void dividePoints(int n, int tasks, int *array) {
    int points = n / tasks;
    for (int i = 0; i < n; i++) {
        array[i] = points;
    }

    int pointsLeft = n % tasks;
    for (int i = 0; pointsLeft > 0; i++) {
        array[i]++;
        pointsLeft--;
    }
}

int findDestination(int id, int NumTasks) {
    if (id == NumTasks - 1)
        return 0;
    else
        return (id + 1);
}

int findSender(int id, int NumTasks) {
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
    double stime = ((double) (end.tv_sec - begin.tv_sec) * 1000 ) +
                   ((double) (end.tv_usec - begin.tv_usec) / 1000 );
    stime = stime / 1000;
    return(stime);
}

void coo2csc(
    uint32_t       * const row,       /*!< CSC row start indices */
    uint32_t       * const col,       /*!< CSC column indices */
    uint32_t const * const row_coo,   /*!< COO row indices */
    uint32_t const * const col_coo,   /*!< COO column indices */
    uint32_t const         nnz,       /*!< Number of nonzero elements */
    uint32_t const         n,         /*!< Number of rows/columns */
    uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
)
{

    for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


    for (uint32_t l = 0; l < nnz; l++)
        col[col_coo[l] - isOneBased]++;

    // ----- cumulative sum
    for (uint32_t i = 0, cumsum = 0; i < n; i++)
    {
        uint32_t temp = col[i];
        col[i] = cumsum;
        cumsum += temp;
    }
    col[n] = nnz;
    // ----- copy the row indices to the correct place
    for (uint32_t l = 0; l < nnz; l++)
    {
        uint32_t col_l;
        col_l = col_coo[l] - isOneBased;

        uint32_t dst = col[col_l];
        row[dst] = row_coo[l]+1;

        col[col_l]++;
    }
    // ----- revert the column pointers
    for (uint32_t i = 0, last = 0; i < n; i++)
    {
        uint32_t temp = col[i];
        col[i] = last;
        last = temp;
    }

}

	
	
	