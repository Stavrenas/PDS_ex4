#ifndef UTILITIES_H
#define UTILITIES_H

void swap(double *n1, double *n2);

void swapInts(int *n1, int *n2);

void dividePoints(int n, int tasks, int *array);

int findDestination(int id, int NumTasks);

int findSender(int id, int NumTasks);

struct timeval tic();

double toc(struct timeval begin);

void coo2csc(
    uint32_t       * const row,       /*!< CSC row start indices */
    uint32_t       * const col,       /*!< CSC column indices */
    uint32_t const * const row_coo,   /*!< COO row indices */
    uint32_t const * const col_coo,   /*!< COO column indices */
    uint32_t const         nnz,       /*!< Number of nonzero elements */
    uint32_t const         n,         /*!< Number of rows/columns */
    uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);

#endif //UTILITIES_H