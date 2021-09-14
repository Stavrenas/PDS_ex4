#ifndef READ_H
#define READ_H

typedef struct
{
    uint32_t *csc_row;
    uint32_t *csc_col;
    uint32_t size;
} Matrix;


void readMatrix(char *file_path, Matrix* M);

void coo2csc(
    uint32_t       * const row,       /*!< CSC row start indices */
    uint32_t       * const col,       /*!< CSC column indices */
    uint32_t const * const row_coo,   /*!< COO row indices */
    uint32_t const * const col_coo,   /*!< COO column indices */
    uint32_t const         nnz,       /*!< Number of nonzero elements */
    uint32_t const         n,         /*!< Number of rows/columns */
    uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);

#endif //READ_H