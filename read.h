#ifndef READ_H
#define READ_H

typedef struct
{
    uint32_t *csc_idx;  //column indexes
    uint32_t *csc_elem;
    uint32_t size;
    uint32_t sizeX;
    uint32_t sizeY;
} Matrix;

typedef struct
{
    uint32_t size;
    uint32_t* offsets;
    Matrix** list;
} BlockedMatrix;

void cscBMM(Matrix *A, Matrix *B, Matrix *C);

void cscBMM2(Matrix *A, Matrix *B, Matrix *C);

void cscBMMparallel(Matrix *A, Matrix *B, Matrix *C);

void readMatrix(char *file_path, Matrix* M);

void saveMatrix(Matrix *res, char* filename);

void printMatrix(Matrix *res);

void printBlockedMatrix(BlockedMatrix *res);

void blockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blocked);

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