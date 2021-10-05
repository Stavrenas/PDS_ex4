#ifndef MPI_UTILITIES_H
#define MPI_UTILITIES_H
#include "types.h"

void multBlockedMatrixMPI(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C);

Matrix *MPI_Mult(BlockedMatrix *A, BlockedMatrix *B);

#endif //MPI_UTILITIES_H
