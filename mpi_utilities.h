#ifndef MPI_UTILITIES_H
#define MPI_UTILITIES_H
#include "types.h"

void multBlockedMatrixMPI(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C);

void multBlockedMatrixMPIMasked(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C, BlockedMatrix *mask);

Matrix *MPI_Mult(BlockedMatrix *A, BlockedMatrix *B);

Matrix *MPI_MultMasked(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *mask);
#endif //MPI_UTILITIES_H
