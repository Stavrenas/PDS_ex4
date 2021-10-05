#ifndef OMP_UTILITIES_H
#define OMP_UTILITIES_H
#include "types.h"

void multMatrixParallel(Matrix *A, Matrix *B, Matrix *C);

void multMatrixParallelMasked(Matrix *A, Matrix *B, Matrix *C,Matrix* mask);

void multBlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C);

void multBlockedMatrixMasked(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C,BlockedMatrix *mask);

void blockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix);

void unblockMatrix(BlockedMatrix *blockedMatrix, Matrix *mtr);

void addBlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C);

#endif //UTILITIES_H
