#ifndef OMP_UTILITIES_H
#define OMP_UTILITIES_H
#include "types.h"

void multMatrixParallel(Matrix *A, Matrix *B, Matrix *C);

void multBlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C);

void blockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix);

void unblockMatrix(BlockedMatrix *blockedMatrix, Matrix *mtr);

void addΒlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C);

#endif //UTILITIES_H
