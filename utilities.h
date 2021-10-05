#ifndef UTILITIES_H
#define UTILITIES_H
#include "types.h"

struct timeval tic();

double toc(struct timeval begin);

void multMatrix(Matrix *A, Matrix *B, Matrix *C);

void multMatrixMasked(Matrix *A, Matrix *B, Matrix *C, Matrix *mask);

void multMatrix2(Matrix *A, Matrix *B, Matrix *C);

void addMatrix(Matrix *A, Matrix *B, Matrix *C);

int binarySearch(uint32_t *list, uint32_t left, uint32_t right, uint32_t index);

int findIndex(BlockedMatrix *mtr, uint32_t indx);

#endif //UTILITIES_H
