#ifndef UTILITIES_H
#define UTILITIES_H
#include "types.h"

void swap(double *n1, double *n2);

void swapInts(int *n1, int *n2);

void dividePoints(int n, int tasks, int *array);

int findDestination(int id, int NumTasks);

int findSender(int id, int NumTasks);

struct timeval tic();

double toc(struct timeval begin);

void multMatrix(Matrix *A, Matrix *B, Matrix *C);

void multMatrix2(Matrix *A, Matrix *B, Matrix *C);

void multMatrixParallel(Matrix *A, Matrix *B, Matrix *C);

void blockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix);

void addMatrix(Matrix *A, Matrix *B, Matrix *C);

uint32_t binarySeach(uint32_t* list, uint32_t left, uint32_t right,uint32_t index );

void unblockMatrix( BlockedMatrix *blockedMatrix,Matrix *mtr);

int binarySearch(uint32_t *list, uint32_t left, uint32_t right, uint32_t index);

int findIndex(BlockedMatrix *mtr, uint32_t indx);

void multBlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C);

void multBlockedMatrixMPI(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C,  uint32_t *rows, uint32_t rows_size);

#endif //UTILITIES_H
