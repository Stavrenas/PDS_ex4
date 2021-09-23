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

void cscBMM(Matrix *A, Matrix *B, Matrix *C);

void cscBMM2(Matrix *A, Matrix *B, Matrix *C);

void cscBMMparallel(Matrix *A, Matrix *B, Matrix *C);

void blockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix);

void addMatrix(Matrix *A, Matrix *B, Matrix *C);

uint32_t binarySeach(uint32_t* list, uint32_t left, uint32_t right,uint32_t index );

void unblockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix);

int binarySearch(uint32_t *list, uint32_t left, uint32_t right, uint32_t index);

int findIndex(BlockedMatrix *mtr, uint32_t indx);

void blockBMM(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C);

#endif //UTILITIES_H
