#ifndef UTILITIES_H
#define UTILITIES_H

void swap(double *n1, double *n2);

void swapInts(int *n1, int *n2);

void dividePoints(int n, int tasks, int *array);

int findDestination(int id, int NumTasks);

int findSender(int id, int NumTasks);

struct timeval tic();

double toc(struct timeval begin);

BlockedMatrix* blockMatrix(Matrix* mtr, uint32_t blockSize);
#endif //UTILITIES_H