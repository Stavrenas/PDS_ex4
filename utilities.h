#ifndef UTILITIES_H
#define UTILITIES_H

void swap(double *n1, double *n2);

void swapInts(int *n1, int *n2);

void dividePoints(int n, int tasks, int *array);

int findDestination(int id, int NumTasks);

int findSender(int id, int NumTasks);

#endif //UTILITIES_H