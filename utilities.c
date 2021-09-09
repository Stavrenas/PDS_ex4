#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>  // sqrt
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include "utilities.h"
#include "mmio.h"


void swap(double *n1, double *n2) {
    double temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void swapInts(int *n1, int *n2) {
    int temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void dividePoints(int n, int tasks, int *array) {
    int points = n / tasks;
    for (int i = 0; i < n; i++) {
        array[i] = points;
    }

    int pointsLeft = n % tasks;
    for (int i = 0; pointsLeft > 0; i++) {
        array[i]++;
        pointsLeft--;
    }
}

int findDestination(int id, int NumTasks) {
    if (id == NumTasks - 1)
        return 0;
    else
        return (id + 1);
}

int findSender(int id, int NumTasks) {
    if (id == 0)
        return (NumTasks - 1);
    else
        return (id - 1);
}

struct timeval tic()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv;
}

double toc(struct timeval begin)
{
    struct timeval end;
    gettimeofday(&end, NULL);
    double stime = ((double) (end.tv_sec - begin.tv_sec) * 1000 ) +
                   ((double) (end.tv_usec - begin.tv_usec) / 1000 );
    stime = stime / 1000;
    return(stime);
}


	
	