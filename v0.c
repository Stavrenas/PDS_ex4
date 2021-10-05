#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utilities.h"
#include "read.h"
#include "mmio.h"

int main(int argc, char **argv)
{

    Matrix *A = (Matrix *)malloc(sizeof(Matrix));
    Matrix *B = (Matrix *)malloc(sizeof(Matrix));
    Matrix *C = (Matrix *)malloc(sizeof(Matrix));
    Matrix *mask = (Matrix *)malloc(sizeof(Matrix));

    char *matrix = (char *)malloc(40 * sizeof(char));

    if (argc == 1)
        sprintf(matrix, "%s", "mycielskian");

    else if (argc == 2)
        sprintf(matrix, "%s", argv[1]);
    else
    {
        printf("Usage: ./v0 matrix_name \n");
        exit(-1);
    }

    // printf("\n\n***Multipling %s ***\n\n", matrix);

    char *filenameA = (char *)malloc(40 * sizeof(char));
    char *filenameB = (char *)malloc(40 * sizeof(char));
    char *name = (char *)malloc(40 * sizeof(char));

    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    struct timeval start = tic();

    multMatrix(A, B, C);
    sprintf(name, "%s_serialMasked.txt", matrix);
    saveMatrix(C, name);
    
    mask = C;
    multMatrixMasked(A, B, C, mask);
    printf("Serial mult time : %f\n", toc(start));

    sprintf(name, "%s_serialMasked.txt", matrix);
    saveMatrix(C, name);
}
