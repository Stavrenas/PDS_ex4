#include <stdio.h>
#include <stdlib.h>
#include <math.h> // sqrt
#include <string.h>
#include "utilities.h"
#include "controller.h"
#include "read.h"
#include "mmio.h"

int main(int argc, char **argv)
{

    Matrix *A = (Matrix *)malloc(sizeof(Matrix));
    Matrix *B = (Matrix *)malloc(sizeof(Matrix));
    Matrix *C = (Matrix *)malloc(sizeof(Matrix));

    BlockedMatrix *blockA = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockB = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockC = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    // BlockedMatrix *blockResult = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    char matrix[] = "12";
    char *filenameA = (char *)malloc(25 * sizeof(char));
    char *filenameB = (char *)malloc(25 * sizeof(char));
    char *name = (char *)malloc(25 * sizeof(char));
    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    int blocksize = 4;
    printf("Elements before blocking: %d\n",A->csc_elem[A->size]);
    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);
    unblockMatrix2(blockA,C);
    printMatrix(C);
    printMatrix(A); 
}
