#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utilities.h"
#include "omp_utilities.h"
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
    BlockedMatrix *temp = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
    BlockedMatrix *blockMask = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));

    char *matrix = (char *)malloc(40 * sizeof(char));
    int blocksize = 1;

    if (argc == 1)
    {
        sprintf(matrix, "%s", "mycielskian");
        blocksize = 250;
    }
    else if (argc == 3)
    {
        sprintf(matrix, "%s", argv[1]);
        blocksize = atoi(argv[2]);
    }

    if (argc != 1 && argc != 3)
        printf("Usage: ./v2 matrix_name blocksize\n");

    char *filenameA = (char *)malloc(40 * sizeof(char));
    char *filenameB = (char *)malloc(40 * sizeof(char));
    char *name = (char *)malloc(40 * sizeof(char));

    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    struct timeval start = tic();

    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);

    blockMask = blockA;

    multBlockedMatrixMasked(blockA, blockA, blockC, blockMask);
    unblockMatrix(blockC, C);

    printf("Block mult time : %f\n", toc(start));
    sprintf(name, "%s_blocked.txt", matrix);
    saveMatrix(C, name);
}
