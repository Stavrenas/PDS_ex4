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
    // else
        // printf("\n\n***Multipling %s with a blocksize of %d***\n\n", matrix, blocksize);

    char *filenameA = (char *)malloc(40 * sizeof(char));
    char *filenameB = (char *)malloc(40 * sizeof(char));
    char *name = (char *)malloc(40 * sizeof(char));

    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    struct timeval start = tic();
    struct timeval total = tic();
    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);
    // printf("Block time : %f\n", toc(start));

    start = tic();
    // multBlockedMatrix(blockA, blockA, temp);
    // unblockMatrix(temp,C);
    // sprintf(name, "%s_blocked.txt", matrix);
    // saveMatrix(C, name);
    // blockMask = temp;

    multBlockedMatrixMasked(blockA, blockA, blockC, blockA);
    // printf("Mult time : %f\n", toc(start));
    // start = tic();

    unblockMatrix(blockC, C);
    // printf("Unblock time : %f\n", toc(start));
    // printf("Total time : %f\n", toc(total));
    printf("%lf", toc(start));
    sprintf(name, "%s_blockedMasked.txt", matrix);
    saveMatrix(C, name);
}
