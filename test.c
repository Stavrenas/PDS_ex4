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

    char matrix[] = "mycielskian";
    char matrix2[] = "12";
    char *filenameA = (char *)malloc(40 * sizeof(char));
    char *filenameB = (char *)malloc(40 * sizeof(char));
    char *name = (char *)malloc(40 * sizeof(char));
    int blocksize = 250;
    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    // blockMatrix(C, blocksize, blockC);
    // printBlockedMatrix(blockC);
    struct timeval start = tic();
    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);
    printf("Block time : %f\n", toc(start));
    // addΒlockedMatrix(blockA, blockB, blockC);
    // printBlockedMatrix(blockC);
    // unblockMatrix(blockC, C);
    // sprintf(name, "%s_block.txt", matrix);
    // saveMatrix(C, name);
    start = tic();
    multBlockedMatrix(blockA, blockA, blockC);
    printf("Mult time : %f\n", toc(start));
    start = tic();
    //printBlockedMatrix(blockC);
    unblockMatrix(blockC, C);
    printf("Unblock time : %f\n", toc(start));

    sprintf(name, "%s_normal.txt", matrix);
    saveMatrix(C, name);
}
