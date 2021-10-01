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

    char matrix[] = "50";
    char matrix2[] = "12";
    char *filenameA = (char *)malloc(25 * sizeof(char));
    char *filenameB = (char *)malloc(25 * sizeof(char));
    char *name = (char *)malloc(25 * sizeof(char));
    int blocksize = 4;
    sprintf(filenameA, "%s.mtx", matrix);
    sprintf(filenameB, "%s.mtx", matrix);

    readMatrix(filenameA, A);
    readMatrix(filenameB, B);

    // blockMatrix(C, blocksize, blockC);
    // printBlockedMatrix(blockC);

    blockMatrix(A, blocksize, blockA);
    blockMatrix(B, blocksize, blockB);

    // addΒlockedMatrix(blockA, blockB, blockC);
    // printBlockedMatrix(blockC);
    // unblockMatrix(blockC, C);
    // sprintf(name, "%s_block.txt", matrix);
    // saveMatrix(C, name);

    multBlockedMatrix(blockA, blockA, blockC);
    printBlockedMatrix(blockC);
    unblockMatrix(blockC, C);
    sprintf(name, "%s_normal.txt", matrix);
    saveMatrix(C, name);
}
