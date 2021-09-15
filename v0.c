#include <stdio.h>
#include <stdlib.h>
#include <math.h> // sqrt
#include <mpi.h>
#include <string.h>
#include "utilities.h"
#include "controller.h"
#include "read.h"
#include "mmio.h"

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    int world_size, world_rank, name_len;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);

    char filename[] = "12.mtx";
    Matrix *mtr = malloc(sizeof(Matrix));
    readMatrix(filename, mtr);

    cscSymmetricBMM(mtr, mtr, mtr);

    // Blocking algorithms: BCSR or CSB

    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);

    MPI_Finalize();
}

void cscSymmetricBMM(Matrix *A, Matrix *B, Matrix *C)
{
    for (int row = 1; row <= A->size; row++)
    { //go to each row of mtr A
        for (int col = 1; col <= B->size; col++)
        { //go to each col of mtr Β

            for (int a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
            { //go to each element in row of mtr A

                //printf("\n(%d ,%d):", A->csc_idx[a], row);

                for (int b = B->csc_elem[A->csc_idx[a] -1 ]; b < B->csc_elem[A->csc_idx[a]]; b++)
                {
                    
                    if (B->csc_idx[b] == col){

                        printf("HIT %d,%d elements %d and %d\n", row, col,a,b);
                        a = A->csc_elem[row];
                        b = B->csc_elem[A->csc_idx[a]];
                        continue;
                    }

                    else if (B->csc_idx[b] > col){
                        A->csc_elem[row];
                        b = B->csc_elem[A->csc_idx[a]];
                        continue;
                    }
                }
            }

        }
    }
}
