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
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);


    char filename[] = "12.mtx";
    Matrix* mtr = malloc(sizeof(Matrix));
    readMatrix(filename, mtr);


    printf("\nSize is %d\n",mtr->size);


    // Blocking algorithms: BCSR or CSB

    for (int row = 1; row <= mtr->size; row++) 
    {
        for (int i1 = mtr->csc_col[row - 1]; i1 < mtr->csc_col[row]; i1++)   //csc_row[i1] -> column
            printf("(%d ,%d) ", mtr->csc_row[i1], row);
            printf("\n");

    }


    // Print off a hello world message
    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}

void csrBMM(Matrix A, Matrix B, Matrix* C){

}

