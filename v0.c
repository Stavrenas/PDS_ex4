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

    printf("\n");
    char filename[] = "12.mtx";
    uint32_t *csc_row, *csc_col;  //csc_col[i] -> total elements up to ith row | csc_row[i] -> column index of ith element
    int n, nnz;

    readMatrix(&csc_row, &csc_col, &n, filename);
    nnz = csc_col[n-1];
    printf("\n");

    for (int row = 1; row <= n; row++) 
    {
        for (int i1 = csc_col[row - 1]; i1 < csc_col[row]; i1++)   //csc_row[i1] -> column
            printf("(%d ,%d) ", csc_row[i1], row);
            printf("\n");

    }


    // Print off a hello world message
    //printf("Hello world from processor %s, rank %d out of %d processors\n",processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}

