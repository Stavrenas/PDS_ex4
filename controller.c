#include <stdio.h> // fprintf, printf
#include <stdlib.h> // EXIT_FAILURE, EXIT_SUCCESS
#include <stdint.h>
#include <mpi.h>


void saveStats(char *name, double time, int n, int d, int k, char *inputName, char *filename) {
    FILE *f = fopen(filename, "wb");
    fprintf(f, "Algorithm: %s\n", name);
    fprintf(f, "Time: %lf sec\n", time);
    fprintf(f, "Input filename: %s\n", inputName);
    fprintf(f, "n=%d, d=%d, k=%d\n", n, d, k);
    fclose(f);
}

void runAndPresentResult(void (*runnable)(double *x, int n, int d, int k), double *x, int n, int d, int k, char *inputName, char *name, char *resultsFilename) {
    int SelfTID, NumTasks;
    MPI_Status mpistat;
    MPI_Request mpireq;
    MPI_Comm_size(MPI_COMM_WORLD, &NumTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &SelfTID);
    double time; // = run function

    if (SelfTID == 0) {
        double maxTime = time;
        for (int i = 1; i < NumTasks; i++) {
            double *incomingTime = malloc(sizeof(double));
            MPI_Recv(incomingTime, 1, MPI_DOUBLE, i, 101, MPI_COMM_WORLD, &mpistat);
            if (*incomingTime > maxTime)
                maxTime = *incomingTime;
            free(incomingTime);
        }
        saveStats(name, maxTime, n, d, k, inputName, resultsFilename);
    } else {
        MPI_Request mpireq;
        MPI_Isend(&time, 1, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD, &mpireq);
    }

}
