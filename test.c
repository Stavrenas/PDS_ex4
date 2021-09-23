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

    printf("\n");
    char filename[] = "12.mtx";
    uint32_t *csc_row, *csc_col;
    int n;

    readMatrix(csc_row, csc_col, &n, filename);

    printf("TEST is %d\n", n);

    for (uint32_t i = 1; i <= n; i++)
    {
        for (uint32_t i1 = csc_col[i - 1]; i1 < csc_col[i]; i1++)
        {
            printf("%d ",csc_row[i1]);
        }
        printf("\n");
    }

  }
