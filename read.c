#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utilities.h"
#include "controller.h"
#include "read.h"
#include "mmio.h"

void processMatrix(uint32_t *csc_rowOut, uint32_t *csc_colOut, int *n,char *file_path)
{

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    uint32_t i, *I, *J;
    double *val;

    if ((f = fopen(file_path, "r")) == NULL)
        exit(1);

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
        exit(1);

    /* reseve memory for matrices */

    I = (uint32_t *)malloc(nz * sizeof(uint32_t));
    J = (uint32_t *)malloc(nz * sizeof(uint32_t));
    val = (double *)malloc(nz * sizeof(double));

    for (i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d", &I[i], &J[i]);
        I[i]--; /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f != stdin)
        fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    //mm_write_banner(stdout, matcode);
    //mm_write_mtx_crd_size(stdout, M, N, nz);
    //Up to this point, I[] cointains row index and J[] column index for the nonzero elements
    /*for (i=0; i<nz; i++)
        {
             fprintf(stdout, "%d %d\n", I[i]+1, J[i]+1);
             printf("%d %d\n", I[i]+1, J[i]+1);
        }
		*/

    const uint32_t nnz = nz;

    printf("M is %d, nnz is %d\n", M, nnz);
    uint32_t *csc_row = (uint32_t *)malloc(nnz * sizeof(uint32_t));
    uint32_t *csc_col = (uint32_t *)malloc((M + 1) * sizeof(uint32_t));
    uint32_t isOneBased = 0;

    // Call coo2csc for isOneBase false
    coo2csc(csc_row, csc_col, I, J, nnz, M, isOneBased);

    *n = M; //from  mm_read_mtx_crd_size(f, &M, &N, &nz)
    csc_rowOut = csc_row;
    csc_colOut = csc_col;
}
