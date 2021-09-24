#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h> // sqrt
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "utilities.h"
#include "mmio.h"


void swap(double *n1, double *n2)
{
    double temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void swapInts(int *n1, int *n2)
{
    int temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void dividePoints(int n, int tasks, int *array)
{
    int points = n / tasks;
    for (int i = 0; i < n; i++)
    {
        array[i] = points;
    }

    int pointsLeft = n % tasks;
    for (int i = 0; pointsLeft > 0; i++)
    {
        array[i]++;
        pointsLeft--;
    }
}

int findDestination(int id, int NumTasks)
{
    if (id == NumTasks - 1)
        return 0;
    else
        return (id + 1);
}

int findSender(int id, int NumTasks)
{
    if (id == 0)
        return (NumTasks - 1);
    else
        return (id - 1);
}

struct timeval tic()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv;
}

double toc(struct timeval begin)
{
    struct timeval end;
    gettimeofday(&end, NULL);
    double stime = ((double)(end.tv_sec - begin.tv_sec) * 1000) +
                   ((double)(end.tv_usec - begin.tv_usec) / 1000);
    stime = stime / 1000;
    return (stime);
}

void multMatrix2(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, *temp, indexB, indexA, last;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    temp = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    uint32_t idx_size = 1;
    uint32_t end = 0;

    c_elem[0] = 0;
    last = -1;
    elements = 0;

    int perc = 0;
    uint32_t sizeA, sizeB, start_a, end_a, start_b, end_b;

    sizeA = A->size;
    sizeB = B->size;
    for (uint32_t row = 1; row <= sizeA; row++)
    { //go to each row of mtr A

        // if (row * 100 / A->size != perc)
        // {
        //     printf("%d%% \n", perc);
        //     perc = row * 100 / A->size;
        // }
        //printf("Row %d of %d\n",row,A->size);

        for (uint32_t col = 1; col <= sizeB; col++)
        { //go to each col of mtr Β

            start_a = A->csc_elem[row - 1];
            end_a = A->csc_elem[row];
            for (uint32_t a = start_a; a < end_a; a++)
            { //go to each element in row "row" of mtr A

                indexA = A->csc_idx[a];
                start_b = B->csc_elem[indexA - 1];
                end_b = B->csc_elem[indexA];

                for (uint32_t b = start_b; b < end_b; b++)
                { //check if there is a match in col "col"of mtr B
                    indexB = B->csc_idx[b];

                    if (indexB > col)
                        break;

                    else if (indexB == col)
                    {
                        if (last == col) //check if the element is already added
                        {
                            //do not add it
                        }

                        else
                        {
                            c_idx[elements] = col;
                            elements++;
                            if (elements == idx_size)
                            {
                                idx_size *= 2;
                                c_idx = realloc(c_idx, idx_size * sizeof(uint32_t));
                            }

                            last = col;
                        }
                        break;
                    }
                }
            }
        }

        c_elem[row] = elements;
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}

void multMatrix(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, *temp, indexB, last;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    temp = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    uint32_t idx_size = 0;
    uint32_t end = 0;

    c_elem[0] = 0;
    last = -1; //the first element is set to -1 and is used to avoid adding the same element twice

    for (uint32_t row = 1; row <= A->size; row++)
    { //go to each row of mtr A

        elements = 0;
        c_elem[row] = c_elem[row - 1] + elements;

        for (uint32_t col = 1; col <= B->size; col++)
        { //go to each col of mtr Β

            for (uint32_t a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
            { //go to each element in row "row" of mtr A

                uint32_t indexA = A->csc_idx[a];

                for (uint32_t b = B->csc_elem[indexA - 1]; b < B->csc_elem[indexA]; b++)
                { //check if there is a match in col "col"of mtr B
                    indexB = B->csc_idx[b];

                    if (indexB > col)
                        break;

                    else if (indexB == col)
                    {
                        if (last == col) //check if the element is already added
                        {
                            continue;
                            //do not add it
                        }

                        else
                        {
                            temp[elements] = col;
                            last = col;
                            elements++;
                        }
                        break;
                    }
                }
            }
        }

        c_elem[row] = c_elem[row - 1] + elements;

        if (elements != 0)
        {
            idx_size += elements * sizeof(uint32_t); //in bytes
            c_idx = realloc(c_idx, idx_size);
            end += elements;

            for (uint32_t i = end - elements; i < end; i++)
                c_idx[i] = temp[i - end + elements];
        }
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}

void multMatrixParallel(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *elements, *c_idx;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    elements = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    c_elem[0] = 0;
    int perc = 0;

#pragma omp parallel
    {
        uint32_t *temp, indexB, indexA, last, localElements, totalElements, tempSize;
        temp = (uint32_t *)malloc(sizeof(uint32_t));
        tempSize = 1;
        last = -1;

        localElements = 0;
        totalElements = 0;

        int nthreads = omp_get_num_threads();
        int id = omp_get_thread_num();

        for (uint32_t row = 1 + id; row <= A->size; row += nthreads)
        { //go to each row of mtr A

            localElements = 0;

            for (uint32_t col = 1; col <= B->size; col++)
            { //go to each col of mtr Β

                for (uint32_t a = A->csc_elem[row - 1]; a < A->csc_elem[row]; a++)
                { //go to each element in row "row" of mtr A

                    indexA = A->csc_idx[a];

                    for (uint32_t b = B->csc_elem[indexA - 1]; b < B->csc_elem[indexA]; b++)
                    { //check if there is a match in col "col"of mtr B
                        indexB = B->csc_idx[b];

                        if (indexB > col)
                            break;

                        else if (indexB == col)
                        {
                            if (last == col) //check if the element is already added
                            {
                                //do not add it
                            }

                            else
                            {
                                temp[totalElements] = col; //replaced local elements with total elements
                                totalElements++;
                                if (totalElements == tempSize)
                                {
                                    tempSize *= 2;
                                    temp = realloc(temp, tempSize * sizeof(uint32_t));
                                }
                                localElements++;
                                last = col;
                            }
                            break;
                        }
                    }
                }
            }

            elements[row] = localElements;
        }

#pragma omp barrier //sync

        if (id == 0)
        {
            for (uint32_t row = 1; row <= A->size; row++)
                c_elem[row] = c_elem[row - 1] + elements[row];
            c_idx = realloc(c_idx, c_elem[A->size] * sizeof(uint32_t));
        }

#pragma omp barrier //sync

        uint32_t index = 0;
        uint32_t start, end;
        for (uint32_t row = 1 + id; row <= A->size; row += nthreads)
        {

            start = c_elem[row - 1];
            end = c_elem[row];

            for (uint32_t j = start; j < end; j++)
                c_idx[j] = temp[j - start];
        }
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}

void blockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix)
{
    uint32_t maxBlocks = ceil(mtr->size / blockSize);
    uint32_t totalBlocks = 0;
    uint32_t listSize = 1; //also equals to offset size

    blockedMatrix->list = (Matrix **)malloc(1 * sizeof(Matrix *));
    blockedMatrix->offsets = (uint32_t *)malloc(1 * sizeof(uint32_t)); //maximum size of blocks
    blockedMatrix->row_ptr = (uint32_t *)malloc(maxBlocks * sizeof(uint32_t));

    for (uint32_t blockY = 1; blockY <= maxBlocks; blockY++)
    {
        // Save first block of current row (blockY)
        blockedMatrix->row_ptr[blockY-1] = totalBlocks;

        for (uint32_t blockX = 1; blockX <= maxBlocks; blockX++)
        {

            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            int *block_idx, *block_elem, elements, idx_size;

            block_idx = (int *)malloc(sizeof(int));
            block_elem = (int *)malloc((blockSize + 1) * sizeof(int));

            idx_size = 1;
            block_elem[0] = 0;
            elements = 0;

            for (int row = (blockY - 1) * blockSize + 1; row < blockY * blockSize + 1; row++)
            {
                // Check if row exceeds matrix size
                if (row <= mtr->size)
                {
                    uint32_t start = mtr->csc_elem[row - 1];
                    uint32_t end = mtr->csc_elem[row];

                    for (int j = start; j < end; j++)
                    {

                        if (mtr->csc_idx[j] > blockX * blockSize) //check if it is worth it
                            break;

                        else if (mtr->csc_idx[j] <= blockX * blockSize && mtr->csc_idx[j] > (blockX - 1) * blockSize)
                        {

                            block_idx[elements] = mtr->csc_idx[j] - (blockX - 1) * blockSize;
                            elements++;

                            if (elements == idx_size)
                            {
                                idx_size++;
                                block_idx = realloc(block_idx, idx_size * sizeof(int));
                            }
                        }
                    }

                    block_elem[row - (blockY - 1) * blockSize] = elements;
                }

                else // zero padding vertically
                    block_elem[row - (blockY - 1) * blockSize] = block_elem[row - (blockY - 1) * blockSize - 1];
            }

            if (elements != 0)
            {
                block->size = blockSize;
                block->csc_idx = block_idx;
                block->csc_elem = block_elem;

                blockedMatrix->list[totalBlocks] = block;
                blockedMatrix->offsets[totalBlocks] = (blockY - 1) * maxBlocks + blockX;

                totalBlocks++;

                if (listSize == totalBlocks)
                {
                    listSize++;
                    blockedMatrix->list = realloc(blockedMatrix->list, listSize * sizeof(Matrix *));
                    blockedMatrix->offsets = realloc(blockedMatrix->offsets, listSize * sizeof(uint32_t *));
                }
            }
        }
    }

    blockedMatrix->size = mtr->size;
    blockedMatrix->totalBlocks = totalBlocks;

    printf("Max blocks are %d and current blocks: %d. Non zero blocks: %f \n", maxBlocks * maxBlocks, totalBlocks, (float)(totalBlocks) / (maxBlocks * maxBlocks));
    printf("Start of each row:\n");
    for (int i = 0; i < maxBlocks; ++i)
    {
        printf("row: %d, block: %d\n", i+1, blockedMatrix->row_ptr[i]);
    }
}

void unblockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix)
{
    // Function to convert blocked matrix into unblocked csc

    mtr->size = blockedMatrix->size;
    // Allocate memory for unblocked matrix
    mtr->csc_idx = (uint32_t *)malloc(sizeof(uint32_t));
    mtr->csc_elem = (uint32_t *)malloc((mtr->size + 1) * sizeof(uint32_t));

    mtr->csc_elem[0] = 0;

    uint32_t maxBlocks = ceil(mtr->size / blockSize);
    // Index to iterate over all blocks
    uint32_t currentBlock = 0; // blockedMatrix->row_ptr[0]

    // Indices showing current row's column indices inside csc_idx
    int row_start = 0;
    int row_end = 0;

    // Total column indices added to mtr->csc_idx
    int elements = 0;
    // Total elements allocated for mtr->csc_idx (more than actual)
    int idx_size = 1;

    // Block column offset
    int blockX = 0;

    // Construct each row of the unblocked matrix
    for (int row = 0; row < mtr->size; row++)
    {
        // Loop through blocks containing current row
        while (blockedMatrix->offsets[currentBlock] >= floor(row/blockSize) + 1 &&
        blockedMatrix->offsets[currentBlock] <= floor(row/blockSize) + maxBlocks)
        {
            // Get column indices of current block's rowls
            row_start = blockedMatrix->list[currentBlock]->csc_elem[row % blockSize];
            row_end = blockedMatrix->list[currentBlock]->csc_elem[(row + 1) % blockSize];

            for (int i = row_start; i < row_end; i++)
            {
                // Find block column offset (blockX)
                blockX = (blockedMatrix->offsets[currentBlock] - 1) % blockSize + 1;
                // Save column index taking blockX offset into account
                mtr->csc_idx[elements] = blockedMatrix->list[currentBlock]->csc_idx[i] + (blockX-1)*blockSize;
                elements++;

                // Reallocate memory if needed
                if (elements == idx_size)
                {
                    idx_size++;
                    mtr->csc_idx = realloc(mtr->csc_idx, idx_size * sizeof(int));
                }
            }

            // Go to next block containing current row
            currentBlock++;
        }

        // Check if current 'block-row' contains next row of mtr
        // if so then iterate the same blocks
        if(ceil((row+1)/blockSize) - ceil(row/blockSize) < 1)
        {
            currentBlock = blockedMatrix->row_ptr[row/blockSize];
        } else // Go to the first block of the next 'block-row'
        {
            currentBlock = blockedMatrix->row_ptr[(row+1)/blockSize];
        }

        // Save non zero elements of current row
        mtr->csc_elem[row] = elements;
    }

    // Set last element
    mtr->csc_elem[mtr->size] = mtr->size;

    printf("%d\n", elements);
}

void addMatrix(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, idx_size;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));

    idx_size = 1;
    c_elem[0] = 0;
    elements = 0;

    uint32_t size, start_a, end_a, start_b, end_b, indexA, indexB;

    size = A->size;

    for (uint32_t row = 1; row <= size; row++)
    { //go to each row of mtr A
        //printf("row %d\n",row);

        start_a = A->csc_elem[row - 1];
        end_a = A->csc_elem[row];

        start_b = B->csc_elem[row - 1];
        end_b = B->csc_elem[row];
        
        //printf("start_a is %d , start_b is %d, end_a is %d and end_b is %d\n",start_a,start_b,end_a,end_b);
        for (uint32_t a = start_a, b = start_b;;)
        { //go to each element in row of mtr A and mtr B

            if (a == end_a && b == end_b)
                break;
            else
            {
                indexA = A->csc_idx[a];
                indexB = B->csc_idx[b];

                //printf("a is %d and b is %d, indexA is %d and indexB is %d\n", a, b, indexA, indexB);

                if (indexA < indexB)
                {
                    if (a != end_a)
                    {
                        c_idx[elements] = indexA;
                        elements++;
                        a++;
                    }
                    else if (b != end_b)
                    {
                        c_idx[elements] = indexB;
                        elements++;
                        b++;
                    }
                }

                else if (indexA > indexB)
                {
                    if (b != end_b)
                    {
                        c_idx[elements] = indexB;
                        elements++;
                        b++;
                    }
                    else if (a != end_a)
                    {
                        c_idx[elements] = indexA;
                        elements++;
                        a++;
                    }
                }

                else
                {
                    c_idx[elements] = indexB;
                    elements++;
                    if (a != end_a)
                        a++;
                    if (b != end_b)
                        b++;
                }

                if (elements == idx_size)
                {
                    idx_size *= 2;
                    c_idx = realloc(c_idx, idx_size * sizeof(uint32_t));
                }
            }
        }

        c_elem[row] = elements;
    }

    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}

int findIndex(BlockedMatrix *mtr, uint32_t indx)
{
    return binarySearch(mtr->offsets, 0, mtr->totalBlocks - 1, indx);
}

int binarySearch(uint32_t *list, uint32_t left, uint32_t right, uint32_t index)
{
    if (right >= left)
    {

        int mid = left + (right - left) / 2;
        // If the element is present at the middle
        // itself
        if (list[mid] == index)
            return mid;

        // If element is smaller than mid, then
        // it can only be present in left subarray
        if (list[mid] > index)
            return binarySearch(list, left, mid - 1, index);

        // Else the element can only be present
        // in right subarray
        return binarySearch(list, mid + 1, right, index);
    }

    // We reach here when element is not
    // present in array
    return -1;
}

void multBlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C)
{
    uint32_t size, totalBlocks;

    uint32_t blockSize = A->list[0]->size;
    uint32_t maxBlocks = ceil(A->size / blockSize);

    size = 1;
    totalBlocks = 0; // LOL LOL KILL ME

    //initialize result matrix
    C->list = (Matrix **)malloc(size * sizeof(Matrix *));
    C->offsets = (uint32_t *)malloc(size * sizeof(uint32_t));

    for (uint32_t blockY = 1; blockY <= maxBlocks; blockY++)
    {
        for (uint32_t blockX = 1; blockX <= maxBlocks; blockX++)
        {
            //Create block: Cp,q (p = BlockY, q = BlockX)

            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            Matrix *result = (Matrix *)malloc(sizeof(Matrix)); //used for mult

            //initialize block
            block->csc_elem = (uint32_t *)malloc((block->size + 1) * sizeof(uint32_t));
            block->csc_idx = (uint32_t *)malloc((0) * sizeof(uint32_t));
            block->size = A->list[0]->size; //get blocksize

            for (int i = 0; i <= block->size; i++)
                block->csc_elem[i] = 0;


            //find indexes of Ap1 and B1q
            uint32_t indexA, indexB;
            for (int i = 1; i <= maxBlocks; i++)
            {
                indexA = findIndex(A, (blockY - 1) * maxBlocks + i);
                if (indexA != -1)
                    break; //stop when we find the first nonzero block in row p
            }

            for (int i = 1; i <= maxBlocks; i++)
            {
                indexB = findIndex(B, maxBlocks * (i - 1) + blockX);
                if (indexB != -1)
                    break; //stop when we find the first nonzero block in column q
            }

            //maxBlocks is the maximum number of mults for Cp,q. Variable s is not used
            for (int s = 1; s <= maxBlocks; s++)
            {

                //if either block does not exist
                if (indexA == -1 || indexB == -1)
                    break;

                uint32_t offsetA = A->offsets[indexA];
                uint32_t offsetB = B->offsets[indexB];

                //break if blocks are out of the desired row/col
                if (offsetA > blockY * maxBlocks || offsetB % maxBlocks > blockX)
                    break;
                //or if we run out of blocks
                if (indexB > B->totalBlocks || indexA > A->totalBlocks)
                    break;

                //check if the blocks match
                uint32_t sA = (offsetA - 1) % maxBlocks;
                uint32_t sB = floor((offsetB - 1) / maxBlocks);

                if (sA == sB)
                {
                    multMatrix2(A->list[indexA], B->list[indexB], result);
                    addMatrix(result, block, block);

                    //find block Bsq
                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks);
                        if (indexB != -1)
                            break; 
                    }

                    indexA++;
                }

                else if (sA > sB)
                {
                    //find block Bsq
                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks);
                        if (indexB != -1)
                            break; 
                    }
                }

                else if (sA < sB)
                    indexA++; //go to the next block in the same line of A
            }


            // if the mult results in a nonzero block, add it to the result matrix
            if (block->csc_elem[block->size] != 0)
            {
                C->list[totalBlocks] = block;
                C->offsets[totalBlocks] = (blockY - 1) * maxBlocks + blockX;
                totalBlocks++;

                if (size == totalBlocks)
                {
                    size++;
                    C->list = realloc(C->list, size * sizeof(Matrix *));
                    C->offsets = realloc(C->offsets, size * sizeof(uint32_t *));
                }
            }
            free(result); //result matrix will not be needed in the future, counter to "block" matrix
        }
    }
    C->size = A->size;
    C->totalBlocks = totalBlocks;
}