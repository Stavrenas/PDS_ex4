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
    //allocate memory for result mult
    uint32_t *c_elem, *c_idx, elements, *temp, indexB, indexA, last;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    temp = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    uint32_t idx_size = 1;
    uint32_t end = 0;

    c_elem[0] = 0;
    last = -1;
    elements = 0;

    //variables used in for loops
    uint32_t sizeA, sizeB, start_a, end_a, start_b, end_b;

    sizeA = A->size;
    sizeB = B->size;

    for (uint32_t row = 1; row <= sizeA; row++)
    { //go to each row of mtr A

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
                { //check if there is a match in col "col" of mtr B
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
                            { //faster realloc than multMatrix
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
    //alocate memory for result matrix
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
            { //go to each element in "row"-th row  of mtr A

                uint32_t indexA = A->csc_idx[a];

                for (uint32_t b = B->csc_elem[indexA - 1]; b < B->csc_elem[indexA]; b++)
                { //check if there is a match in col "col" of mtr B
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
        { //add the indexes and realloc

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

    //allocate memory for result mult
    uint32_t *c_elem, *elements, *c_idx;

    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    elements = (uint32_t *)malloc((A->size) * sizeof(uint32_t));

    c_elem[0] = 0;
    int perc = 0;

#pragma omp parallel
    {
        //allocate memory for local matrices used by each thread
        uint32_t *temp, indexB, indexA, last, localElements, totalElements, tempSize;
        temp = (uint32_t *)malloc(sizeof(uint32_t));
        tempSize = 1;
        last = -1;

        localElements = 0; //elements for each row
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
                            { //add it

                                temp[totalElements] = col;
                                totalElements++;
                                if (totalElements == tempSize)
                                { //realloc if needed
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

            elements[row] = localElements; //this matrix contains the number of elements in each row
        }

#pragma omp barrier //sync

        if (id == 0)
        {
            for (uint32_t row = 1; row <= A->size; row++)
                c_elem[row] = c_elem[row - 1] + elements[row]; //c_elem contains the total number of elements up to row[i]
            c_idx = realloc(c_idx, c_elem[A->size] * sizeof(uint32_t));
        }

#pragma omp barrier //sync

        uint32_t index = 0;
        uint32_t start, end;

        //each thread saves the indices in the final array
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

    //initialize result matrix
    blockedMatrix->list = (Matrix **)malloc(1 * sizeof(Matrix *));
    blockedMatrix->offsets = (uint32_t *)malloc(1 * sizeof(uint32_t)); //maximum size of blocks
    blockedMatrix->row_ptr = (uint32_t *)malloc(maxBlocks * sizeof(uint32_t));

    for (uint32_t blockY = 1; blockY <= maxBlocks; blockY++)
    {
        // Save first block of current row (blockY)
        blockedMatrix->row_ptr[blockY - 1] = totalBlocks;

        for (uint32_t blockX = 1; blockX <= maxBlocks; blockX++)
        {

            //initialize block to be created
            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            int *block_idx, *block_elem, elements, idx_size;

            block_idx = (int *)malloc(sizeof(int));
            block_elem = (int *)malloc((blockSize + 1) * sizeof(int));

            idx_size = 1;
            block_elem[0] = 0;
            elements = 0;

            for (int row = (blockY - 1) * blockSize + 1; row < blockY * blockSize + 1; row++) //calculate row from vertical and horizontal coordinates
            {
                // Check if row exceeds matrix size
                if (row <= mtr->size)
                {
                    uint32_t start = mtr->csc_elem[row - 1];
                    uint32_t end = mtr->csc_elem[row];

                    for (int j = start; j < end; j++)
                    {
                        //break if we exceed block boundaries
                        if (mtr->csc_idx[j] > blockX * blockSize) //check if it is worth it
                            break;

                        //check if the element is inside the block
                        else if (mtr->csc_idx[j] <= blockX * blockSize && mtr->csc_idx[j] > (blockX - 1) * blockSize)
                        { //add it

                            block_idx[elements] = mtr->csc_idx[j] - (blockX - 1) * blockSize;
                            elements++;

                            if (elements == idx_size)
                            {
                                idx_size++;
                                block_idx = realloc(block_idx, idx_size * sizeof(int));
                            }
                        }
                    }
                    //update csc_elem
                    block_elem[row - (blockY - 1) * blockSize] = elements;
                }

                else // if a row of a block is outside the matrix, we assume it contains only zero elements so we only update csc_elem
                    block_elem[row - (blockY - 1) * blockSize] = block_elem[row - (blockY - 1) * blockSize - 1];
            }

            if (elements != 0)
            { //create and add block to the list

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

    // printf("Start of each row:\n");
    // for (int i = 0; i < maxBlocks; ++i)
    //     printf("row: %d, block: %d\n", i + 1, blockedMatrix->row_ptr[i]);
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
    for (int row = 1; row <= mtr->size; row++)
    {
        // Loop through blocks containing current row
        // check offset to see if a block contains current row
        while (blockedMatrix->offsets[currentBlock] >= (floor((row-1) / blockSize) * maxBlocks + 1) &&
               blockedMatrix->offsets[currentBlock] <= (floor((row-1) / blockSize) * maxBlocks + maxBlocks))
        {
            // Get column indices of current block's rowls
            row_start = blockedMatrix->list[currentBlock]->csc_elem[(row-1) % blockSize];
            row_end = blockedMatrix->list[currentBlock]->csc_elem[(row) % blockSize];

            for (int i = row_start; i < row_end; i++)
            {
                // Find block column offset (blockX)
                blockX = (blockedMatrix->offsets[currentBlock] - 1) % blockSize + 1;
                // Save column index taking blockX offset into account
                mtr->csc_idx[elements] = blockedMatrix->list[currentBlock]->csc_idx[i] + (blockX - 1) * blockSize;
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
        if (ceil(row  / blockSize) - ceil((row-1) / blockSize) < 1)
        {
            currentBlock = blockedMatrix->row_ptr[row / blockSize];
        }
        else // Go to the first block of the next 'block-row'
        {
            currentBlock = blockedMatrix->row_ptr[(row + 1) / blockSize];
        }

        // Save non zero elements of current row
        mtr->csc_elem[row] = elements;
    }

    // Set last element
    mtr->csc_elem[mtr->size] = mtr->size;
}

void addMatrix(Matrix *A, Matrix *B, Matrix *C)
{
    uint32_t *c_elem, *c_idx, elements, idx_size;

    //initialize result mult
    c_idx = (uint32_t *)malloc(sizeof(uint32_t));
    c_elem = (uint32_t *)malloc((A->size + 1) * sizeof(uint32_t));
    idx_size = 1; //used in realloc
    c_elem[0] = 0;
    elements = 0;

    //variables used in for loops
    uint32_t size, start_a, end_a, start_b, end_b, indexA, indexB;
    size = A->size;

    for (uint32_t row = 1; row <= size; row++)
    { //go to each row of mtr C

        start_a = A->csc_elem[row - 1];
        end_a = A->csc_elem[row];

        start_b = B->csc_elem[row - 1];
        end_b = B->csc_elem[row];

        for (uint32_t a = start_a, b = start_b;;)
        { //  merge  current row of A and  B

            if (a == end_a && b == end_b)
                break;
            else
            {
                indexA = A->csc_idx[a];
                indexB = B->csc_idx[b];

                if (indexA < indexB)
                {
                    //check if row has elements
                    if (a != end_a)
                    {
                        c_idx[elements] = indexA;
                        elements++;
                        a++;
                    }
                    //check if row has elements
                    else if (b != end_b)
                    {
                        c_idx[elements] = indexB;
                        elements++;
                        b++;
                    }
                }

                else if (indexA > indexB)
                {
                    //check if row has elements
                    if (b != end_b)
                    {
                        c_idx[elements] = indexB;
                        elements++;
                        b++;
                    }
                    //check if row has elements
                    else if (a != end_a)
                    {
                        c_idx[elements] = indexA;
                        elements++;
                        a++;
                    }
                }

                else
                { //if indexA = indexB
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
    //check if a block with offset = indx exists
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
    totalBlocks = 0;

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

                    indexA++; //go to the next block in the same line of A
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
