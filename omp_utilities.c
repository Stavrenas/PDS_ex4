#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "utilities.h"
#include "omp_utilities.h"
#include "mmio.h"
#include "read.h"

void multMatrixParallel(Matrix *A, Matrix *B, Matrix *C)
{
    //allocate memory for result mult
    uint32_t *c_elem, *elements, *c_idx, sizeA, sizeB;

    sizeA = A->size;
    sizeB = B->size;
    c_elem = (uint32_t *)malloc((sizeA + 1) * sizeof(uint32_t));
    elements = (uint32_t *)malloc((sizeA + 1) * sizeof(uint32_t)); //elements in each row

#pragma omp parallel
    {
        //allocate memory for local matrices used by each thread
        uint32_t *temp, indexB, indexA, last, localElements, totalElements, tempSize, start_a, start_b, end_a, end_b;
        temp = (uint32_t *)malloc(sizeof(uint32_t));
        tempSize = 1;
        last = -1;

        totalElements = 0;

        int nthreads = omp_get_num_threads();
        int id = omp_get_thread_num();

        for (uint32_t row = 1 + id; row <= sizeA; row += nthreads)
        { //go to each row of mtr A
            last = -1;
            localElements = 0; //elements in each row
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
                                continue;
                                //do not add it
                            }

                            else
                            {
                                temp[totalElements] = col;
                                totalElements++;
                                localElements++;

                                if (totalElements == tempSize)
                                { //realloc if needed
                                    tempSize++;
                                    uint32_t *tmp = realloc(temp, tempSize * sizeof(uint32_t));

                                    if (tmp != NULL)
                                        temp = tmp;
                                    else
                                    {
                                        tmp = malloc(tempSize * sizeof(uint32_t));
                                        for (int i = 0; i < tempSize - 1; i++)
                                            tmp[i] = temp[i];
                                        //free(temp);
                                        temp = tmp;
                                    }
                                }

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
#pragma omp single

        c_elem[0] = 0;
        for (uint32_t row = 1; row <= sizeA; row++)
            c_elem[row] = c_elem[row - 1] + elements[row]; //c_elem contains the total number of elements up to row[i]

        c_idx = malloc(c_elem[sizeA] * sizeof(uint32_t));

#pragma omp barrier //sync

        uint32_t start, end, index;
        index = 0;

        //each thread saves the indices in the final array
        for (uint32_t row = 1 + id; row <= sizeA; row += nthreads)
        {
            start = c_elem[row - 1];
            end = c_elem[row];
            for (uint32_t j = start; j < end; j++, index++)
                c_idx[j] = temp[index];
        }
    }
    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}

void multMatrixParallelMasked(Matrix *A, Matrix *B, Matrix *C, Matrix *mask)
{
    //allocate memory for result mult
    uint32_t *c_elem, *elements, *c_idx, sizeA, sizeB;

    sizeA = A->size;
    sizeB = B->size;
    c_elem = (uint32_t *)malloc((sizeA + 1) * sizeof(uint32_t));
    elements = (uint32_t *)malloc((sizeA + 1) * sizeof(uint32_t)); //elements in each row

#pragma omp parallel
    {
        //allocate memory for local matrices used by each thread
        uint32_t *temp, indexB, indexA, last, localElements, totalElements, tempSize, start_a, start_b, end_a, end_b;
        temp = (uint32_t *)malloc(sizeof(uint32_t));
        tempSize = 1;
        last = -1;

        totalElements = 0;

        int nthreads = omp_get_num_threads();
        int id = omp_get_thread_num();

        for (uint32_t row = 1 + id; row <= sizeA; row += nthreads)
        { //go to each row of mtr A
            last = -1;
            localElements = 0; //elements in each row

            uint32_t maskStart, maskEnd;
            maskStart = mask->csc_elem[row - 1];
            maskEnd = mask->csc_elem[row];

            for (uint32_t index = maskStart; index < maskEnd; index++)
            { //go to each col of mtr Mask
                uint32_t col = mask->csc_idx[index];
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
                                continue;
                                //do not add it
                            }

                            else
                            {
                                temp[totalElements] = col;
                                totalElements++;
                                localElements++;

                                if (totalElements == tempSize)
                                { //realloc if needed
                                    tempSize++;
                                    uint32_t *tmp = realloc(temp, tempSize * sizeof(uint32_t));

                                    if (tmp != NULL)
                                        temp = tmp;
                                    else
                                    {
                                        tmp = malloc(tempSize * sizeof(uint32_t));
                                        for (int i = 0; i < tempSize - 1; i++)
                                            tmp[i] = temp[i];
                                        //free(temp);
                                        temp = tmp;
                                    }
                                }

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
#pragma omp single

        c_elem[0] = 0;
        for (uint32_t row = 1; row <= sizeA; row++)
            c_elem[row] = c_elem[row - 1] + elements[row]; //c_elem contains the total number of elements up to row[i]

        c_idx = malloc(c_elem[sizeA] * sizeof(uint32_t));

#pragma omp barrier //sync

        uint32_t start, end, index;
        index = 0;

        //each thread saves the indices in the final array
        for (uint32_t row = 1 + id; row <= sizeA; row += nthreads)
        {
            start = c_elem[row - 1];
            end = c_elem[row];
            for (uint32_t j = start; j < end; j++, index++)
                c_idx[j] = temp[index];
        }
    }
    C->csc_idx = c_idx;
    C->csc_elem = c_elem;
    C->size = A->size;
}

void multBlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C)
{
    uint32_t size, totalBlocks, maxBlocks, blockSize;

    blockSize = A->list[0]->size;
    maxBlocks = floor(A->size / blockSize) + 1;
    if (A->size % blockSize == 0)
        maxBlocks--;
    size = 1;
    totalBlocks = 0;

    //initialize result matrix
    C->list = (Matrix **)malloc(size * sizeof(Matrix *));
    C->offsets = (uint32_t *)malloc(size * sizeof(uint32_t));
    C->row_ptr = (uint32_t *)malloc(maxBlocks * sizeof(uint32_t));

    for (uint32_t blockY = 1; blockY <= maxBlocks; blockY++)
    {
        C->row_ptr[blockY - 1] = totalBlocks;

        for (uint32_t blockX = 1; blockX <= maxBlocks; blockX++)
        {
            //Create block: Cp,q (p = BlockY, q = BlockX)

            Matrix *block = (Matrix *)malloc(sizeof(Matrix));
            Matrix *result = (Matrix *)malloc(sizeof(Matrix)); //used for mult

            //initialize block
            block->size = blockSize;
            block->csc_elem = (uint32_t *)malloc((blockSize + 1) * sizeof(uint32_t));
            block->csc_idx = (uint32_t *)malloc((0) * sizeof(uint32_t));

            for (int i = 0; i <= blockSize; i++)
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
                    //choose best algorithm for speeeeed
                    if (blockSize <= 40)
                        multMatrix2(A->list[indexA], B->list[indexB], result);
                    else
                        multMatrixParallel(A->list[indexA], B->list[indexB], result);

                    addMatrix(result, block, block);

                    //find block Bsq
                    for (int i = 1; i <= maxBlocks; i++)
                    {
                        indexB = findIndex(B, offsetB + maxBlocks * i);
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
                        indexB = findIndex(B, offsetB + maxBlocks * i);
                        if (indexB != -1)
                            break;
                    }
                }

                else if (sA < sB)
                    indexA++; //go to the next block in the same line of A
            }

            // if the mult results in a nonzero block, add it to the result matrix
            if (block->csc_elem[blockSize] != 0)
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
    C->blockSize = A->blockSize;
    C->totalBlocks = totalBlocks;
}

void multBlockedMatrixMasked(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C, BlockedMatrix *mask)
{
    uint32_t size, totalBlocks, maxBlocks, blockSize;

    blockSize = A->list[0]->size;
    maxBlocks = floor(A->size / blockSize) + 1;
    if (A->size % blockSize == 0)
        maxBlocks--;
    size = 1;
    totalBlocks = 0;

    //initialize result matrix
    C->list = (Matrix **)malloc(size * sizeof(Matrix *));
    C->offsets = (uint32_t *)malloc(size * sizeof(uint32_t));
    C->row_ptr = (uint32_t *)calloc(maxBlocks, sizeof(uint32_t));

    //go to every block of the mask matrix and calculate the corresponding Cpq
    for (uint32_t index = 0; index < mask->totalBlocks; index++)
    {
        uint32_t offset = mask->offsets[index];
        uint32_t blockX = (offset - 1) % maxBlocks + 1;
        uint32_t blockY = (offset - 1) / maxBlocks + 1;

        //Create block: Cp,q (p = BlockY, q = BlockX)
        Matrix *block = (Matrix *)malloc(sizeof(Matrix));
        Matrix *result = (Matrix *)malloc(sizeof(Matrix)); //used for mult

        //initialize block and result
        block->size = blockSize;
        block->csc_elem = (uint32_t *)calloc((blockSize + 1), sizeof(uint32_t));
        block->csc_idx = (uint32_t *)malloc((0) * sizeof(uint32_t));

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
            uint32_t sB = (offsetB - 1) / maxBlocks;

            if (sA == sB)
            {
                //choose best algorithm for speeeeed
                if (blockSize <= 40)
                    multMatrixMasked(A->list[indexA], B->list[indexB], result, mask->list[index]);
                else
                    multMatrixParallelMasked(A->list[indexA], B->list[indexB], result, mask->list[index]);

                addMatrix(result, block, block);

                //find block Bsq
                for (int i = 1; i <= maxBlocks; i++)
                {
                    indexB = findIndex(B, offsetB + maxBlocks * i);
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
                    indexB = findIndex(B, offsetB + maxBlocks * i);
                    if (indexB != -1)
                        break;
                }
            }

            else if (sA < sB)
                indexA++; //go to the next block in the same line of A
        }

        // if the mult results in a nonzero block, add it to the result matrix
        if (block->csc_elem[blockSize] != 0)
        {
            C->list[totalBlocks] = block;
            C->offsets[totalBlocks] = offset;
            totalBlocks++;
            if (size == totalBlocks)
            {
                size++;
                C->list = realloc(C->list, size * sizeof(Matrix *));
                C->offsets = realloc(C->offsets, size * sizeof(uint32_t *));
            }
        }

        C->row_ptr[blockY] = totalBlocks;
        free(result);
    }

    for (int i = 1; i < maxBlocks; i++)
    {
        if (C->row_ptr[i] == 0)
            C->row_ptr[i] = C->row_ptr[i - 1];
    }

    C->size = A->size;
    C->blockSize = A->blockSize;
    C->totalBlocks = totalBlocks;
}

void blockMatrix(Matrix *mtr, uint32_t blockSize, BlockedMatrix *blockedMatrix)
{
    uint32_t maxBlocks = floor(mtr->size / blockSize) + 1;
    if (mtr->size % blockSize == 0)
        maxBlocks--;
    uint32_t totalBlocks = 0;
    uint32_t listSize = 1; //also equals to offset size

    //initialize result matrix
    blockedMatrix->list = (Matrix **)malloc(1 * sizeof(Matrix *));
    blockedMatrix->offsets = (uint32_t *)malloc(1 * sizeof(uint32_t)); //maximum size of blocks
    blockedMatrix->row_ptr = (uint32_t *)malloc(maxBlocks * sizeof(uint32_t));

    if (blockSize > mtr->size / 2)
    {
        printf("Invalid blocksize\n");
        exit(-90);
    }

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
    blockedMatrix->blockSize = blockSize;
    blockedMatrix->totalBlocks = totalBlocks;
}

void unblockMatrix(BlockedMatrix *blockedMatrix, Matrix *mtr)
{
    mtr->size = blockedMatrix->size;
    mtr->csc_elem = (uint32_t *)calloc((mtr->size + 1), sizeof(uint32_t));
    //initialize csc_elem

    uint32_t blockSize = blockedMatrix->blockSize;
    uint32_t maxBlocks = floor(mtr->size / blockSize) + 1;
    if (mtr->size % blockSize == 0)
        maxBlocks--;

    uint32_t idx_size = 0;
    //add total block elements to find idx_size
    for (int i = 0; i < blockedMatrix->totalBlocks; i++)
        idx_size += blockedMatrix->list[i]->csc_elem[blockSize];

    mtr->csc_idx = (uint32_t *)malloc(idx_size * sizeof(uint32_t));

    uint32_t row_elements, block_row, block_start, block_end, row;
    uint32_t elements = 0;

    for (uint32_t mtr_row = 1; mtr_row <= mtr->size; mtr_row++)
    {
        block_row = (mtr_row - 1) / blockSize; //current row of the blocked matrix
        block_start = blockedMatrix->row_ptr[block_row];

        if (block_row + 1 == maxBlocks)
            block_end = blockedMatrix->totalBlocks;

        else
            block_end = blockedMatrix->row_ptr[block_row + 1];

        for (uint32_t block_idx = block_start; block_idx < block_end; block_idx++)
        {
            uint32_t upper, lower; //bounds for the corresponding block
            uint32_t blockX;       //horizontal block coordinate

            uint32_t offset = blockedMatrix->offsets[block_idx];

            blockX = (offset - 1) % maxBlocks + 1;
            row = (mtr_row - 1) % blockSize + 1; //current row of the block

            uint32_t start = blockedMatrix->list[block_idx]->csc_elem[row - 1];
            uint32_t end = blockedMatrix->list[block_idx]->csc_elem[row];

            for (int i = start; i < end; i++, elements++)
                mtr->csc_idx[elements] = blockedMatrix->list[block_idx]->csc_idx[i] + (blockX - 1) * blockSize; //compensate for the offset
            mtr->csc_elem[mtr_row] = elements;
        }
    }

    //update csc_elem matrix to the corresponding format
    for (int i = 1; i <= mtr->size; i++)
    {
        if (mtr->csc_elem[i] == 0)
            mtr->csc_elem[i] = mtr->csc_elem[i - 1];
    }
}

void addBlockedMatrix(BlockedMatrix *A, BlockedMatrix *B, BlockedMatrix *C)
{

    uint32_t blockSize = A->blockSize;
    uint32_t maxBlocks = floor(A->size / blockSize) + 1;
    if (A->size % blockSize == 0)
        maxBlocks--;
    uint32_t Blocks, totalBlocksA, totalBlocksB;
    totalBlocksA = A->totalBlocks;
    totalBlocksB = B->totalBlocks;
    Blocks = totalBlocksA + totalBlocksB;

    //initialize result matrix
    C->size = A->size;
    C->blockSize = blockSize;
    C->list = (Matrix **)malloc(Blocks * sizeof(Matrix *));
    C->offsets = (uint32_t *)malloc(Blocks * sizeof(uint32_t)); //maximum size of blocks
    C->row_ptr = (uint32_t *)malloc(maxBlocks * sizeof(uint32_t));
    uint32_t offsetA, offsetB;
    bool added = false;

    //merge block lists of A and B matrices
    for (int i = 0, iA = 0, iB = 0; iA != totalBlocksA && iB != totalBlocksB;)
    {
        added = false;
        offsetA = offsetB = 0;
        if (iA < totalBlocksA)
            offsetA = A->offsets[iA];
        if (iB < totalBlocksB)
            offsetB = B->offsets[iB];
        //printf("OffsetA is %d and offsetB is %d\n", offsetA, offsetB);
        if (offsetA > offsetB)
        {
            if (iB != totalBlocksB)
            {
                C->list[i] = B->list[iB];
                C->offsets[i] = offsetB;
                iB++;
                added = true;
            }
            else if (iA != totalBlocksA)
            {
                C->list[i] = A->list[iA];
                C->offsets[i] = offsetA;
                iA++;
                added = true;
            }
        }

        else if (offsetA < offsetB)
        {
            if (iA != totalBlocksA)
            {
                C->list[i] = A->list[iA];
                C->offsets[i] = offsetA;
                iA++;
                added = true;
            }
            else if (iB != totalBlocksB)
            {
                C->list[i] = B->list[iB];
                C->offsets[i] = offsetB;
                iB++;
                added = true;
            }
        }

        else if (offsetA == offsetB)
        {
            Matrix *result = (Matrix *)malloc(sizeof(Matrix)); //used for mult
            addMatrix(A->list[iA], B->list[iB], result);

            C->list[i] = result;
            C->offsets[i] = offsetA;

            if (iA != totalBlocksA)
                iA++;
            if (iB != totalBlocksB)
                iB++;

            added = true;
        }

        if (added)
        {
            uint32_t row = (C->offsets[i] - 1) / maxBlocks + 1;
            //printf("row is %d and i is %d\n",row,i);
            i++;
            C->row_ptr[row] = i;
            C->totalBlocks = i;
        }
    }
}
