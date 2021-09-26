#include <stdint.h>

#ifndef TYPES_H
#define TYPES_H

typedef struct
{
    uint32_t *csc_idx;  // column indexes
    uint32_t *csc_elem; // cumulative elements of each row
    uint32_t size;      // Matrix size (assuming square matrices only)
} Matrix;

typedef struct
{
    uint32_t size;        // Matrix size (assuming square matrices only)
    uint32_t blockSize;   // Block size (assuming square blocks only)
    uint32_t totalBlocks; // Number of blocks (with non-zero elements)
    uint32_t *offsets;    // Block offset with relative to first block
    uint32_t *row_ptr;    // Index showing first block of each row
    Matrix **list;        // List containing blocks (with non-zero elements)
} BlockedMatrix;

#endif //TYPES_H
