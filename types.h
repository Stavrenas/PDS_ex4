#include <stdint.h>

#ifndef TYPES_H
#define TYPES_H

typedef struct
{
    uint32_t *csc_idx; //column indexes
    uint32_t *csc_elem;
    uint32_t size;
} Matrix;

typedef struct
{
    uint32_t size;
    uint32_t totalBlocks;
    uint32_t *offsets;
    Matrix **list;
} BlockedMatrix;

#endif //TYPES_H
