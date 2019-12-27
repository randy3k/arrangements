#ifndef PERMUTATIONS_UTILS_H__
#define PERMUTATIONS_UTILS_H__

#include <stdlib.h>

static void swap(unsigned int *ar, unsigned int first, unsigned int second)
{
    unsigned int temp = ar[first];
    ar[first] = ar[second];
    ar[second] = temp;
}

static void reverse(unsigned int *ar, size_t len)
{
    unsigned int i, j;

    for (i = 0, j = len - 1; i < j; i++, j--) {
        swap(ar, i, j);
    }
}


#endif /* end of include guard: PERMUTATIONS_UTILS_H__ */
