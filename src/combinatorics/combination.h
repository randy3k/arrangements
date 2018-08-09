#ifndef _COMB_H
#define _COMB_H 1

#include <stddef.h>

// mirror of python itertools.combinations
// https://docs.python.org/3/library/itertools.html#itertools.combinations
unsigned int next_combination(unsigned int *ar, size_t n, unsigned int k)
{
    // ar = [0, 1, ..., r-1]
    unsigned int i, j, temp;

    for (i = k-1; ; i--) {
        if (ar[i] != i + n - k) {
            break;
        }
        else if (i == 0) {
            return 0;
        }
    }

    // turn the elements after it into a linear sequence
    temp = ar[i];
    for (j = i; j < k; j++) {
        temp += 1;
        ar[j] = temp;
    }
    return 1;
}

#endif
