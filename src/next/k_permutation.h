#ifndef _K_PERM_H
#define _K_PERM_H 1

#include <stddef.h>
#include "utils.h"

// mirror of python itertools.permutations
// https://docs.python.org/3/library/itertools.html#itertools.permutations
unsigned int next_k_permutation(unsigned int *ar, unsigned int *cycle, size_t n, size_t k)
{
    // ar = [0, 1, ..., n], cycle = [n, n-1, ..., n-r+1]
    unsigned int i, j, temp;

    for (i = k-1; ; i--) {
        cycle[i] -= 1;
        if (cycle[i] == 0) {
            temp = ar[i];
            for (j = i; j <= n - 2; j++) {
                ar[j] = ar[j + 1];
            }
            ar[n - 1] = temp;
            cycle[i] = n - i;
        } else {
            swap(ar, i, n - cycle[i]);
            return 1;
        }
        if (i == 0) {
            return 0;
        }
    }
}

#endif
