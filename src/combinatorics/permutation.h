#ifndef _PERM_H
#define _PERM_H 1

#include <stddef.h>
#include "utils.h"


unsigned int next_permutation(unsigned int *ar, size_t n)
{
    unsigned int k, j;
    unsigned int result = 0;

    // trival for only one element
    if (n == 1) {
        return result;
    }
    // find the largest k such that a[k] < a[k + 1]
    for (k = n - 1; k && ar[k - 1] >= ar[k]; k--);

    if (!k--) {
        // if not found, array is highest permutation
        reverse(ar, n);
    } else {
        // fnd the largest index j such that a[k] < a[l]
        for (j = n - 1; ar[j] <= ar[k]; j--);
        // swap it with the first one
        swap(ar, k, j);
        // reverse the remainder
        reverse(ar + k + 1, n - k - 1);
        result = 1;
    }
    return result;
}

#endif
