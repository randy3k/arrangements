#include "k_permutation.h"

static void swap(unsigned int *ar, unsigned int first, unsigned int second)
{
    unsigned int temp = ar[first];
    ar[first] = ar[second];
    ar[second] = temp;
}

// this is the algorithm used in python itertools.permutations
unsigned int next_k_permutation(unsigned int *ar, unsigned int *cycle, size_t n, size_t k)
{
    long i, j;
    unsigned int temp;
    for (i = k-1; i >= 0; i--) {
        cycle[i] -= 1;
        if (cycle[i] == 0) {
            temp = ar[i];
            for (j = i; j <= n - 2; j++) {
                ar[j] = ar[j + 1];
            }
            ar[n - 1] = temp;
            cycle[i] = n - i;
        } else {
            j = cycle[i];
            swap(ar, i, n - j);
            return 1;
        }
    }
    return 0;
}
