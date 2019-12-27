#ifndef PERMUTATIONS_UTILS_H__
#define PERMUTATIONS_UTILS_H__

#include <stdlib.h>

static double fact(int n) {
    double out;
    int i;
    out = 1;
    for(i=0; i<n; i++) {
        out = out * (n - i);
    }
    return out;
}

static double fallfact(int n, int k) {
    double out;
    int i;
    if (n < k) {
        return 0;
    }
    out = 1;
    for(i=0; i<k; i++) {
        out = out * (n - i);
    }
    return out;
}


static double multichoose(int* freq, size_t flen) {
    double out = 1;
    size_t h, i, j;
    h = 0;
    for (i=0; i<flen; i++) {
        for (j=1; j<=freq[i]; j++) {
            h++;
            out = out * h / j;
        }
    }
    return out;
}


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
