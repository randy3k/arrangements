#ifndef _K_COMPOSITION_H
#define _K_COMPOSITION_H 1

#include <stddef.h>

unsigned int next_asc_k_composition(unsigned int *ar, size_t n, unsigned int k) {
    int i, j;
    unsigned int a;
    for (i = k-1; i>= 0; i--) {
        if (ar[i] > 1) {
            a = ar[i];
            break;
        }
    }
    if (i == 0) return 0;
    ar[i - 1] = ar[i - 1] + 1;

    for (j = i; j < k -1; j++) {
        ar[j] = 1;
    }
    ar[k - 1] = a - 1;
    return 1;
}


unsigned int next_desc_k_composition(unsigned int *ar, size_t n, unsigned int k, int* tp) {
    int i, j;
    int t = *tp;

    for (i = k-1; i>= 1; i--) {
        if (ar[i - 1] > 1) {
            break;
        }
    }
    if (i == 0) {
        return 0;
    }
    if (t > 0) {
        t--;
    } else {
        t = 0;
        for (j = i; j < k - 1; j++) {
            t += ar[j];
        }
    }
    ar[i - 1] = ar[i - 1] - 1;
    ar[i] = t + ar[k - 1] - (k - i - 1) + 1;
    for (j = i+1; j < k; j++) ar[j] = 1;

    *tp = t;
    return 1;
}

#endif
