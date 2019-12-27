#ifndef _COMPOSITION_H
#define _COMPOSITION_H 1

#include <stddef.h>

unsigned int next_asc_composition(unsigned int *ar, int* kp) {
    // ar = [1, 1, 1,...]
    // kp = n - 1
    int i, j;
    int k = *kp;
    if (k == 0) return 0;

    ar[k - 1] += 1;
    if (ar[k] == 1) {
        ar[k] -= 1;
        k--;
    } else {
        j = k + ar[k] - 2;
        for (i=k; i<=j ; i++) ar[i] = 1;
        k = j;
    }
    *kp = k;
    return 1;
}

unsigned int next_desc_composition(unsigned int *ar, int* kp) {
    // ar = [n, 0, 0, ...]
    // kp = 1
    int i, j;
    int k = *kp;

    for (j = k-1; j >= 0; j--) {
        if (ar[j] > 1) {
            break;
        } else if (j == 0) {
            return 0;
        }
    }
    ar[j] = ar[j] - 1;
    ar[j + 1] = k - j;
    for (i = j + 2; i < k; i++) ar[i] = 0;
    *kp = j + 2;
    return 1;
}

#endif
