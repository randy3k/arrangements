#include "k_partition.h"

unsigned int next_asc_k_partition(unsigned int *ar, size_t n, unsigned int k) {
    // Hindenburg algorithm
    unsigned int a, t;
    unsigned int i, j;
    a = ar[k-1];

    for (j = k-1; j && a - ar[j-1] < 2; j--);

    if (j == 0) {
        for(j=0; j<k-1; j++) ar[j] = 1;
        ar[k-1] = n - k + 1;
        return 0;
    }
    j = j - 1;
    a = ar[j];
    for (i=j; i<k-1; i++) ar[i] = a + 1;

    t = 0;
    for (i=0; i<k-1; i++) t += ar[i];
    ar[k - 1] = n - t;
    return 1;
}


unsigned int next_desc_k_partition(unsigned int *ar, size_t n, unsigned int k) {
    // an adhoc algorithm
    unsigned int a, b, t;
    unsigned int i, j;
    a = ar[k-1];
    t = a;
    for (j = k-1; j > 0; j--) {
        b = ar[j-1];
        if (b - a >= 2) break;
        t += b;
    }
    if (j == 0) {
        for(j=0; j<k-1; j++) ar[j] = 1;
        ar[0] = n - k + 1;
        return 0;
    }
    j = j - 1;
    b -= 1;
    ar[j] = b;
    t = t + 1 - (k - j - 1);
    while (t > b - 1) {
        j += 1;
        ar[j] = b;
        t -= b - 1;
    }
    if (j + 1 < k) {
        ar[j + 1] = t + 1;
        for(i=j+2; i<k; i++) {
            if (ar[i] == 1) {
                break;
            } else {
                ar[i] = 1;
            }
        }
    }

    return 1;
}
