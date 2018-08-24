#ifndef _PART_H
#define _PART_H 1

#include <stddef.h>

unsigned int next_asc_partition(unsigned int *ar, int* kp) {
    // by J Kellehers 2005 Encoding Partitions As Ascending Compositions
    // ar = [1, 1, 1,....], *k = n - 1
    // or ar = [m-1, n-m+1, 0,....], *k = 1 where m is the initial part

    unsigned int x, y;
    int k = *kp;
    if (k == 0) {
        x = ar[0];
        for (y = 0; y < x; y++) ar[y] = 1;
        *kp = x - 1;
        return 0;
    }
    y = ar[k] - 1;
    k--;
    x = ar[k] + 1;
    while (x <= y) {
        ar[k] = x;
        y -= x;
        k++;
    }
    ar[k] = x + y;
    *kp = k;
    return 1;
}

unsigned int next_desc_partition(unsigned int *ar, int* hp, int* kp) {
    // by A Zoghbi and I Stojmenovic 1994 Fast Algorithms for Generating Integer Partitions
    // ar = [n, 1, 1,....], *h = 0, *k = 1

    unsigned int x, r, t;
    int h = *hp;
    int k = *kp;

    if (ar[0] == 1) {
        for (x = 0; x < k; x++) ar[x] = 1;
        ar[0] = k;
        return 0;
    }

    if (ar[h] == 2) {
        k += 1;
        ar[h] = 1;
        h -= 1;
    } else {
        r = ar[h] - 1;
        t = k - h;
        ar[h] = r;
        while (t >= r) {
            h += 1;
            ar[h] = r;
            t -= r;
        }
        if (t == 0) {
            k = h + 1;
        } else {
            k = h + 2;
            if (t > 1) {
                h += 1;
                ar[h] = t;
            }
        }
    }
    *hp = h;
    *kp = k;
    return 1;
}

#endif
