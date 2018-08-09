#ifndef _M_COMB_H
#define _M_COMB_H 1

#include <stddef.h>

unsigned int next_multiset_combination(
    const unsigned int *m, unsigned int *ar, size_t n, unsigned int k)
{
    unsigned int finished = 0;
    unsigned int changed = 0;
    unsigned int i, j, l;

    for (i = k - 1; !finished && !changed; i--) {
        if (ar[i] < m[i + (n - k)]) {
            // find the successor
            for (j = 0; m[j] <= ar[i]; j++);
            // replace this element with it
            ar[i] = m[j];
            if (i < k - 1) {
                // make the elements after it the same as this part of the m
                for (l = i + 1, j = j + 1; l < k; l++, j++) {
                    ar[l] = m[j];
                }
            }
            changed = 1;
        }
        finished = i == 0;
    }
    if (!changed) {
        // reset to first combination
        for (i = 0; i < k; i++) {
            ar[i] = m[i];
        }
    }
    return changed;
}


#endif
