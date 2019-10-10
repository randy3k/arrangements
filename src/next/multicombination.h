#ifndef _MT_COMB_H
#define _MT_COMB_H 1

#include <stddef.h>

unsigned int next_multicombination(unsigned int *ar, size_t n, unsigned int k)
{
    unsigned int i;
    unsigned int j;
    unsigned int temp;

    for (i = k - 1; ; i--) {
        if (ar[i] < n - 1) {
            // increment this element
            ar[i]++;
            if (i < k - 1) {
                // make the elements after it the same
                temp = ar[i];
                for (j = i + 1; j < k; j++) {
                    ar[j] = temp;
                }
            }
            return 1;
        } else if (i == 0) {
            return 0;
        }
    }
}


#endif
