#ifndef _CART_PROD_H_
#define _CART_PROD_H_ 1

#include <stddef.h>

unsigned int next_cartesian_product(unsigned int *ar, size_t len, const size_t *sizes)
{
    unsigned int changed = 0;
    unsigned int finished = 0;
    unsigned int i;

    for (i = len - 1; !changed && !finished; i--) {
        if (ar[i] < sizes[i] - 1) {
            ar[i]++;
            changed = 1;
        }
        else {
            ar[i] = 0;
        }
        finished = i == 0;
    }

    return changed;
}

#endif
