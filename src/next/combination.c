#include "combination.h"

unsigned int next_combination(unsigned int *ar, size_t n, unsigned int k)
{
    unsigned int finished = 0;
    unsigned int changed = 0;
    unsigned int i;

    if (k > 0) {
        for (i = k - 1; !finished && !changed; i--) {
            if (ar[i] < (n - 1) - (k - 1) + i) {
                // increment this element
                ar[i]++;
                if (i < k - 1) {
                    // turn the elements after it into a linear sequence
                    unsigned int j;
                    for (j = i + 1; j < k; j++) {
                        ar[j] = ar[j - 1] + 1;
                    }
                }
                changed = 1;
            }
            finished = i == 0;
        }
        if (!changed) {
            // reset to first combination
            for (i = 0; i < k; i++) {
                ar[i] = i;
            }
        }
    }
    return changed;
}
