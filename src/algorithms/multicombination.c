#include "multicombination.h"

unsigned int next_multicombination(unsigned int *ar, size_t n, unsigned int k)
{
    unsigned int changed = 0;
    int i;

    for (i = k - 1; i >= 0 && !changed; i--) {
        if (ar[i] < n - 1) {
            // increment this element
            ar[i]++;
            if (i < k - 1) {
                // make the elements after it the same
                unsigned int j;
                for (j = i + 1; j < k; j++) {
                    ar[j] = ar[j - 1];
                }
            }
            changed = 1;
        }
    }
    if (!changed) {
        // reset to first combination
        for (i = 0; i < k; i++) {
            ar[i] = 0;
        }
    }
    return changed;
}
