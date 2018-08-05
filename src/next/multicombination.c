#include "multicombination.h"

unsigned int next_multicombination(unsigned int *ar, size_t n, unsigned int k)
{
    unsigned int finished = 0;
    unsigned int changed = 0;
    unsigned int i;
    unsigned int j;
    unsigned int temp;

    for (i = k - 1; !finished && !changed; i--) {
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
            changed = 1;
        }
        finished = i == 0;
    }
    if (!changed) {
        // reset to first combination
        for (i = 0; i < k; i++) {
            ar[i] = 0;
        }
    }
    return changed;
}
