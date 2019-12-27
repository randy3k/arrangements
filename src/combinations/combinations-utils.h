#ifndef COMBINATIONS_UTILS_H
#define COMBINATIONS_UTILS_H

static double choose(int n, int k) {
    double out = 1;
    int i;
    if (n >= 0 && n < k) {
        return 0;
    } else if (k == 0) {
        return 1;
    }

    if (k > n / 2) k = n - k;

    for (i = 1; i <= k; i++) {
        out *= n - k + i;
        out /= i;
    }

    return out;
}

#endif /* end of include guard: COMBINATIONS_UTILS_H */
