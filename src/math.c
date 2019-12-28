#include "math.h"

double choose(int n, int k) {
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


double fact(int n) {
    double out;
    int i;
    out = 1;
    for(i=0; i<n; i++) {
        out = out * (n - i);
    }
    return out;
}

double fallfact(int n, int k) {
    double out;
    int i;
    if (n < k) {
        return 0;
    }
    out = 1;
    for(i=0; i<k; i++) {
        out = out * (n - i);
    }
    return out;
}


double multichoose(int* freq, size_t flen) {
    double out = 1;
    size_t h, i, j;
    h = 0;
    for (i=0; i<flen; i++) {
        for (j=1; j<=freq[i]; j++) {
            h++;
            out = out * h / j;
        }
    }
    return out;
}
