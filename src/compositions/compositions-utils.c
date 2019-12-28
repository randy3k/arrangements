#include "compositions-utils.h"
#include "../math.h"

double n_compositions(int n) {
    return n ? pow(2, n - 1) : 1;
}


void n_compositions_bigz(mpz_t z, int n) {
    if (n > 0) {
        mpz_ui_pow_ui(z, 2, n - 1);
    } else {
        mpz_set_ui(z, 1);
    }
}


double n_k_compositions(int n, int k) {
    if (n < k) {
        return 0;
    } else if (k == 0) {
        return n == 0;
    } else {
        return choose(n - 1, k - 1);
    }
}


void n_k_compositions_bigz(mpz_t z, int n, int k) {
    if (n < k) {
        mpz_set_ui(z, 0);
    } else if (k == 0) {
        mpz_set_ui(z, n == 0);
    } else {
        mpz_bin_uiui(z, n - 1, k - 1);
    }
}
