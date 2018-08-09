#include "gmp_utils.h"

static int raw_size(mpz_t z) {
    int numb = 8 * sizeof(int);
    return sizeof(int) * (2 + (mpz_sizeinbase(z, 2) + numb - 1) / numb);
}

// mirror of biginteger::as_raw
// https://github.com/cran/gmp/blob/cb2935d6a7948e85e42b1402c6f5a3450547f0cf/src/biginteger.cc#L75
SEXP mpz_to_bigz1(mpz_t z) {

    int size = raw_size(z);
    SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size + 1));
    unsigned char* raw = RAW(ans);
    int* r = (int*) raw;
    r[0] = 1; // scalar RAWSXP
    r[1] = size / sizeof(int) - 2;
    r[2] = (int) mpz_sgn(z);
    mpz_export(&r[3], 0, 1, sizeof(int), 0, 0, z);

    Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));
    UNPROTECT(1);
    return ans;
}
