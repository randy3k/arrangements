#include "gmp_utils.h"
#include "utils.h"

static int raw_size(mpz_t z) {
    int numb = 8 * sizeof(int);
    return sizeof(int) * (2 + (mpz_sizeinbase(z, 2) + numb - 1) / numb);
}

// from biginteger::as_raw
// https://github.com/cran/gmp/blob/cb2935d6a7948e85e42b1402c6f5a3450547f0cf/src/biginteger.cc#L75
SEXP mpz_to_bigz1(mpz_t z) {

    int size = raw_size(z);
    SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size + sizeof(int)));
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

// from biginteger::biginteger
// https://github.com/cran/gmp/blob/cb2935d6a7948e85e42b1402c6f5a3450547f0cf/src/biginteger.cc#L25
int as_mpz_array(mpz_t* a, size_t n, SEXP x) {
    size_t i;
    if (TYPEOF(x) == RAWSXP && Rf_inherits(x, "bigz")) {
        int* v = ((int*) RAW(x)) + 1;
        for (i = 0; i < n; i++) {
            if (v[0] > 0) {
                mpz_import(a[i], v[0], 1, sizeof(int) , 0, 0, &(v[2]));
                if(v[1] == -1) {
                    mpz_neg(a[i], a[i]);
                }
                v = v + v[0] + 2;
            } else {
                mpz_set_ui(a[i], 0);
                v = v + 1;
            }
        }
        return 0;
    } else if (TYPEOF(x) == INTSXP) {
        int* xp = INTEGER(x);
        for (i = 0; i < n; i++) {
            mpz_set_ui(a[i], abs(xp[i]));
            if (xp[i] < 0) {
                mpz_neg(a[i], a[i]);
            }
        }
        return 0;
    } else if (TYPEOF(x) == REALSXP) {
        double* xp = REAL(x);
        int w;
        for (i = 0; i < n; i++) {
            w = (int) fabs(xp[i]);
            mpz_set_ui(a[i], w == xp[i] ? w : 0);
            if (xp[i] < 0) {
                mpz_neg(a[i], a[i]);
            }
        }
        return 0;
    } else if (TYPEOF(x) == STRSXP) {
        for (i = 0; i < n; i++) {
            if (mpz_set_str(a[i], CHAR(STRING_ELT(x, i)), 10) < 0) {
                return -1;
            }
        }
        return 0;
    }
    return -1;
}
