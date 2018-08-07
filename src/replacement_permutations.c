#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/cartesian_product.h"
#include "utils.h"
#include "macros.h"


SEXP next_replacement_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout) {
    int i, j;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    int n = as_uint(_n);
    int k = as_uint(_k);
    char layout = check_layout(_layout);

    double dd = Rf_asInteger(_d) == -1 ? pow(n, k) : as_uint(_d);
    int d = check_dimension(dd, k, layout);

    size_t *sizes;
    sizes = (size_t*) R_alloc(k, sizeof(*sizes));
    for(i=0; i<k; i++) sizes[i] = n;

    unsigned int* ap;

    if (!variable_exist(state, "a", INTSXP, k, (void**) &ap)) {
        for(i=0; i<k; i++) ap[i] = 0;
        status = 0;
    }

    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_cartesian_product(ap, k, sizes)) { \
            status = 0; \
            break; \
        }

    int labels_type = TYPEOF(labels);
    if (labels_type == NILSXP) {
        RESULT_NILSXP(k);
    } else if (labels_type == INTSXP) {
        RESULT_INTSXP(k);
    } else if (labels_type == REALSXP) {
        RESULT_REALSXP(k);
    } else if (labels_type == STRSXP) {
        RESULT_STRSXP(k);
    }

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    check_factor(result, labels);
    UNPROTECT(nprotect);
    return result;
}

void ith_replacement_permutation(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    unsigned int i, j;

    for (i = 0; i < k; i++) {
        j = pow(n, k - 1 - i);
        ar[i] = index / j;
        index = index % j;
    }
}

void ith_replacement_permutation_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    unsigned int i;

    mpz_t q;
    mpz_init(q);
    mpz_t p;
    mpz_init(p);

    for (i = 0; i < k; i++) {
        mpz_ui_pow_ui(p, n, k - 1 - i);
        mpz_tdiv_qr(q, index, index, p);
        ar[i] = mpz_get_ui(q);
    }

    mpz_clear(q);
    mpz_clear(p);
}

SEXP get_ith_replacement_permutation(SEXP _n, SEXP _k, SEXP _index) {
    unsigned int i;
    int n = as_uint(_n);
    int k = as_uint(_k);
    SEXP as = PROTECT(Rf_allocVector(INTSXP, k));
    unsigned int* ar = (unsigned int*) INTEGER(as);

    if (TYPEOF(_index) == STRSXP || pow(n, k) > INT_MAX) {
        mpz_t z;
        mpz_init(z);

        if (TYPEOF(_index) == STRSXP) {
            mpz_set_str(z, CHAR(STRING_ELT(_index, 0)), 10);
            mpz_sub_ui(z, z, 1);
        } else {
            mpz_set_ui(z, as_uint(_index) - 1);
        }
        ith_replacement_permutation_bigz(ar, n, k, z);
        mpz_clear(z);
    } else {
        ith_replacement_permutation(ar, n, k, as_uint(_index) - 1);
    }

    for (i = 0; i < k; i++) {
        ar[i]++;
    }
    UNPROTECT(1);
    return as;
}
