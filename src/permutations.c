#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/permutation.h"
#include "utils.h"
#include "macros.h"

SEXP next_permutations(SEXP _n, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout) {
    int i, j, h;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    int n = as_uint(_n);
    int* fp;
    if (freq != R_NilValue) {
        fp = INTEGER(freq);
    }
    char layout = check_layout(_layout);

    double dd;
    if (Rf_asInteger(_d) == -1) {
        if (freq == R_NilValue) {
            dd = fact(n);
        } else {
            dd = multichoose(fp, Rf_length(freq));
        }
    } else {
        dd = as_uint(_d);
    }
    int d = check_dimension(dd, n, layout);

    unsigned int* ap;

    if (!variable_exist(state, "a", INTSXP, n, (void**) &ap)) {
        if (freq == R_NilValue) {
            for(i=0; i<n; i++) ap[i] = i;
        } else {
            h = 0;
            for (i = 0; i< Rf_length(freq); i++) {
                for (j = 0; j< fp[i]; j++) {
                    ap[h++] = i;
                }
            }
        }
        status = 0;
    }

    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_permutation(ap, n)) { \
            status = 0; \
            break; \
        }

    int labels_type = TYPEOF(labels);
    if (labels_type == NILSXP) {
        RESULT_NILSXP(n);
    } else if (labels_type == INTSXP) {
        RESULT_INTSXP(n);
    } else if (labels_type == REALSXP) {
        RESULT_REALSXP(n);
    } else if (labels_type == STRSXP) {
        RESULT_STRSXP(n);
    }

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    check_factor(result, labels);
    UNPROTECT(nprotect);
    return result;
}

SEXP num_multiset_n_permutations(SEXP freq) {
    int* fp = INTEGER(as_uint_array(freq));
    size_t flen = Rf_length(freq);
    return Rf_ScalarReal(multichoose(fp, flen));
}

void n_multiset_n_permutations_bigz(mpz_t z, int* freq, size_t flen) {
    mpz_set_ui(z, 1);
    size_t i, j, h;
    h = 0;
    for (i=0; i<flen; i++) {
        for (j=1; j<=freq[i]; j++) {
            h++;
            mpz_mul_ui(z, z, h);
            mpz_cdiv_q_ui(z, z, j);
        }
    }
}

SEXP num_multiset_n_permutations_bigz(SEXP freq) {
    int* fp = INTEGER(as_uint_array(freq));
    size_t flen = Rf_length(freq);
    mpz_t z;
    mpz_init(z);
    n_multiset_n_permutations_bigz(z, fp, flen);
    char* c = mpz_get_str(NULL, 10, z);
    SEXP out = Rf_mkString(c);
    mpz_clear(z);
    free(c);
    return out;
}

void ith_permutation(unsigned int* ar, unsigned int n, unsigned int index) {
    unsigned int i, j;
    unsigned int* fact = (unsigned int*) malloc(n * sizeof(unsigned int));

    fact[0] = 1;
    for (i = 1; i < n; i++) {
        fact[i] = fact[i-1] * i;
    }
    for (i = 0; i < n; i++) {
        ar[i] = index / fact[n - 1 - i];
        index = index % fact[n - 1 - i];
    }

    for (i = n - 1; i > 0; i--) {
        j = i;
        while (j-- > 0) {
            if (ar[j] <= ar[i]) {
                ar[i]++;
            }
        }
    }

    free(fact);
}

void ith_permutation_bigz(unsigned int* ar, unsigned int n, mpz_t index) {
    unsigned int i, j;

    mpz_t q;
    mpz_init(q);

    mpz_t* fact = (mpz_t*) malloc(n * sizeof(mpz_t));
    for (i=0; i< n; i++) mpz_init(fact[i]);

    mpz_set_ui(fact[0], 1);
    for (i=1; i< n; i++) mpz_mul_ui(fact[i], fact[i-1], i);

    for (i = 0; i < n; i++) {
        mpz_tdiv_qr(q, index, index, fact[n - 1 - i]);
        ar[i] = mpz_get_ui(q);
    }

    for (i = n - 1; i > 0; i--) {
        j = i;
        while (j-- > 0) {
            if (ar[j] <= ar[i]) {
                ar[i]++;
            }
        }
    }

    mpz_clear(q);
    for (i=0; i< n; i++) mpz_clear(fact[i]);
    free(fact);
}

SEXP get_ith_permutation(SEXP _n, SEXP _index) {
    unsigned int i;
    int n = as_uint(_n);
    SEXP as = PROTECT(Rf_allocVector(INTSXP, n));
    unsigned int* ar = (unsigned int*) INTEGER(as);

    if (TYPEOF(_index) == STRSXP || n > 12) {
        mpz_t z;
        mpz_init(z);

        if (TYPEOF(_index) == STRSXP) {
            mpz_set_str(z, CHAR(STRING_ELT(_index, 0)), 10);
            mpz_sub_ui(z, z, 1);
        } else {
            mpz_set_ui(z, as_uint(_index) - 1);
        }
        ith_permutation_bigz(ar, n, z);
        mpz_clear(z);
    } else {
        ith_permutation(ar, n, as_uint(_index) - 1);
    }

    for (i = 0; i < n; i++) {
        ar[i]++;
    }
    UNPROTECT(1);
    return as;
}
