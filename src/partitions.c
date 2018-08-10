#include <gmp.h>
#include "gmp_utils.h"
#include "partitions/partitions.h"
#include "partitions/k_partitions.h"


SEXP npartitions(SEXP _n, SEXP _k, SEXP _bigz) {
    SEXP ans;

    int n = as_uint(_n);
    int k = Rf_isNull(_k) ? -1 : as_uint(_k);

    if (Rf_asLogical(_bigz)) {
        mpz_t z;
        mpz_init(z);
        if (k == -1) {
            n_partitions_bigz(z, n);
        } else {
            n_k_partitions_bigz(z, n, k);
        }
        ans = mpz_to_bigz1(z);
        mpz_clear(z);
    } else {
        double d;
        if (k == -1) {
            d = n_partitions(n);
        } else {
            d = n_k_partitions(n, k);
        }
        if (d > INT_MAX) {
            Rf_error("integer overflow: use bigz instead");
        }
        ans = Rf_ScalarInteger((int) d);
    }

    return ans;
}


SEXP get_partitions(SEXP _n, SEXP _k, SEXP _descending, SEXP _layout, SEXP _d,
                    SEXP _index, SEXP _nsample, SEXP state, SEXP _skip) {

    SEXP ans = R_NilValue;

    int n = as_uint(_n);
    int k = Rf_isNull(_k) ? -1 : as_uint(_k);
    int descending = Rf_asInteger(_descending);
    char layout = layout_flag(_layout);
    int d = Rf_asInteger(_d);

    if (k == -1) {
        if (n == 0) {
            if (layout == 'r') {
                ans = Rf_allocMatrix(INTSXP, 1, 0);
            } else if (layout == 'c') {
                ans = Rf_allocMatrix(INTSXP, 0, 1);
            } else if (layout == 'l') {
                ans = PROTECT(Rf_allocVector(VECSXP, 1));
                SEXP ansi = PROTECT(Rf_allocVector(INTSXP, 0));
                SET_VECTOR_ELT(ans, 0, ansi);
                UNPROTECT(2);
            }
        } else if (descending) {
            ans = next_desc_partitions(n, layout, d, state);
        } else {
            ans = next_asc_partitions(n, layout, d, state);
        }
    } else {
        if (n == 0 && k == 0) {
            if (layout == 'r') {
                ans = Rf_allocMatrix(INTSXP, 1, 0);
            } else if (layout == 'c') {
                ans = Rf_allocMatrix(INTSXP, 0, 1);
            } else if (layout == 'l') {
                ans = PROTECT(Rf_allocVector(VECSXP, 1));
                SEXP ansi = PROTECT(Rf_allocVector(INTSXP, 0));
                SET_VECTOR_ELT(ans, 0, ansi);
                UNPROTECT(2);
            }
        } else if (k > n || k == 0) {
            if (layout == 'r') {
                ans = Rf_allocMatrix(INTSXP, 0, k);
            } else if (layout == 'c') {
                ans = Rf_allocMatrix(INTSXP, k, 0);
            } else if (layout == 'l') {
                ans = Rf_allocVector(VECSXP, 0);
            }
        } else if (descending) {
            ans = next_desc_k_partitions(n, k, layout, d, state);
        } else {
            ans = next_asc_k_partitions(n, k, layout, d, state);
        }
    }
    return ans;
}
