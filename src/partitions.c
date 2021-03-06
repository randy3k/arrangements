#include <gmp.h>
#include "math.h"
#include "utils.h"
#include "partitions.h"
#include "partitions/partitions-asc.c"
#include "partitions/partitions-desc.c"
#include "partitions/partitions-asc-k.c"
#include "partitions/partitions-desc-k.c"
#include "partitions/partitions-asc-distinct.c"
#include "partitions/partitions-desc-distinct.c"
#include "partitions/partitions-asc-k-distinct.c"
#include "partitions/partitions-desc-k-distinct.c"
#include "partitions/partitions-utils.c"


SEXP npartitions(SEXP _n, SEXP _k, SEXP _distinct, SEXP _bigz) {
    SEXP ans;

    int n = as_uint(_n);
    int k = Rf_isNull(_k) ? -1 : as_uint(_k);

    if (Rf_asLogical(_bigz)) {
        mpz_t z;
        mpz_init(z);
        if (Rf_asLogical(_distinct)) {
            if (k == -1) {
                n_distinct_partitions_bigz(z, n);
            } else {
                n_k_distinct_partitions_bigz(z, n, k);
            }
        } else {
            if (k == -1) {
                n_partitions_bigz(z, n);
            } else {
                n_k_partitions_bigz(z, n, k);
            }
        }

        ans = mpz_to_bigz1(z);
        mpz_clear(z);
    } else {
        double d;
        if (Rf_asLogical(_distinct)) {
            if (k == -1) {
                d = n_distinct_partitions(n);
            } else {
                d = n_k_distinct_partitions(n, k);
            }
        } else {
            if (k == -1) {
                d = n_partitions(n);
            } else {
                d = n_k_partitions(n, k);
            }
        }
        if (d > INT_MAX) {
            Rf_error("integer overflow: use bigz instead");
        }
        ans = Rf_ScalarInteger((int) d);
    }

    return ans;
}


SEXP get_partitions(SEXP _n, SEXP _k, SEXP _distinct, SEXP _descending, SEXP _layout, SEXP _d,
                    SEXP _index, SEXP _nsample, SEXP state, SEXP _skip, SEXP _drop) {

    SEXP ans = R_NilValue;

    int n = as_uint(_n);
    int k = Rf_isNull(_k) ? -1 : as_uint(_k);
    int distinct = Rf_asInteger(_distinct);
    int descending = Rf_asInteger(_descending);
    char layout = layout_flag(_layout);
    int d = Rf_asInteger(_d);

    if (Rf_isNull(_index) && Rf_isNull(_nsample)) {
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
            } else if (distinct) {
                if (descending) {
                    ans = next_desc_distinct_partitions(n, layout, d, _skip, state);
                } else {
                    ans = next_asc_distinct_partitions(n, layout, d, _skip, state);
                }
            } else {
                if (descending) {
                    ans = next_desc_partitions(n, layout, d, _skip, state);
                } else {
                    ans = next_asc_partitions(n, layout, d, _skip, state);
                }
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
            } else if (k > n || k == 0 || (distinct && n - choose(k, 2) < k)) {
                if (layout == 'r') {
                    ans = Rf_allocMatrix(INTSXP, 0, k);
                } else if (layout == 'c') {
                    ans = Rf_allocMatrix(INTSXP, k, 0);
                } else if (layout == 'l') {
                    ans = Rf_allocVector(VECSXP, 0);
                }
            } else if (distinct) {
                if (descending) {
                    ans = next_desc_k_distinct_partitions(n, k, layout, d, _skip, state);
                } else {
                    ans = next_asc_k_distinct_partitions(n, k, layout, d, _skip, state);
                }
            } else {
                if (descending) {
                    ans = next_desc_k_partitions(n, k, layout, d, _skip, state);
                } else {
                    ans = next_asc_k_partitions(n, k, layout, d, _skip, state);
                }
            }
        }
    } else {
        if (distinct) {
            if (k == -1) {
                if (descending) {
                    ans = draw_desc_distinct_partitions(n, layout, _index, _nsample);
                } else {
                    ans = draw_asc_distinct_partitions(n, layout, _index, _nsample);
                }
            } else {
                if (descending) {
                    ans = draw_desc_k_distinct_partitions(n, k, layout, _index, _nsample);
                } else {
                    ans = draw_asc_k_distinct_partitions(n, k, layout, _index, _nsample);
                }
            }
        } else {
            if (k == -1) {
                if (descending) {
                    ans = draw_desc_partitions(n, layout, _index, _nsample);
                } else {
                    ans = draw_asc_partitions(n, layout, _index, _nsample);
                }
            } else {
                if (descending) {
                    ans = draw_desc_k_partitions(n, k, layout, _index, _nsample);
                } else {
                    ans = draw_asc_k_partitions(n, k, layout, _index, _nsample);
                }
            }
        }
    }

    PROTECT(ans);
    if (d > 0 && !Rf_isNull(state)) {
        if ((layout == 'r' && (Rf_nrows(ans) == 0)) ||
                (layout == 'c' && Rf_ncols(ans) == 0) ||
                (layout == 'l' && Rf_length(ans) == 0)) {
            ans = R_NilValue;
        } else if ((n == 0 && k <= 0) ||
                (layout == 'r' && (Rf_nrows(ans) < d)) ||
                (layout == 'c' && Rf_ncols(ans) < d) ||
                (layout == 'l' && Rf_length(ans) < d)) {
            Rf_defineVar(Rf_install("null_pending"), Rf_ScalarLogical(1), state);
        }
    }
    if ((!Rf_isNull(_drop) && Rf_asLogical(_drop)) ||
            ((Rf_isNull(_drop) || Rf_asLogical(_drop)) && ((d == 1 && Rf_isNull(_layout)) ||
            (!Rf_isNull(_index) && index_length(_index) == 1 && Rf_isNull(_layout)) ||
            (!Rf_isNull(_nsample) && as_uint(_nsample) == 1 && Rf_isNull(_layout))) )) {
        if (layout == 'r' && Rf_nrows(ans) == 1) {
            Rf_setAttrib(ans, R_DimSymbol, R_NilValue);
        } else if (layout == 'c' && Rf_ncols(ans) == 1) {
            Rf_setAttrib(ans, R_DimSymbol, R_NilValue);
        } else if (layout == 'l' && Rf_length(ans) == 1) {
            ans = VECTOR_ELT(ans, 0);
        }
    }
    UNPROTECT(1);
    return ans;
}
