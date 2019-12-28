#include <gmp.h>
#include "macros.h"
#include "utils.h"
#include "combinations.h"
#include "combinations/combinations-ordinary.c"
#include "combinations/combinations-multiset.c"
#include "combinations/combinations-replacement.c"


SEXP ncombinations(SEXP _x, SEXP _k, SEXP _n, SEXP _v, SEXP _freq, SEXP _replace, SEXP _bigz) {

    int i;
    SEXP ans;

    int has_vector = !Rf_isNull(_v);
    int multiset = !Rf_isNull(_freq);

    int n, k;
    int* fp;
    int flen;
    int replace = Rf_asLogical(_replace);

    VALIDATE_ARGUMENTS();

    if (multiset && replace) {
        n = flen;
    }

    if (Rf_asLogical(_bigz)) {
        mpz_t z;
        mpz_init(z);
        if (replace) {
            mpz_bin_uiui(z, n + k - 1, k);
        } else if (n < k) {
            mpz_set_ui(z, 0);
        } else if (multiset) {
            n_multiset_combinations_bigz(z, fp, flen, k);
        } else {
            mpz_bin_uiui(z, n, k);
        }
        ans = mpz_to_bigz1(z);
        mpz_clear(z);
    } else {
        double d;
        if (replace) {
            d = choose(n + k - 1, k);
        } else if (n < k) {
            d = 0;
        } else if (multiset) {
            d = n_multiset_combinations(fp, flen, k);
        } else {
            d = choose(n ,k);
        }
        if (d > INT_MAX) {
            Rf_error("integer overflow: use bigz instead");
        }
        ans = Rf_ScalarInteger((int) d);
    }

    return ans;
}


SEXP get_combinations(SEXP _x, SEXP _k, SEXP _n, SEXP _v, SEXP _freq, SEXP _replace,
                      SEXP _layout, SEXP _d, SEXP _index, SEXP _nsample, SEXP state, SEXP _skip, SEXP _drop) {
    int i;
    SEXP ans = R_NilValue;

    int has_vector = !Rf_isNull(_v);
    int multiset = !Rf_isNull(_freq);

    int n, k;
    int* fp;
    int flen;
    int replace = Rf_asLogical(_replace);

    VALIDATE_ARGUMENTS();

    char layout = layout_flag(_layout);
    int d = Rf_asInteger(_d) == -1 ? -1 : as_uint(_d);

    if (Rf_isNull(_index) && Rf_isNull(_nsample)) {
        if (k == 0) {
            if (layout == 'r') {
                ans = Rf_allocMatrix(has_vector ? TYPEOF(_v) : INTSXP, 1, 0);
            } else if (layout == 'c') {
                ans = Rf_allocMatrix(has_vector ? TYPEOF(_v) : INTSXP, 0, 1);
            } else if (layout == 'l') {
                if (n == 0) {
                    ans = PROTECT(Rf_allocVector(VECSXP, 1));
                    SEXP ansi = PROTECT(Rf_allocVector(has_vector ? TYPEOF(_v) : INTSXP, 0));
                    SET_VECTOR_ELT(ans, 0, ansi);
                    UNPROTECT(2);
                } else {
                    ans = Rf_allocVector(VECSXP, 0);
                }
            }
        } else if (k > n && (!replace || n == 0)) {
            if (layout == 'r') {
                ans = Rf_allocMatrix(has_vector ? TYPEOF(_v) : INTSXP, 0, k);
            } else if (layout == 'c') {
                ans = Rf_allocMatrix(has_vector ? TYPEOF(_v) : INTSXP, k, 0);
            } else if (layout == 'l') {
                ans = Rf_allocVector(VECSXP, 0);
            }
        } else if (replace) {
            ans = next_replacement_combinations(n, k, _v, layout, d, _skip, state);
        } else if (multiset) {
            ans = next_multiset_combinations(fp, flen, k, _v, layout, d, _skip, state);
        } else {
            ans = next_ordinary_combinations(n, k, _v, layout, d, _skip, state);
        }
    } else {
        if (replace) {
            ans = draw_replacement_combinations(n, k, _v, layout, _index, _nsample);
        } else if (multiset) {
            ans = draw_multiset_combinations(fp, flen, k, _v, layout, _index, _nsample);
        } else {
            ans = draw_ordinary_combinations(n, k, _v, layout, _index, _nsample);
        }
    }

    PROTECT(ans);
    attach_factor_levels(ans, _v);
    if (d > 0 && !Rf_isNull(state)) {
        if ((layout == 'r' && (Rf_nrows(ans) == 0)) ||
                (layout == 'c' && Rf_ncols(ans) == 0) ||
                (layout == 'l' && Rf_length(ans) == 0)) {
            ans = R_NilValue;
        } else if (k == 0 ||
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
