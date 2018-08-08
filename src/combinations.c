#define R_NO_REMAP
#define R_NO_REMAP_RMATH
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/combination.h"
#include "utils.h"
#include "macros.h"


SEXP next_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout) {
    int i, j;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    int n = as_uint(_n);
    int k = as_uint(_k);
    char layout = layout_flag(_layout);

    double dd = Rf_asInteger(_d) == -1 ? choose(n, k) : as_uint(_d);
    int d = verify_dimension(dd, k, layout);

    unsigned int* ap;

    if (!variable_exist(state, "a", INTSXP, k, (void**) &ap)) {
        for(i=0; i<k; i++) ap[i] = i;
        status = 0;
    }

    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_combination(ap, n, k)) { \
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
    attach_factor_levels(result, labels);
    UNPROTECT(nprotect);
    return result;
}


void ith_combination(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    unsigned int i, j;
    unsigned int start = 0;
    unsigned int count, this_count;

    for (i = 0; i < k; i++) {
        count = 0;
        for (j = start; j < n; j++) {
            this_count = count + choose(n - j - 1, k - i - 1);
            if (this_count > index) {
                ar[i] = j;
                start = j + 1;
                index -= count;
                break;
            }
            count = this_count;
        }
    }
}

void ith_combination_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    unsigned int i, j;
    unsigned int start = 0;
    mpz_t count, this_count;
    mpz_init(count);
    mpz_init(this_count);

    for (i = 0; i < k; i++) {
        mpz_set_ui(count, 0);
        for (j = start; j < n; j++) {
            mpz_bin_uiui(this_count, n - j - 1, k - i - 1);
            mpz_add(this_count, this_count, count);
            if (mpz_cmp(this_count, index) > 0) {
                ar[i] = j;
                start = j + 1;
                mpz_sub(index, index, count);
                break;
            }
            mpz_set(count, this_count);
        }
    }

    mpz_clear(count);
    mpz_clear(this_count);
}

SEXP get_combinations(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample) {
    int i, j;
    int nprotect = 0;
    SEXP result = R_NilValue;

    int n = as_uint(_n);
    int k = as_uint(_k);
    char layout = layout_flag(_layout);

    double dd = _index == R_NilValue ? as_uint(_nsample) : Rf_length(_index);
    int d = verify_dimension(dd, k, layout);

    unsigned int* ap;
    ap = (unsigned int*) R_alloc(k, sizeof(int));

    double max = choose(n, k);
    int sampling = _index == R_NilValue;
    int bigz = TYPEOF(_index) == STRSXP || max > INT_MAX;

    if (bigz) {
        gmp_randstate_t randstate;
        mpz_t z;
        mpz_t maxz;
        mpz_init(z);

        if (sampling) {
            GetRNGstate();
            set_gmp_randstate(randstate);
            mpz_init(maxz);
            mpz_bin_uiui(maxz, n , k);
        } else {
            if (TYPEOF(_index) != STRSXP) {
                _index = Rf_coerceVector(_index, STRSXP);
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                mpz_urandomm(z, randstate, maxz); \
            } else { \
                mpz_set_str(z, CHAR(STRING_ELT(_index, j)), 10); \
                mpz_sub_ui(z, z, 1); \
            } \
            ith_combination_bigz(ap, n, k, z);

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

        mpz_clear(z);
        if (sampling){
            mpz_clear(maxz);
            gmp_randclear(randstate);
            PutRNGstate();
        }

    } else {
        int* index;
        if (sampling) {
            GetRNGstate();
        } else {
            index = INTEGER(_index);
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                ith_combination(ap, n, k, floor(max * unif_rand())); \
            } else { \
                ith_combination(ap, n, k, index[j] - 1); \
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

        if (sampling){
            PutRNGstate();
        }
    }

    attach_factor_levels(result, labels);
    UNPROTECT(nprotect);
    return result;
}
