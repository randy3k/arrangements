#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "../combinatorics.h"
#include "../utils.h"
#include "../macros.h"


void identify_ordinary_permutation(unsigned int* ar, unsigned int n, unsigned int index) {
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


void identify_ordinary_permutation_bigz(unsigned int* ar, unsigned int n, mpz_t index) {
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


SEXP next_ordinary_permutations(int n, int k, SEXP labels, SEXP freq, char layout, int d, SEXP state) {
    int i, j, h;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    int* fp;
    int flen;
    if (freq != R_NilValue) {
        fp = as_uint_array(freq);
        flen = Rf_length(freq);
    }

    double dd;
    if (d == -1) {
        if (freq == R_NilValue) {
            dd = fact(n);
        } else {
            dd = multichoose(fp, flen);
        }
    } else {
        dd = d;
    }
    d = verify_dimension(dd, n, layout);

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

    #undef NEXT
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
    attach_factor_levels(result, labels);
    UNPROTECT(nprotect);
    return result;
}


SEXP obtain_ordinary_permutations(int n, int k, SEXP labels, char layout, SEXP _index, SEXP _nsample) {
    int i, j;
    int nprotect = 0;
    SEXP result = R_NilValue;

    double dd = _index == R_NilValue ? as_uint(_nsample) : Rf_length(_index);
    int d = verify_dimension(dd, n, layout);

    unsigned int* ap;
    ap = (unsigned int*) R_alloc(n, sizeof(int));

    double max = fact(n);
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
            mpz_fac_ui(maxz, n);
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
            identify_ordinary_permutation_bigz(ap, n, z);

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
                identify_ordinary_permutation(ap, n, floor(max * unif_rand())); \
            } else { \
                identify_ordinary_permutation(ap, n, index[j] - 1); \
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

        if (sampling){
            PutRNGstate();
        }
    }

    attach_factor_levels(result, labels);
    UNPROTECT(nprotect);
    return result;
}
