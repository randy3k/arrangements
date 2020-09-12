#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "stdlib.h"
#include "../macros.h"
#include "../math.h"
#include "../utils.h"


unsigned int next_permutation(unsigned int *ar, size_t n) {
    unsigned int k, j;
    unsigned int result = 0;

    // trival for only one element
    if (n == 1) {
        return result;
    }
    // find the largest k such that a[k] < a[k + 1]
    for (k = n - 1; k && ar[k - 1] >= ar[k]; k--);

    if (!k--) {
        // if not found, array is highest permutation
        reverse(ar, n);
    } else {
        // fnd the largest index j such that a[k] < a[l]
        for (j = n - 1; ar[j] <= ar[k]; j--);
        // swap it with the first one
        swap(ar, k, j);
        // reverse the remainder
        reverse(ar + k + 1, n - k - 1);
        result = 1;
    }
    return result;
}

void nth_ordinary_permutation(unsigned int* ar, unsigned int n, unsigned int index) {
    int i, j;
    if (n == 0) return;

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


void nth_ordinary_permutation_bigz(unsigned int* ar, unsigned int n, mpz_t index) {
    unsigned int i, j;
    if (n == 0) return;

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


void nth_multiset_permutation(unsigned int* ar, int* freq, size_t flen, size_t k, unsigned int index);

void nth_multiset_permutation_bigz(unsigned int* ar, int* freq, size_t flen, size_t k, mpz_t index);

void n_multiset_n_permutations_bigz(mpz_t z, int* freq, size_t flen);


SEXP next_ordinary_permutations(int n, int k, SEXP labels, SEXP freq, char layout, int d, SEXP _skip, SEXP state) {
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
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        if (freq == R_NilValue) {
            maxd = fact(n);
        } else {
            maxd = multichoose(fp, flen);
        }
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;

    if (!variable_exists(state, (char*)"a", INTSXP, n, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
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
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);

                if (freq == R_NilValue) {
                    mpz_fac_ui(maxz, n);
                } else {
                    n_multiset_n_permutations_bigz(maxz, fp, flen);
                }

                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set_ui(skipz, 0);
                }
                mpz_clear(maxz);
                if (freq == R_NilValue) {
                    nth_ordinary_permutation_bigz(ap, n, skipz);
                } else {
                    nth_multiset_permutation_bigz(ap, fp, flen, n, skipz);
                }
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                if (freq == R_NilValue) {
                    nth_ordinary_permutation(ap, n, skip);
                } else {
                    nth_multiset_permutation(ap, fp, flen, n, skip);
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
    } else if (labels_type == LGLSXP) {
        RESULT_LGLSXP(n);
    } else {
        Rf_error("label type not supported");
    }

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    UNPROTECT(nprotect);
    return result;
}


SEXP draw_ordinary_permutations(int n, SEXP labels, char layout, SEXP _index, SEXP _nsample) {
    int i, j;
    int nprotect = 0;
    int bigz = 0;
    int sampling = _index == R_NilValue;
    SEXP result = R_NilValue;

    double dd;
    if (sampling) {
        dd = as_uint(_nsample);
    } else if (TYPEOF(_index) == RAWSXP && Rf_inherits(_index, "bigz")) {
        dd = *((int* ) RAW(_index));
        bigz = 1;
    } else {
        dd = Rf_length(_index);
    }
    int d = verify_dimension(dd, n, layout);

    double maxd;
    if (!bigz) {
        maxd = fact(n);
        bigz = maxd > INT_MAX;
    }

    unsigned int* ap;
    ap = (unsigned int*) R_alloc(n, sizeof(int));

    if (bigz) {
        mpz_t* index;
        gmp_randstate_t randstate;
        mpz_t z;
        mpz_t maxz;
        mpz_init(z);
        mpz_init(maxz);
        mpz_fac_ui(maxz, n);

        if (sampling) {
            GetRNGstate();
            set_gmp_randstate(randstate);
        } else {
            index = (mpz_t*) R_alloc(d, sizeof(mpz_t));
            for (i = 0; i < d; i++) mpz_init(index[i]);
            int status = as_mpz_array(index, d, _index);
            for(i = 0; i < d; i++) {
                if (status < 0 || mpz_sgn(index[i]) <= 0 || mpz_cmp(index[i], maxz) > 0) {
                    for (i = 0; i < d; i++) mpz_clear(index[i]);
                    mpz_clear(maxz);
                    mpz_clear(z);
                    Rf_error("invalid index");
                }
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                mpz_urandomm(z, randstate, maxz); \
            } else { \
                mpz_sub_ui(z, index[j], 1); \
            } \
            nth_ordinary_permutation_bigz(ap, n, z);

        int labels_type = TYPEOF(labels);
        if (labels_type == NILSXP) {
            RESULT_NILSXP(n);
        } else if (labels_type == INTSXP) {
            RESULT_INTSXP(n);
        } else if (labels_type == REALSXP) {
            RESULT_REALSXP(n);
        } else if (labels_type == STRSXP) {
            RESULT_STRSXP(n);
        } else if (labels_type == LGLSXP) {
            RESULT_LGLSXP(n);
        } else {
            Rf_error("label type not supported");
        }

        mpz_clear(z);
        mpz_clear(maxz);
        if (sampling){
            gmp_randclear(randstate);
            PutRNGstate();
        } else {
            for (i = 0; i < d; i++) mpz_clear(index[i]);
        }

    } else {
        int* index;
        if (sampling) {
            GetRNGstate();
        } else {
            index = as_uint_index(_index);
            for (i = 0; i < d; i++) {
                if (index[i] <= 0 || index[i] > maxd) {
                    Rf_error("invalid index");
                }
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                nth_ordinary_permutation(ap, n, floor(maxd * unif_rand())); \
            } else { \
                nth_ordinary_permutation(ap, n, index[j] - 1); \
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
        } else if (labels_type == LGLSXP) {
            RESULT_LGLSXP(n);
        } else {
            Rf_error("label type not supported");
        }

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}
