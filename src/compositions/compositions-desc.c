#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include <math.h>
#include "stdlib.h"
#include "../utils.h"
#include "../macros.h"
#include "compositions-utils.h"


unsigned int next_desc_composition(unsigned int *ar, int* kp) {
    // ar = [n, 0, 0, ...]
    // kp = 1
    int i, j;
    int k = *kp;

    for (j = k-1; j >= 0; j--) {
        if (ar[j] > 1) {
            break;
        } else if (j == 0) {
            return 0;
        }
    }
    ar[j] = ar[j] - 1;
    ar[j + 1] = k - j;
    for (i = j + 2; i < k; i++) ar[i] = 0;
    *kp = j + 2;
    return 1;
}

void nth_desc_composition(unsigned int* ar, unsigned int n, unsigned int index) {
    int i, j;
    unsigned int s;
    unsigned int n1 = n - 1;
    int* bs;

    if (n == 0) return;

    // convert index to binary
    bs = (int*) malloc(n1 * sizeof(int));
    for (j = 0; j < n1; j++) {
        bs[j] = (index >> j) & 1;
    }
    s = 0;
    i = 0;
    // compute successive diff
    for (j = 0; j < n1; j++) {
        if (bs[n1 - j - 1] == 0) continue;
        ar[i] = j + 1 - s;
        s = j + 1;
        i++;
    }
    ar[i] = n1 - s + 1;
    for (j = i + 1; j < n; j++) ar[j] = 0;
    free(bs);
}


void nth_desc_composition_bigz(unsigned int* ar, unsigned int n, mpz_t index) {
    int i, j;
    unsigned int s;
    int n1 = n - 1;
    int* bs;

    if (n == 0) return;

    // convert index to binary
    bs = (int*) malloc(n1 * sizeof(int));
    for (j = 0; j < n1; j++) {
        bs[j] = mpz_tstbit(index, j);
    }
    s = 0;
    i = 0;
    // compute successive diff
    for (j = 0; j < n1; j++) {
        if (bs[n1 - j - 1] == 0) continue;
        ar[i] = j + 1 - s;
        s = j + 1;
        i++;
    }
    ar[i] = n1 - s + 1;
    for (j = i + 1; j < n; j++) ar[j] = 0;
    free(bs);
}

SEXP next_desc_compositions(int n, char layout, int d, SEXP _skip, SEXP state) {
    int i, j, k;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_compositions(n);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;
    int* kp;

    if (!variable_exists(state, (char*)"a", INTSXP, n, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            ap[0] = n;
            for(i=1; i<n; i++) ap[i] = 1;
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_compositions_bigz(maxz, n);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set_ui(skipz, 0);
                }
                mpz_clear(maxz);
                nth_desc_composition_bigz(ap, n, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_desc_composition(ap, n, skip);
            }
        }
        status = 0;
    }

    if (!variable_exists(state, (char*)"k", INTSXP, 1, (void**) &kp)) {
        if (Rf_isNull(_skip)) {
            kp[0] = 1;
        } else  {
            for (i = 0; i < n; i++) {
                if (ap[i] == 0) {
                    break;
                }
            }
            kp[0] = i;
        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_desc_composition(ap, kp)) { \
            status = 0; \
            break; \
        } \
        k = kp[0];

    RESULT_PART(n, k);

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    UNPROTECT(nprotect);
    return result;
}


SEXP draw_desc_compositions(int n, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = n_compositions(n);
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
        n_compositions_bigz(maxz, n);

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

        int k;

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                mpz_urandomm(z, randstate, maxz); \
            } else { \
                mpz_sub_ui(z, index[j], 1); \
            } \
            nth_desc_composition_bigz(ap, n, z); \
            for (i = 0; i < n; i++) { \
                if (ap[i] == 0) { \
                    break; \
                } \
            } \
            k = i;

        RESULT_PART(n, k);

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

        int k;

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                nth_desc_composition(ap, n, floor(maxd * unif_rand())); \
            } else { \
                nth_desc_composition(ap, n, index[j] - 1); \
            } \
            for (i = 0; i < n; i++) { \
                if (ap[i] == 0) { \
                    break; \
                } \
            } \
            k = i;

        RESULT_PART(n, k);

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}
