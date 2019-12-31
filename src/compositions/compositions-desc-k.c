#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "stdlib.h"
#include "../utils.h"
#include "../macros.h"
#include "compositions-utils.h"


unsigned int next_desc_k_composition(unsigned int *ar, size_t n, unsigned int k, int* tp) {
    int i, j;
    int t = *tp;

    for (i = k-1; i>= 1; i--) {
        if (ar[i - 1] > 1) {
            break;
        }
    }
    if (i == 0) {
        return 0;
    }
    if (t > 0) {
        t--;
    } else {
        t = 0;
        for (j = i; j < k - 1; j++) {
            t += ar[j];
        }
    }
    ar[i - 1] = ar[i - 1] - 1;
    ar[i] = t + ar[k - 1] - (k - i - 1) + 1;
    for (j = i+1; j < k; j++) ar[j] = 1;

    *tp = t;
    return 1;
}


void nth_desc_k_composition(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    unsigned int i, j;
    unsigned int count, this_count;

    for (i = 0; i < k; i++) {
        count = 0;
        for (j = n; j >= 1; j--) {
            this_count = count + n_k_compositions(n - j, k - i - 1);
            if (this_count > index) {
                ar[i] = j;
                n -= j;
                index -= count;
                break;
            }
            count = this_count;
        }
    }
}


void nth_desc_k_composition_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    unsigned int i, j;
    mpz_t count, this_count;
    mpz_init(count);
    mpz_init(this_count);

    for (i = 0; i < k; i++) {
        mpz_set_ui(count, 0);
        for (j = n; j >= 1; j--) {
            n_k_compositions_bigz(this_count, n - j, k - i - 1);
            mpz_add(this_count, this_count, count);
            if (mpz_cmp(this_count, index) > 0) {
                ar[i] = j;
                n -= j;
                mpz_sub(index, index, count);
                break;
            }
            mpz_set(count, this_count);
        }
    }

    mpz_clear(count);
    mpz_clear(this_count);
}


SEXP next_desc_k_compositions(int n, int k, char layout, int d, SEXP _skip, SEXP state) {
    int i, j;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_k_compositions(n, k);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;
    int* tp;
    int t;

    if (!variable_exists(state, (char*)"a", INTSXP, k, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            for(i=1; i<k; i++) ap[i] = 1;
            ap[0] = n - k + 1;
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_k_compositions_bigz(maxz, n, k);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set_ui(skipz, 0);
                }
                mpz_clear(maxz);
                nth_desc_k_composition_bigz(ap, n, k, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_desc_k_composition(ap, n, k, skip);
            }
        }
        status = 0;
    }
    if (!variable_exists(state, (char*)"t", INTSXP, 1, (void**) &tp)) {
        if (Rf_isNull(_skip)) {
            tp[0] = k - 1;
        } else  {
            t = 0;
            for (i = k-1; i>= 1; i--) {
                if (ap[i - 1] > 1) {
                    break;
                }
                t += ap[i - 1];
            }
            tp[0] = t + 1;
        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_desc_k_composition(ap, n, k, tp)) { \
            status = 0; \
            break; \
        }

    RESULT_K_PART();

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    UNPROTECT(nprotect);
    return result;
}


SEXP draw_desc_k_compositions(int n, int k, char layout, SEXP _index, SEXP _nsample) {
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
    int d = verify_dimension(dd, k, layout);

    double maxd;
    if (!bigz) {
        maxd = n_k_compositions(n, k);
        bigz = maxd > INT_MAX;
    }

    unsigned int* ap;
    ap = (unsigned int*) R_alloc(k, sizeof(int));

    if (bigz) {
        mpz_t* index;
        gmp_randstate_t randstate;
        mpz_t z;
        mpz_t maxz;
        mpz_init(z);
        mpz_init(maxz);
        n_k_compositions_bigz(maxz, n, k);

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
            nth_desc_k_composition_bigz(ap, n, k, z);

        RESULT_K_PART();

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
                nth_desc_k_composition(ap, n, k, floor(maxd * unif_rand())); \
            } else { \
                nth_desc_k_composition(ap, n, k, index[j] - 1); \
            }

        RESULT_K_PART();

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}
