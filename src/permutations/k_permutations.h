#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "../combinatorics.h"
#include "../utils.h"
#include "../macros.h"


void n_k_permutations_bigz(mpz_t p, size_t n, size_t k) {
    size_t i;
    if (n < k) {
        mpz_set_ui(p, 0);
        return;
    }
    mpz_set_ui(p, 1);
    for(i=0; i<k; i++) {
        mpz_mul_ui(p, p, n - i);
    }
}


void identify_k_permutation(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    unsigned int i, j;

    for (i = 0; i < k; i++) {
        j = fallfact(n - 1 - i, k - 1 - i);
        ar[i] = index / j;
        index = index % j;
    }

    for (i = k - 1; i > 0; i--) {
        j = i;
        while (j-- > 0) {
            if (ar[j] <= ar[i]) {
                ar[i]++;
            }
        }
    }
}

void identify_k_permutation_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    unsigned int i, j;

    mpz_t q;
    mpz_init(q);
    mpz_t p;
    mpz_init(p);

    for (i = 0; i < k; i++) {
        n_k_permutations_bigz(p, n - 1 - i, k - 1 - i);
        mpz_tdiv_qr(q, index, index, p);
        ar[i] = mpz_get_ui(q);
    }

    for (i = k - 1; i > 0; i--) {
        j = i;
        while (j-- > 0) {
            if (ar[j] <= ar[i]) {
                ar[i]++;
            }
        }
    }

    mpz_clear(q);
    mpz_clear(p);
}


SEXP next_k_permutations(int n, int k, SEXP labels, char layout, int d, SEXP state) {
    int i, j;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd = d == -1 ? fallfact(n, k) : d;
    d = verify_dimension(dd, k, layout);

    unsigned int* ap;
    unsigned int* cyclep;

    if (!variable_exist(state, "a", INTSXP, n, (void**) &ap)) {
        for(i=0; i<n; i++) ap[i] = i;
        status = 0;
    }
    if (!variable_exist(state, "cycle", INTSXP, k, (void**) &cyclep)) {
        for(i=0; i<k; i++) cyclep[i] = n - i;;
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_k_permutation(ap, cyclep, n, k)) { \
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


SEXP obtain_k_permutations(int n, int k, SEXP labels, char layout, SEXP _index, SEXP _nsample) {
    int i, j;
    int nprotect = 0;
    int bigz = 0;
    int sampling = _index == R_NilValue;
    SEXP result = R_NilValue;

    double dd;
    if (sampling) {
        dd = as_uint(_nsample);
    } else if (TYPEOF(_index) == RAWSXP || Rf_inherits(_index, "bigz")) {
        dd = *((int* ) RAW(_index));
        bigz = 1;
    } else {
        dd = Rf_length(_index);
    }
    int d = verify_dimension(dd, k, layout);

    double max;
    if (!bigz) {
        max = fallfact(n, k);
        bigz = max > INT_MAX;
    }

    unsigned int* ap;
    ap = (unsigned int*) R_alloc(k, sizeof(int));

    if (bigz) {
        mpz_t* index;
        gmp_randstate_t randstate;
        mpz_t z;
        mpz_t maxz;
        mpz_init(z);

        if (sampling) {
            GetRNGstate();
            set_gmp_randstate(randstate);
            mpz_init(maxz);
            n_k_permutations_bigz(maxz, n , k);
        } else {
            index = (mpz_t*) R_alloc(d, sizeof(mpz_t));
            for (i = 0; i < d; i++) mpz_init(index[i]);
            int status = as_mpz_array(index, d, _index);
            for(i = 0; i < d; i++) {
                if (status < 0 || mpz_sgn(index[i]) <= 0) {
                    for (i = 0; i < d; i++) mpz_clear(index[i]);
                    Rf_error("expect integer");
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
            identify_k_permutation_bigz(ap, n, k, z);

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
        } else {
            for (i = 0; i < d; i++) mpz_clear(index[i]);
        }

    } else {
        int* index;
        if (sampling) {
            GetRNGstate();
        } else {
            index = as_uint_array(_index);
            for (i = 0; i < d; i++) {
                if (index[0] <= 0) {
                    Rf_error("expect integer");
                }
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                identify_k_permutation(ap, n, k, floor(max * unif_rand())); \
            } else { \
                identify_k_permutation(ap, n, k, index[j] - 1); \
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
