#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "stdlib.h"
#include "../macros.h"
#include "../math.h"
#include "../utils.h"


// mirror of python itertools.permutations
// https://docs.python.org/3/library/itertools.html#itertools.permutations
unsigned int next_k_permutation(unsigned int *ar, unsigned int *cycle, size_t n, size_t k) {
    // ar = [0, 1, ..., n], cycle = [n, n-1, ..., n-r+1]
    unsigned int i, j, temp;

    for (i = k-1; ; i--) {
        cycle[i] -= 1;
        if (cycle[i] == 0) {
            temp = ar[i];
            for (j = i; j <= n - 2; j++) {
                ar[j] = ar[j + 1];
            }
            ar[n - 1] = temp;
            cycle[i] = n - i;
        } else {
            swap(ar, i, n - cycle[i]);
            return 1;
        }
        if (i == 0) {
            return 0;
        }
    }
}


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


void nth_k_permutation(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    unsigned int i, j;

    for (i = 0; i < k; i++) {
        j = fallfact(n - 1 - i, k - 1 - i);
        ar[i] = index / j;
        index = index % j;
    }

    if (k > 0) {
        for (i = k - 1; i > 0; i--) {
            j = i;
            while (j-- > 0) {
                if (ar[j] <= ar[i]) {
                    ar[i]++;
                }
            }
        }
    }
}

void nth_k_permutation_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
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

    if (k > 0) {
        for (i = k - 1; i > 0; i--) {
            j = i;
            while (j-- > 0) {
                if (ar[j] <= ar[i]) {
                    ar[i]++;
                }
            }
        }
    }

    mpz_clear(q);
    mpz_clear(p);
}


SEXP next_k_permutations(int n, int k, SEXP labels, char layout, int d, SEXP _skip, SEXP state) {
    int i, j;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = fallfact(n, k);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;
    unsigned int* cyclep;

    if (!variable_exists(state, (char*)"a", INTSXP, n, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            for(i=0; i<n; i++) ap[i] = i;
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_k_permutations_bigz(maxz, n, k);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set_ui(skipz, 0);
                }
                mpz_clear(maxz);
                nth_k_permutation_bigz(ap, n, k, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_k_permutation(ap, n, k, skip);
            }
            int* count = (int*) malloc(n * sizeof(int));
            for(i = 0; i < n; i++) count[i] = 1;
            for(i = 0; i < k; i++) count[ap[i]] = 0;
            j = 0;
            for (i = k; i < n; i++) {
                while (count[j] == 0) j++;
                ap[i] = j++;
            }
            free(count);
        }
        status = 0;
    }
    if (!variable_exists(state, (char*)"cycle", INTSXP, k, (void**) &cyclep)) {
        if (Rf_isNull(_skip)) {
            for(i=0; i<k; i++) cyclep[i] = n - i;;
        } else {
            for(i=0; i<k; i++) {
                cyclep[i] = n - ap[i];
                for (j = 0; j < i; j++) {
                    if (ap[j] > ap[i]) {
                        cyclep[i]--;
                    }
                }
            }
        }
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
    UNPROTECT(nprotect);
    return result;
}


SEXP draw_k_permutations(int n, int k, SEXP labels, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = fallfact(n, k);
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
        n_k_permutations_bigz(maxz, n , k);

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
            nth_k_permutation_bigz(ap, n, k, z);

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
                nth_k_permutation(ap, n, k, floor(maxd * unif_rand())); \
            } else { \
                nth_k_permutation(ap, n, k, index[j] - 1); \
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

    UNPROTECT(nprotect);
    return result;
}
