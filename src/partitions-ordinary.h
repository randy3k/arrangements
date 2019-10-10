#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "next.h"
#include "utils.h"
#include "macros.h"
#include "partitions_utils.h"


void nth_asc_partition(unsigned int* ar, unsigned int n, unsigned int index) {
    unsigned int i, j;
    unsigned int start = 1;
    unsigned int sum = n;
    unsigned int count, this_count;

    for (i = 0; i < n; i++) {
        count = 0;
        if (sum > 0 && i < n - 1) {
            for (j = start; j <= n; j++) {
                this_count = count + n_min_partitions(sum - j, j);
                if (this_count > index) {
                    ar[i] = j;
                    start = j;
                    sum -= j;
                    index -= count;
                    break;
                }
                count = this_count;
            }
        } else if (i == n - 1) {
            ar[i] = sum;
        } else {
            ar[i] = 0;
        }
    }
}


void nth_asc_partition_bigz(unsigned int* ar, unsigned int n, mpz_t index) {
    unsigned int i, j;
    unsigned int start = 1;
    unsigned int sum = n;
    mpz_t count, this_count;
    mpz_init(count);
    mpz_init(this_count);

    for (i = 0; i < n; i++) {
        mpz_set_ui(count, 0);
        if (sum > 0 && i < n - 1) {
            for (j = start; j <= n; j++) {
                n_min_partitions_bigz(this_count, sum - j, j);
                mpz_add(this_count, this_count, count);
                if (mpz_cmp(this_count, index) > 0) {
                    ar[i] = j;
                    start = j;
                    sum -= j;
                    mpz_sub(index, index, count);
                    break;
                }
                mpz_set(count, this_count);
            }
        } else if (i == n - 1) {
            ar[i] = sum;
        } else {
            ar[i] = 0;
        }
    }

    mpz_clear(count);
    mpz_clear(this_count);
}


SEXP next_asc_partitions(int n, char layout, int d, SEXP _skip, SEXP state) {
    int i, j, k;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_partitions(n);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;
    int* kp;

    if (!variable_exists(state, "a", INTSXP, n, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            for(i=0; i<n; i++) ap[i] = 1;
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_partitions_bigz(maxz, n);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set(skipz, 0);
                }
                mpz_clear(maxz);
                nth_asc_partition_bigz(ap, n, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_asc_partition(ap, n, skip);
            }
        }
        status = 0;
    }

    if (!variable_exists(state, "k", INTSXP, 1, (void**) &kp)) {
        if (Rf_isNull(_skip)) {
            kp[0] = n - 1;
        } else  {
            for (i = 0; i < n; i++) {
                if (ap[i] == 0) {
                    break;
                }
            }
            kp[0] = i - 1;
        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_asc_partition(ap, kp)) { \
            status = 0; \
            break; \
        } \
        k = kp[0] + 1;

    RESULT_PART();

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    UNPROTECT(nprotect);
    return result;
}


SEXP catch_asc_partitions(int n, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = n_partitions(n);
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
        n_partitions_bigz(maxz, n);

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
            nth_asc_partition_bigz(ap, n, z); \
            for (i = 0; i < n; i++) { \
                if (ap[i] == 0) { \
                    break; \
                } \
            } \
            k = i;

        RESULT_PART();

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
                nth_asc_partition(ap, n, floor(maxd * unif_rand())); \
            } else { \
                nth_asc_partition(ap, n, index[j] - 1); \
            } \
            for (i = 0; i < n; i++) { \
                if (ap[i] == 0) { \
                    break; \
                } \
            } \
            k = i;

        RESULT_PART();

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}


void nth_desc_partition(unsigned int* ar, unsigned int n, unsigned int index) {
    unsigned int i, j;
    unsigned int start = n;
    unsigned int sum = n;
    unsigned int count, this_count;

    for (i = 0; i < n; i++) {
        count = 0;
        if (sum > 0 && i < n - 1) {
            for (j = start; j >= 1; j--) {
                this_count = count + n_max_partitions(sum - j, j);
                if (this_count > index) {
                    ar[i] = j;
                    start = j;
                    sum -= j;
                    index -= count;
                    break;
                }
                count = this_count;
            }
        } else if (i == n - 1) {
            ar[i] = sum;
        } else {
            ar[i] = 0;
        }
    }
}


void nth_desc_partition_bigz(unsigned int* ar, unsigned int n, mpz_t index) {
    unsigned int i, j;
    unsigned int start = n;
    unsigned int sum = n;
    mpz_t count, this_count;
    mpz_init(count);
    mpz_init(this_count);

    for (i = 0; i < n; i++) {
        mpz_set_ui(count, 0);
        if (sum > 0 && i < n - 1) {
            for (j = start; j >= 1; j--) {
                n_max_partitions_bigz(this_count, sum - j, j);
                mpz_add(this_count, this_count, count);
                if (mpz_cmp(this_count, index) > 0) {
                    ar[i] = j;
                    start = j;
                    sum -= j;
                    mpz_sub(index, index, count);
                    break;
                }
                mpz_set(count, this_count);
            }
        } else if (i == n - 1) {
            ar[i] = sum;
        } else {
            ar[i] = 0;
        }
    }

    mpz_clear(count);
    mpz_clear(this_count);
}

SEXP next_desc_partitions(int n, char layout, int d, SEXP _skip, SEXP state) {
    int i, j, k;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_partitions(n);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;
    int* hp;
    int* kp;

    if (!variable_exists(state, "a", INTSXP, n, (void**) &ap)) {
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
                n_partitions_bigz(maxz, n);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set(skipz, 0);
                }
                mpz_clear(maxz);
                nth_desc_partition_bigz(ap, n, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_desc_partition(ap, n, skip);
            }
        }
        status = 0;
    }

    if (!variable_exists(state, "h", INTSXP, 1, (void**) &hp)) {
        if (Rf_isNull(_skip)) {
            hp[0] = 0;
        } else  {
            for (i = 0; i < n; i++) {
                if (ap[i] <= 1) {
                    break;
                }
            }
            hp[0] = i - 1;
        }
        status = 0;
    }

    if (!variable_exists(state, "k", INTSXP, 1, (void**) &kp)) {
        if (Rf_isNull(_skip)) {
            kp[0] = 1;
        } else  {
            for (i = 0; i < n; i++) {
                if (ap[i] == 0) {
                    break;
                }
            }
            kp[0] = i;
            // for the algorithm
            for (j = i; j < n; j++) {
                ap[j] = 1;
            }
        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_desc_partition(ap, hp, kp)) { \
            status = 0; \
            break; \
        } \
        k = kp[0];

    RESULT_PART();

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    UNPROTECT(nprotect);
    return result;
}


SEXP catch_desc_partitions(int n, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = n_partitions(n);
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
        n_partitions_bigz(maxz, n);

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
            nth_desc_partition_bigz(ap, n, z); \
            for (i = 0; i < n; i++) { \
                if (ap[i] == 0) { \
                    break; \
                } \
            } \
            k = i;

        RESULT_PART();

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
                nth_desc_partition(ap, n, floor(maxd * unif_rand())); \
            } else { \
                nth_desc_partition(ap, n, index[j] - 1); \
            } \
            for (i = 0; i < n; i++) { \
                if (ap[i] == 0) { \
                    break; \
                } \
            } \
            k = i;

        RESULT_PART();

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}
