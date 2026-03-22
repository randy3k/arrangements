#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "stdlib.h"
#include "../utils.h"
#include "../macros.h"
#include "partitions-utils.h"


unsigned int next_asc_k_partition(unsigned int *ar, size_t n, unsigned int k) {
    // Hindenburg algorithm
    unsigned int a, t;
    unsigned int i, j;
    a = ar[k-1];

    for (j = k-1; j && a - ar[j-1] < 2; j--);

    if (j == 0) {
        for(j=0; j<k-1; j++) ar[j] = 1;
        ar[k-1] = n - k + 1;
        return 0;
    }
    j = j - 1;
    a = ar[j];
    for (i=j; i<k-1; i++) ar[i] = a + 1;

    t = 0;
    for (i=0; i<k-1; i++) t += ar[i];
    ar[k - 1] = n - t;
    return 1;
}

void nth_asc_k_partition(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    double* table = (double*) malloc((n + 1) * (k + 1) * sizeof(double));
    make_k_partition_table(table, n, k);
    nth_asc_k_partition_table(ar, n, k, index, table, k);
    free(table);
}


void nth_asc_k_partition_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    mpz_t* table = (mpz_t*) malloc((n + 1) * (k + 1) * sizeof(mpz_t));
    make_k_partition_table_bigz(table, n, k);
    nth_asc_k_partition_table_bigz(ar, n, k, index, table, k);
    int i;
    for (i = 0; i < (n + 1) * (k + 1); i++) mpz_clear(table[i]);
    free(table);
}


void nth_asc_k_partition_table(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index, double* table, int table_k) {
    unsigned int i, j;
    unsigned int start = 1;
    unsigned int count, this_count;

    for (i = 0; i < k; i++) {
        count = 0;
        unsigned int k_rem = k - i - 1;
        for (j = start; j <= n; j++) {
            double val = 0;
            if (k_rem == 0) {
                val = (n == j);
            } else {
                int n_new = n - j - k_rem * (j - 1);
                if (n_new >= (int) k_rem) val = table[n_new * (table_k + 1) + k_rem];
            }
            this_count = count + val;
            if (this_count > index) {
                ar[i] = j;
                start = j;
                n -= j;
                index -= count;
                break;
            }
            count = this_count;
        }
    }
}


void nth_asc_k_partition_table_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index, mpz_t* table, int table_k) {
    unsigned int i, j;
    unsigned int start = 1;
    mpz_t count, this_count;
    mpz_init(count);
    mpz_init(this_count);

    for (i = 0; i < k; i++) {
        mpz_set_ui(count, 0);
        unsigned int k_rem = k - i - 1;
        for (j = start; j <= n; j++) {
            if (k_rem == 0) {
                mpz_add_ui(this_count, count, n == j);
            } else {
                int n_new = n - j - k_rem * (j - 1);
                if (n_new >= (int) k_rem) {
                    mpz_add(this_count, count, table[n_new * (table_k + 1) + k_rem]);
                } else {
                    mpz_set(this_count, count);
                }
            }
            if (mpz_cmp(this_count, index) > 0) {
                ar[i] = j;
                start = j;
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


SEXP next_asc_k_partitions(int n, int k, char layout, int d, SEXP _skip, SEXP state) {
    int i, j;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_k_partitions(n, k);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;

    if (!variable_exists(state, (char*)"a", INTSXP, k, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            for(i=0; i<k-1; i++) ap[i] = 1;
            ap[k-1] = n - k + 1;
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_k_partitions_bigz(maxz, n, k);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set_ui(skipz, 0);
                }
                mpz_clear(maxz);
                nth_asc_k_partition_bigz(ap, n, k, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_asc_k_partition(ap, n, k, skip);
            }
        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_asc_k_partition(ap, n, k)) { \
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


SEXP draw_asc_k_partitions(int n, int k, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = n_k_partitions(n, k);
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

        mpz_t* table = (mpz_t*) R_alloc((n + 1) * (k + 1), sizeof(mpz_t));
        make_k_partition_table_bigz(table, n, k);
        mpz_set(maxz, table[n * (k + 1) + k]);

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
                    for (i = 0; i < (n+1)*(k+1); i++) mpz_clear(table[i]);
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
            nth_asc_k_partition_table_bigz(ap, n, k, z, table, k);

        RESULT_K_PART();

        mpz_clear(z);
        mpz_clear(maxz);
        for (i = 0; i < (n+1)*(k+1); i++) mpz_clear(table[i]);
        if (sampling){
            gmp_randclear(randstate);
            PutRNGstate();
        } else {
            for (i = 0; i < d; i++) mpz_clear(index[i]);
        }

    } else {
        int* index;
        double* table = (double*) R_alloc((n + 1) * (k + 1), sizeof(double));
        make_k_partition_table(table, n, k);

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
                nth_asc_k_partition_table(ap, n, k, floor(maxd * unif_rand()), table, k); \
            } else { \
                nth_asc_k_partition_table(ap, n, k, index[j] - 1, table, k); \
            }

        RESULT_K_PART();

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}
