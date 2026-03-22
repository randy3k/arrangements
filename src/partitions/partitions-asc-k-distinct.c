#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include <stdlib.h>
#include "../utils.h"
#include "../macros.h"
#include "partitions-utils.h"

void nth_asc_k_partition_table(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index, double* table, int table_k);
void nth_asc_k_partition_table_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index, mpz_t* table, int table_k);
void nth_asc_k_partition(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index);
void nth_asc_k_partition_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index);

void nth_asc_k_distinct_partition(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    double k2 = choose(k, 2);
    if (n < k2) return;
    double* table = (double*) malloc((n - (int)k2 + 1) * (k + 1) * sizeof(double));
    make_k_partition_table(table, n - (int)k2, k);
    nth_asc_k_distinct_partition_table(ar, n, k, index, table, k);
    free(table);
}


void nth_asc_k_distinct_partition_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    double k2 = choose(k, 2);
    if (n < k2) return;
    mpz_t* table = (mpz_t*) malloc((n - (int)k2 + 1) * (k + 1) * sizeof(mpz_t));
    make_k_partition_table_bigz(table, n - (int)k2, k);
    nth_asc_k_distinct_partition_table_bigz(ar, n, k, index, table, k);
    int i;
    for (i = 0; i < (n - (int)k2 + 1) * (k + 1); i++) mpz_clear(table[i]);
    free(table);
}


void nth_asc_k_distinct_partition_table(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index, double* table, int table_k) {
    int i;
    double k2 = choose(k, 2);
    if (n < k2) return;
    nth_asc_k_partition_table(ar, n - k2, k, index, table, table_k);
    for (i = 1; i < k; i++) ar[i] += i;
}


void nth_asc_k_distinct_partition_table_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index, mpz_t* table, int table_k) {
    int i;
    double k2 = choose(k, 2);
    if (n < k2) return;
    nth_asc_k_partition_table_bigz(ar, n - k2, k, index, table, table_k);
    for (i = 1; i < k; i++) ar[i] += i;
}


SEXP next_asc_k_distinct_partitions(int n, int k, char layout, int d, SEXP _skip, SEXP state) {
    int i, j;
    int nprotect = 0;
    unsigned int* ap;
    SEXP result = R_NilValue;
    double k2 = choose(k, 2);
    int status = 1;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_k_distinct_partitions(n, k);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, k, layout);

    if (!variable_exists(state, (char*)"a", INTSXP, k, (void**) &ap)) {
        if (n < k2) {
             if (layout == 'r' || layout == 'c') {
                 result = Rf_allocMatrix(INTSXP, 0, k);
             } else {
                 result = Rf_allocVector(VECSXP, 0);
             }
             return result;
        }

        if (Rf_isNull(_skip)) {
            nth_asc_k_distinct_partition(ap, n, k, 0);
        } else {
            if (bigz) {
                mpz_t maxz, skipz;
                mpz_init(maxz);
                mpz_init(skipz);
                n_k_distinct_partitions_bigz(maxz, n, k);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set_ui(skipz, 0);
                }
                mpz_clear(maxz);
                nth_asc_k_distinct_partition_bigz(ap, n, k, skipz);
                mpz_clear(skipz);
            } else {
                int skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_asc_k_distinct_partition(ap, n, k, skip);
            }
        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else { \
            for (i = 1; i < k; i++) ap[i] -= i; \
            if (next_asc_k_partition(ap, n - (int)k2, k) == 0) { \
                for (i = 1; i < k; i++) ap[i] += i; \
                status = 0; \
                break; \
            } \
            for (i = 1; i < k; i++) ap[i] += i; \
        }

    RESULT_K_PART();

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }

    UNPROTECT(nprotect);
    return result;
}


SEXP draw_asc_k_distinct_partitions(int n, int k, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = n_k_distinct_partitions(n, k);
        bigz = maxd > INT_MAX;
    }

    unsigned int* ap;
    ap = (unsigned int*) R_alloc(k, sizeof(int));

    double k2 = choose(k, 2);

    if (bigz) {
        mpz_t* index;
        gmp_randstate_t randstate;
        mpz_t z;
        mpz_t maxz;
        mpz_init(z);
        mpz_init(maxz);

        if (n < k2) {
            mpz_set_ui(maxz, 0);
            if (!sampling) {
                index = (mpz_t*) R_alloc(d, sizeof(mpz_t));
                for (i = 0; i < d; i++) mpz_init(index[i]);
                as_mpz_array(index, d, _index);
                for(i = 0; i < d; i++) {
                    if (mpz_sgn(index[i]) <= 0 || mpz_cmp(index[i], maxz) > 0) {
                        for (i = 0; i < d; i++) mpz_clear(index[i]);
                        mpz_clear(maxz);
                        mpz_clear(z);
                        Rf_error("invalid index");
                    }
                }
            }
        } else {
            mpz_t* table = (mpz_t*) R_alloc((n - (int)k2 + 1) * (k + 1), sizeof(mpz_t));
            make_k_partition_table_bigz(table, n - (int)k2, k);
            mpz_set(maxz, table[(n - (int)k2) * (k + 1) + k]);

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
                        for (i = 0; i < (n - (int)k2 + 1) * (k + 1); i++) mpz_clear(table[i]);
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
                nth_asc_k_distinct_partition_table_bigz(ap, n, k, z, table, k);

            RESULT_K_PART();

            for (i = 0; i < (n - (int)k2 + 1) * (k + 1); i++) mpz_clear(table[i]);
            if (sampling){
                gmp_randclear(randstate);
                PutRNGstate();
            } else {
                for (i = 0; i < d; i++) mpz_clear(index[i]);
            }
        }
        mpz_clear(z);
        mpz_clear(maxz);

    } else {
        if (n < k2) {
             if (layout == 'r' || layout == 'c') {
                 result = Rf_allocMatrix(INTSXP, d, k);
             } else {
                 result = Rf_allocVector(VECSXP, d);
             }
        } else {
            int* index;
            double* table = (double*) R_alloc((n - (int)k2 + 1) * (k + 1), sizeof(double));
            make_k_partition_table(table, n - (int)k2, k);

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
                    nth_asc_k_distinct_partition_table(ap, n, k, floor(maxd * unif_rand()), table, k); \
                } else { \
                    nth_asc_k_distinct_partition_table(ap, n, k, index[j] - 1, table, k); \
                }

            RESULT_K_PART();

            if (sampling){
                PutRNGstate();
            }
        }
    }

    UNPROTECT(nprotect);
    return result;
}
