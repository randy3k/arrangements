#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/partition.h"
#include "utils.h"
#include "macros.h"

double n_partitions(int n);

SEXP next_asc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _layout) {
    int i, j, k;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    int n = as_uint(_n);
    char layout = layout_flag(_layout);

    double dd = Rf_asInteger(_d) == -1 ? n_partitions(n) : as_uint(_d);
    int d = verify_dimension(dd, n, layout);

    unsigned int* ap;
    size_t* kp;

    if (!variable_exist(state, "a", INTSXP, n, (void**) &ap)) {
        for(i=0; i<n; i++) ap[i] = 1;
        status = 0;
    }

    if (!variable_exist(state, "k", INTSXP, 1, (void**) &kp)) {
        kp[0] = n - 1;
        status = 0;
    }

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

#undef NEXT

SEXP next_desc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _layout) {
    int i, j, k;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    int n = as_uint(_n);
    char layout = layout_flag(_layout);

    double dd = Rf_asInteger(_d) == -1 ? n_partitions(n) : as_uint(_d);
    int d = verify_dimension(dd, n, layout);

    unsigned int* ap;
    size_t* hp;
    size_t* kp;

    if (!variable_exist(state, "a", INTSXP, n, (void**) &ap)) {
        ap[0] = n;
        for(i=1; i<n; i++) ap[i] = 1;
        status = 0;
    }

    if (!variable_exist(state, "h", INTSXP, 1, (void**) &hp)) {
        hp[0] = 0;
        status = 0;
    }

    if (!variable_exist(state, "k", INTSXP, 1, (void**) &kp)) {
        kp[0] = 1;
        status = 0;
    }

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

double n_partitions(int n) {
    if (n == 0) return 1;
    // find P(1),...,P(n) sequentially
    int i, j, k, s;
    double out;
    double* p = (double*) malloc((n+1) * sizeof(double));
    p[0] = p[1] = 1;
    for(i=2 ; i<=n ; i++){
        p[i] = 0;
        for (j=1, k=1, s=1; i-j>=0; k+=3, j+=k, s=-s) {
            p[i] += s*p[i-j];
        }
        for (j=2, k=2, s=1; i-j>=0; k+=3, j+=k, s=-s) {
            p[i] += s*p[i-j];
        }
    }
    out = p[n];
    free(p);
    return out;
}

SEXP num_partitions(SEXP _n) {
    int n = as_uint(_n);
    return Rf_ScalarReal(n_partitions(n));
}

void n_partitions_bigz(mpz_t z, int n) {
    // find P(1),...,P(n) sequentially
    if (n == 0) {
        mpz_set_ui(z, 1);
        return;
    }
    int i, j, h, s;
    mpz_t* p = (mpz_t*) malloc((n+1) * sizeof(mpz_t));
    for (i=0; i<n+1; i++) mpz_init(p[i]);

    mpz_set_ui(p[0], 1);
    mpz_set_ui(p[1], 1);
    for(i=2 ; i<=n ; i++){
        for (j=1, h=1, s=1; i-j>=0; h+=3, j+=h, s=-s) {
            if (s > 0){
                mpz_add(p[i], p[i], p[i-j]);
            } else {
                mpz_sub(p[i], p[i], p[i-j]);
            }
        }
        for (j=2, h=2, s=1; i-j>=0; h+=3, j+=h, s=-s) {
            if (s > 0){
                mpz_add(p[i], p[i], p[i-j]);
            } else {
                mpz_sub(p[i], p[i], p[i-j]);
            }
        }
    }
    mpz_set(z, p[n]);
    for (i=0; i<n+1; i++) mpz_clear(p[i]);
    free(p);
}

SEXP num_partitions_bigz(SEXP _n) {
    int n = as_uint(_n);
    mpz_t z;
    mpz_init(z);
    n_partitions_bigz(z, n);
    char* c = mpz_get_str(NULL, 10, z);
    SEXP out = Rf_mkString(c);
    mpz_clear(z);
    free(c);
    return out;
}
