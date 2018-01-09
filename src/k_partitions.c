#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "algorithms/k_partition.h"
#include "utils.h"

double npartitions_k(int n, int k);

SEXP next_asc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _type) {
    size_t i, j;

    int n = as_uint(_n);
    int k = as_uint(_k);
    int d;
    double dd;
    if (Rf_asInteger(_d) == -1) {
        dd = npartitions_k(n, k);
    } else {
        dd = as_uint(_d);
    }

    char type;
    if (_type == R_NilValue) {
        type = 'r';
    } else {
        type = CHAR(Rf_asChar(_type))[0];
        if (type != 'r' && type != 'c' && type != 'l') type = 'r';
    }

    if (type == 'l') {
        if (dd > INT_MAX) Rf_error("too many results");
    } else {
        if (dd * k > INT_MAX) Rf_error("too many results");
    }
    d = round(dd);

    SEXP as;
    unsigned int* ap;

    int status = 0;

    if (state == R_NilValue) {
        as = R_UnboundValue;
    } else {
        as = Rf_findVarInFrame(state, Rf_install("a"));
    }
    if (as == R_UnboundValue) {
        if (state == R_NilValue) {
            ap = (unsigned int*) R_alloc(n, sizeof(int));
        } else {
            as = PROTECT(Rf_allocVector(INTSXP, n));
            Rf_defineVar(Rf_install("a"), as, state);
            UNPROTECT(1);
            ap = (unsigned int*) INTEGER(as);
        }
        for(i=0; i<k-1; i++) ap[i] = 1;
        ap[k-1] = n - k + 1;
    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    SEXP rdim;
    SEXP result, resulti;
    int* resultp;
    int nprotect = 0;

    if (type == 'r') {
        result = PROTECT(Rf_allocVector(INTSXP, k*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_k_partition(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<k; i++) {
                resultp[j + i*d] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_row(result, k, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = j;
        INTEGER(rdim)[1] = k;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'c') {
        result = PROTECT(Rf_allocVector(INTSXP, k*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_k_partition(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<k; i++) {
                resultp[j*k + i] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_col(result, k, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = k;
        INTEGER(rdim)[1] = j;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'l') {
        result = PROTECT(Rf_allocVector(VECSXP, d));
        nprotect++;
        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_k_partition(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            resulti = Rf_allocVector(INTSXP, k);
            SET_VECTOR_ELT(result, j, resulti);
            resultp = INTEGER(resulti);
            for (i=0; i<k; i++) {
                resultp[i] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_list(result, d, j));
            nprotect++;
        }
    }

    UNPROTECT(nprotect);
    return result;
}

SEXP next_desc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _type) {
    int n = as_uint(_n);
    int k = as_uint(_k);
    double dd;
    int d = Rf_asInteger(_d);
    if (d == -1) {
        dd = npartitions_k(n, k);
    } else {
        dd = as_uint(_d);
    }

    char type;
    if (_type == R_NilValue) {
        type = 'r';
    } else {
        type = CHAR(Rf_asChar(_type))[0];
        if (type != 'r' && type != 'c' && type != 'l') type = 'r';
    }

    if (type == 'l') {
        if (dd > INT_MAX) Rf_error("too many results");
    } else {
        if (dd * k > INT_MAX) Rf_error("too many results");
    }
    d = round(dd);

    size_t i, j;

    SEXP as;
    unsigned int* ap;

    int status = 0;

    if (state == R_NilValue) {
        as = R_UnboundValue;
    } else {
        as = Rf_findVarInFrame(state, Rf_install("a"));
    }
    if (as == R_UnboundValue) {
        if (state == R_NilValue) {
            ap = (unsigned int*) R_alloc(n, sizeof(int));
        } else {
            as = PROTECT(Rf_allocVector(INTSXP, n));
            Rf_defineVar(Rf_install("a"), as, state);
            UNPROTECT(1);
            ap = (unsigned int*) INTEGER(as);
        }
        for(i=1; i<k; i++) ap[i] = 1;
        ap[0] = n - k + 1;
    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    SEXP rdim;
    SEXP result, resulti;
    int* resultp;
    int nprotect = 0;

    if (type == 'r') {
        result = PROTECT(Rf_allocVector(INTSXP, k*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_desc_k_partition(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<k; i++) {
                resultp[j + i*d] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_row(result, k, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = j;
        INTEGER(rdim)[1] = k;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'c') {
        result = PROTECT(Rf_allocVector(INTSXP, k*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_desc_k_partition(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<k; i++) {
                resultp[j*k + i] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_col(result, k, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = k;
        INTEGER(rdim)[1] = j;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'l') {
        result = PROTECT(Rf_allocVector(VECSXP, d));
        nprotect++;
        for (j=0; j<d; j++) {
            if (status) {
                if (!next_desc_k_partition(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            resulti = Rf_allocVector(INTSXP, k);
            SET_VECTOR_ELT(result, j, resulti);
            resultp = INTEGER(resulti);
            for (i=0; i<k; i++) {
                resultp[i] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_list(result, d, j));
            nprotect++;
        }
    }

    UNPROTECT(nprotect);
    return result;
}

double npartitions_k(int n, int k) {
    if (n < k) {
        return 0;
    } else if (n == 0 && k == 0) {
        return 1;
    } else if (k == 0) {
        return 0;
    }
    int n1 = n-k+1;
    double* p = (double*) malloc(n1*k * sizeof(double));
    int i, j, h;

    for (j=0; j<k; j++) {
        p[j] = 1;
    }
    for (i=1; i<n1; i++) {
        p[i*k] = 1;
        for (j=1; j<k; j++) {
            h = i*k + j;
            if (i > j) {
                p[h] =  p[h - 1] + p[h - (j + 1)*k];
            } else {
                p[h] =  p[h - 1];
            }
        }
    }
    double out = p[n1*k - 1];
    free(p);
    return out;
}


SEXP npart_k(SEXP _n, SEXP _k) {
    int n = as_uint(_n);
    int k = as_uint(_k);
    return Rf_ScalarReal(npartitions_k(n, k));
}


void npartitions_k_bigz(mpz_t z, int n, int k) {
    if (n < k) {
        mpz_set_ui(z, 0);
        return;
    } else if (n == 0 && k == 0) {
        mpz_set_ui(z, 1);
        return;
    } else if (k == 0) {
        mpz_set_ui(z, 0);
        return;
    }

    int n1 = n-k+1;
    int i, j, h;

    mpz_t* p = (mpz_t*) malloc(n1*k * sizeof(mpz_t));
    for (i=0; i<n1*k; i++) mpz_init(p[i]);

    for (j=0; j<k; j++) {
        mpz_set_ui(p[j], 1);
    }
    for (i=1; i<n1; i++) {
        mpz_set_ui(p[i*k], 1);
        for (j=1; j<k; j++) {
            h = i*k + j;
            if (i > j) {
                mpz_add(p[h], p[h - 1], p[h - (j + 1)*k]);
            } else {
                mpz_set(p[h], p[h - 1]);
            }
        }
    }
    mpz_set(z, p[n1*k - 1]);
    for (i=0; i<n1*k; i++) mpz_clear(p[i]);
    free(p);
}

SEXP npart_k_bigz(SEXP _n, SEXP _k) {
    int n = as_uint(_n);
    int k = as_uint(_k);
    mpz_t z;
    mpz_init(z);
    npartitions_k_bigz(z, n, k);
    char* c = mpz_get_str(NULL, 10, z);
    SEXP out = Rf_mkString(c);
    mpz_clear(z);
    free(c);
    return out;
}
