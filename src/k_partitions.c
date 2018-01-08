#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "algorithms/k_partition.h"
#include "utils.h"


SEXP next_asc_k_partitions(SEXP _n, SEXP _m, SEXP _d, SEXP state, SEXP _type) {
    int n = as_uint(_n);
    int m = as_uint(_m);
    int d = as_uint(_d);
    char type = CHAR(Rf_asChar(_type))[0];

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
            ap = (unsigned int*) malloc(n * sizeof(int));
        } else {
            as = PROTECT(Rf_allocVector(INTSXP, n));
            Rf_defineVar(Rf_install("a"), as, state);
            UNPROTECT(1);
            ap = (unsigned int*) INTEGER(as);
        }
        for(i=0; i<m-1; i++) ap[i] = 1;
        ap[m-1] = n - m + 1;
    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    SEXP rdim;
    SEXP result, resulti;
    int* resultp;
    int nprotect = 0;

    if (type == 'r') {
        result = PROTECT(Rf_allocVector(INTSXP, m*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_k_partition(ap, n, m)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<m; i++) {
                resultp[j + i*d] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_row(result, m, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = j;
        INTEGER(rdim)[1] = m;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'c') {
        result = PROTECT(Rf_allocVector(INTSXP, m*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_k_partition(ap, n, m)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<m; i++) {
                resultp[j*m + i] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_col(result, m, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = m;
        INTEGER(rdim)[1] = j;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'l') {
        result = PROTECT(Rf_allocVector(VECSXP, d));
        nprotect++;
        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_k_partition(ap, n, m)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            resulti = Rf_allocVector(INTSXP, m);
            SET_VECTOR_ELT(result, j, resulti);
            resultp = INTEGER(resulti);
            for (i=0; i<m; i++) {
                resultp[i] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_list(result, d, j));
            nprotect++;
        }
    }

    if (state == R_NilValue) {
        free(ap);
    }

    UNPROTECT(nprotect);
    return result;
}

SEXP next_desc_k_partitions(SEXP _n, SEXP _m, SEXP _d, SEXP state, SEXP _type) {
    int n = as_uint(_n);
    int m = as_uint(_m);
    int d = as_uint(_d);
    char type = CHAR(Rf_asChar(_type))[0];

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
            ap = (unsigned int*) malloc(n * sizeof(int));
        } else {
            as = PROTECT(Rf_allocVector(INTSXP, n));
            Rf_defineVar(Rf_install("a"), as, state);
            UNPROTECT(1);
            ap = (unsigned int*) INTEGER(as);
        }
        for(i=1; i<m; i++) ap[i] = 1;
        ap[0] = n - m + 1;
    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    SEXP rdim;
    SEXP result, resulti;
    int* resultp;
    int nprotect = 0;

    if (type == 'r') {
        result = PROTECT(Rf_allocVector(INTSXP, m*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_desc_k_partition(ap, n, m)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<m; i++) {
                resultp[j + i*d] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_row(result, m, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = j;
        INTEGER(rdim)[1] = m;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'c') {
        result = PROTECT(Rf_allocVector(INTSXP, m*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_desc_k_partition(ap, n, m)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<m; i++) {
                resultp[j*m + i] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_col(result, m, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = m;
        INTEGER(rdim)[1] = j;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'l') {
        result = PROTECT(Rf_allocVector(VECSXP, d));
        nprotect++;
        for (j=0; j<d; j++) {
            if (status) {
                if (!next_desc_k_partition(ap, n, m)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            resulti = Rf_allocVector(INTSXP, m);
            SET_VECTOR_ELT(result, j, resulti);
            resultp = INTEGER(resulti);
            for (i=0; i<m; i++) {
                resultp[i] = ap[i];
            }
        }
        if (status == 0) {
            result = PROTECT(resize_list(result, d, j));
            nprotect++;
        }
    }

    if (state == R_NilValue) {
        free(ap);
    }

    UNPROTECT(nprotect);
    return result;
}

int _nfixedpart(int n, int m) {
    if (n < m) return 0;
    int n1 = n-m+1;
    int* p = (int*) malloc(n1*m * sizeof(int));
    int i, j, k;

    for (j=0; j<m; j++) {
        p[j] = 1;
    }
    for (i=1; i<n1; i++) {
        p[i*m] = 1;
        for (j=1; j<m; j++) {
            k = i*m + j;
            if (i > j) {
                p[k] =  p[k - 1] + p[k - (j + 1)*m];
            } else {
                p[k] =  p[k - 1];
            }
        }
    }
    int out = p[n1*m - 1];
    free(p);
    return out;
}


SEXP nfixedpart(SEXP _n, SEXP _m) {
    int n = as_uint(_n);
    int m = as_uint(_m);
    return Rf_ScalarInteger(_nfixedpart(n, m));
}


char* _nfixedpart_bigz(int n, int m) {
    char* out;
    if (n < m) {
        out = (char*) malloc(sizeof(char));
        out[0] = '0';
    }

    int n1 = n-m+1;
    int i, j, k;

    mpz_t* p = (mpz_t*) malloc(n1*m * sizeof(mpz_t));
    for (i=0; i<n1*m; i++) mpz_init(p[i]);

    for (j=0; j<m; j++) {
        mpz_set_ui(p[j], 1);
    }
    for (i=1; i<n1; i++) {
        mpz_set_ui(p[i*m], 1);
        for (j=1; j<m; j++) {
            k = i*m + j;
            if (i > j) {
                mpz_add(p[k], p[k - 1], p[k - (j + 1)*m]);
            } else {
                mpz_set(p[k], p[k - 1]);
            }
        }
    }
    out = mpz_get_str(NULL, 10, p[n1*m - 1]);
    for (i=0; i<n1*m; i++) mpz_clear(p[i]);
    free(p);
    return out;
}

SEXP nfixedpart_bigz(SEXP _n, SEXP _m) {
    int n = as_uint(_n);
    int m = as_uint(_m);
    char* c = _nfixedpart_bigz(n, m);
    SEXP out = Rf_mkString(c);
    free(c);
    return out;
}
