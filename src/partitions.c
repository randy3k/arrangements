#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "algorithms/partition.h"
#include "utils.h"


SEXP next_asc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _type) {
    int n = as_uint(_n);
    int d = as_uint(_d);
    char type = CHAR(Rf_asChar(_type))[0];

    size_t k, i, j;

    SEXP as, ks;
    unsigned int* ap;

    int status = 0;

    if (state == R_NilValue) {
        as = R_UnboundValue;
        ks = R_UnboundValue;
    } else {
        as = Rf_findVarInFrame(state, Rf_install("a"));
        ks = Rf_findVarInFrame(state, Rf_install("k"));
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
        for(i=0; i<n; i++) ap[i] = 1;

    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    k = (ks == R_UnboundValue) ? n-1 : Rf_asInteger(ks);

    SEXP rdim;
    SEXP result, resulti;
    int* resultp;
    int nprotect = 0;

    if (type == 'r') {
        result = PROTECT(Rf_allocVector(INTSXP, n*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_partition(ap, &k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<=k; i++) {
                resultp[j + i*d] = ap[i];
            }
            for (i=k+1; i<n; i++) {
                resultp[j + i*d] = 0;
            }
        }
        if (status == 0) {
            result = PROTECT(resize_row(result, n, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = j;
        INTEGER(rdim)[1] = n;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'c') {
        result = PROTECT(Rf_allocVector(INTSXP, n*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_partition(ap, &k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<=k; i++) {
                resultp[j*n + i] = ap[i];
            }
            for (i=k+1; i<n; i++) {
                resultp[j*n + i] = 0;
            }
        }
        if (status == 0) {
            result = PROTECT(resize_col(result, n, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = n;
        INTEGER(rdim)[1] = j;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'l') {
        result = PROTECT(Rf_allocVector(VECSXP, d));
        nprotect++;
        for (j=0; j<d; j++) {
            if (status) {
                if (!next_asc_partition(ap, &k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            resulti = Rf_allocVector(INTSXP, k+1);
            SET_VECTOR_ELT(result, j, resulti);
            resultp = INTEGER(resulti);

            for (i=0; i<=k; i++) {
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
    } else {
        Rf_defineVar(Rf_install("k"), Rf_ScalarInteger((int) k), state);
    }

    UNPROTECT(nprotect);
    return result;
}


SEXP next_desc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _type) {
    int n = as_uint(_n);
    int d = as_uint(_d);
    char type = CHAR(Rf_asChar(_type))[0];

    size_t h, k, i, j;

    SEXP as, hs, ks;
    unsigned int* ap;

    int status = 0;

    if (state == R_NilValue) {
        as = R_UnboundValue;
        hs = R_UnboundValue;
        ks = R_UnboundValue;
    } else {
        as = Rf_findVarInFrame(state, Rf_install("a"));
        hs = Rf_findVarInFrame(state, Rf_install("h"));
        ks = Rf_findVarInFrame(state, Rf_install("k"));
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
        ap[0] = n;
        for(i=1; i<n; i++) ap[i] = 1;
    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    h = (hs == R_UnboundValue) ? 0 : Rf_asInteger(hs);
    k = (ks == R_UnboundValue) ? 1 : Rf_asInteger(ks);

    SEXP rdim;
    SEXP result, resulti;
    int* resultp;
    int nprotect = 0;

    if (type == 'r') {
        result = PROTECT(Rf_allocVector(INTSXP, n*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if(!next_desc_partition(ap, &h, &k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<k; i++) {
                resultp[j + i*d] = ap[i];
            }
            for (i=k; i<n; i++) {
                resultp[j + i*d] = 0;
            }
        }
        if (status == 0) {
            result = PROTECT(resize_row(result, n, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = j;
        INTEGER(rdim)[1] = n;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'c') {
        result = PROTECT(Rf_allocVector(INTSXP, n*d));
        nprotect++;
        resultp = INTEGER(result);

        for (j=0; j<d; j++) {
            if (status) {
                if(!next_desc_partition(ap, &h, &k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            for (i=0; i<k; i++) {
                resultp[j*n + i] = ap[i];
            }
            for (i=k; i<n; i++) {
                resultp[j*n + i] = 0;
            }
        }
        if (status == 0) {
            result = PROTECT(resize_col(result, n, d, j));
            nprotect++;
        }
        PROTECT(rdim = Rf_allocVector(INTSXP, 2));
        INTEGER(rdim)[0] = n;
        INTEGER(rdim)[1] = j;
        Rf_setAttrib(result, R_DimSymbol, rdim);
        UNPROTECT(1);

    } else if (type == 'l') {
        result = PROTECT(Rf_allocVector(VECSXP, d));
        nprotect++;
        for (j=0; j<d; j++) {
            if (status) {
                if (!next_desc_partition(ap, &h, &k)) {
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

    if (state == R_NilValue) {
        free(ap);
    } else {
        Rf_defineVar(Rf_install("h"), Rf_ScalarInteger((int) h), state);
        Rf_defineVar(Rf_install("k"), Rf_ScalarInteger((int) k), state);
    }

    UNPROTECT(nprotect);
    return result;
}

int _npart(int n) {
    if (n == 0) return 0;
    // find P(1),...,P(n) sequentially
    int i, j, k, s;
    int out;
    int* p = (int*) malloc((n+1) * sizeof(int));
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

SEXP npart(SEXP _n) {
    int n = as_uint(_n);
    return Rf_ScalarInteger(_npart(n));
}

char* _npart_bigz(int n) {
    // find P(1),...,P(n) sequentially
    char* out;
    if (n == 0) {
        out = (char*) malloc(sizeof(char));
        out[0] = '0';
    }
    int i, j, k, s;
    mpz_t* p = (mpz_t*) malloc((n+1) * sizeof(mpz_t));
    for (i=0; i<n+1; i++) mpz_init(p[i]);

    mpz_set_ui(p[0], 1);
    mpz_set_ui(p[1], 1);
    for(i=2 ; i<=n ; i++){
        mpz_set_ui(p[i], 0);
        for (j=1, k=1, s=1; i-j>=0; k+=3, j+=k, s=-s) {
            if (s > 0){
                mpz_add(p[i], p[i], p[i-j]);
            } else {
                mpz_sub(p[i], p[i], p[i-j]);
            }
        }
        for (j=2, k=2, s=1; i-j>=0; k+=3, j+=k, s=-s) {
            if (s > 0){
                mpz_add(p[i], p[i], p[i-j]);
            } else {
                mpz_sub(p[i], p[i], p[i-j]);
            }
        }
    }
    out = mpz_get_str(NULL, 10, p[n]);
    for (i=0; i<n+1; i++) mpz_clear(p[i]);
    free(p);
    return out;
}

SEXP npart_bigz(SEXP _n) {
    int n = as_uint(_n);
    char* c = _npart_bigz(n);
    SEXP out = Rf_mkString(c);
    free(c);
    return out;
}
