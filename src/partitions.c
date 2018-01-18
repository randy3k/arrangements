#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "algorithms/partition.h"
#include "utils.h"

double npartitions(int n);

SEXP next_asc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _type) {
    size_t i, j, k;

    int n = as_uint(_n);
    int d;
    double dd;
    if (Rf_asInteger(_d) == -1) {
        dd = npartitions(n);
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
        if (dd * n > INT_MAX) Rf_error("too many results");
    }
    d = round(dd);

    SEXP as, ks;
    unsigned int* ap;
    int nprotect = 0;

    int status = 0;

    if (state == R_NilValue) {
        as = R_UnboundValue;
        ks = R_UnboundValue;
    } else {
        as = PROTECT(Rf_findVarInFrame(state, Rf_install("a")));
        ks = PROTECT(Rf_findVarInFrame(state, Rf_install("k")));
        nprotect += 2;
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
        for(i=0; i<n; i++) ap[i] = 1;

    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    k = (ks == R_UnboundValue) ? n-1 : Rf_asInteger(ks);

    SEXP rdim;
    SEXP result, resulti;
    int* resultp;

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

    if (state != R_NilValue) {
        ks = PROTECT(Rf_ScalarInteger((int) k));
        nprotect++;
        Rf_defineVar(Rf_install("k"), ks, state);
    }

    UNPROTECT(nprotect);
    return result;
}


SEXP next_desc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _type) {
    size_t h, k, i, j;

    int n = as_uint(_n);
    double dd;
    int d = Rf_asInteger(_d);
    if (d == -1) {
        dd = npartitions(n);
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
        if (dd * n > INT_MAX) Rf_error("too many results");
    }
    d = round(dd);

    SEXP as, hs, ks;
    unsigned int* ap;
    int nprotect = 0;

    int status = 0;

    if (state == R_NilValue) {
        as = R_UnboundValue;
        hs = R_UnboundValue;
        ks = R_UnboundValue;
    } else {
        as = PROTECT(Rf_findVarInFrame(state, Rf_install("a")));
        hs = PROTECT(Rf_findVarInFrame(state, Rf_install("h")));
        ks = PROTECT(Rf_findVarInFrame(state, Rf_install("k")));
        nprotect += 3;
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

    if (state != R_NilValue) {
        hs = PROTECT(Rf_ScalarInteger((int) h));
        nprotect++;
        ks = PROTECT(Rf_ScalarInteger((int) k));
        nprotect++;
        Rf_defineVar(Rf_install("h"), hs, state);
        Rf_defineVar(Rf_install("k"), ks, state);
    }

    UNPROTECT(nprotect);
    return result;
}

double npartitions(int n) {
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

SEXP npart(SEXP _n) {
    int n = as_uint(_n);
    return Rf_ScalarReal(npartitions(n));
}

void npartitions_bigz(mpz_t z, int n) {
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

SEXP npart_bigz(SEXP _n) {
    int n = as_uint(_n);
    mpz_t z;
    mpz_init(z);
    npartitions_bigz(z, n);
    char* c = mpz_get_str(NULL, 10, z);
    SEXP out = Rf_mkString(c);
    mpz_clear(z);
    free(c);
    return out;
}
