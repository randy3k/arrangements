#include <R.h>
#include <Rinternals.h>
#include "utils.h"

SEXP resize_row(SEXP x, size_t n, size_t m, size_t d) {
    if (TYPEOF(x) == INTSXP) {
        SEXP y = PROTECT(Rf_allocVector(INTSXP, n*d));
        int* yp = INTEGER(y);
        int* xp = INTEGER(x);
        size_t i, j;
        for (j=0; j<d; j++) {
            for (i=0; i<n; i++) {
                yp[j + i*d] = xp[j + i*m];
            }
        }
        UNPROTECT(1);
        return y;
    } else if (TYPEOF(x) == REALSXP) {
        SEXP y = PROTECT(Rf_allocVector(REALSXP, n*d));
        double* yp = REAL(y);
        double* xp = REAL(x);
        size_t i, j;
        for (j=0; j<d; j++) {
            for (i=0; i<n; i++) {
                yp[j + i*d] = xp[j + i*m];
            }
        }
        UNPROTECT(1);
        return y;
    } else if (TYPEOF(x) == STRSXP) {
        SEXP y = PROTECT(Rf_allocVector(STRSXP, n*d));
        size_t i, j;
        for (j=0; j<d; j++) {
            for (i=0; i<n; i++) {
                SET_STRING_ELT(y, j + i*d, STRING_ELT(x, j + i*m));
            }
        }
        UNPROTECT(1);
        return y;
    }
    return R_NilValue;
}

SEXP resize_col(SEXP x, size_t n, size_t m, size_t d) {
    if (TYPEOF(x) == INTSXP) {
        SEXP y = PROTECT(Rf_allocVector(INTSXP, n*d));
        int* yp = INTEGER(y);
        int* xp = INTEGER(x);
        size_t i;
        for (i=0; i<n*d; i++) yp[i] = xp[i];
        UNPROTECT(1);
        return y;
    } else if (TYPEOF(x) == REALSXP) {
        SEXP y = PROTECT(Rf_allocVector(REALSXP, n*d));
        double* yp = REAL(y);
        double* xp = REAL(x);
        size_t i;
        for (i=0; i<n*d; i++) yp[i] = xp[i];
        UNPROTECT(1);
        return y;
    } else if (TYPEOF(x) == STRSXP) {
        SEXP y = PROTECT(Rf_allocVector(REALSXP, n*d));
        size_t i;
        for (i=0; i<n*d; i++) SET_STRING_ELT(y, i, STRING_ELT(x, i));
        UNPROTECT(1);
        return y;
    }
    return R_NilValue;
}

SEXP resize_list(SEXP x, size_t m, size_t d) {
    PROTECT(x);
    SEXP y = PROTECT(Rf_allocVector(VECSXP, d));
    size_t i;
    for(i=0; i<d; i++) {
        SET_VECTOR_ELT(y, i, VECTOR_ELT(x, i));
    }
    UNPROTECT(2);
    return y;
}

int as_uint(SEXP x) {
    double y = Rf_asReal(x);
    int z = (int) y;
    if (y != z || z < 0) Rf_error("expect non-negative integer");
    return z;
}

SEXP as_uint_array(SEXP x) {
    SEXP y;
    size_t i, n;
    double w;
    int z;
    int* x_intp;
    double* x_doublep;
    int* yp;

    if (TYPEOF(x) == INTSXP) {
        n = Rf_length(x);
        x_intp = INTEGER(x);
        for (i=0; i<n; i++) {
            z = x_intp[i];
            if (z < 0) Rf_error("expect non-negative integer");
        }
        return x;
    } else if (TYPEOF(x) == REALSXP) {
        n = Rf_length(x);
        PROTECT(y = Rf_allocVector(INTSXP, n));
        yp = INTEGER(y);
        x_doublep = REAL(x);
        for (i=0; i<n; i++) {
            w = x_doublep[i];
            z = (int) w;
            if (w != z || w < 0) Rf_error("expect non-negative integer");
            yp[i] = (int) x_doublep[i];
        }
        UNPROTECT(1);
        return y;
    }
    return x;
}

double npr(size_t n, size_t r) {
    double out;
    size_t i;
    if (n < r) {
        return 0;
    }
    out = 1;
    for(i=0; i<r; i++) {
        out = out * (n - i);
    }
    return out;
}


double ncr(size_t n, size_t r) {
    double out = 1;
    size_t i, k;
    if (n < r) {
        return 0;
    }
    k = 0;
    for (i=1; i<=r; i++){
        k++;
        out = out * k / i;
    }
    for (i=1; i<=n-r; i++){
        k++;
        out = out * k / i;
    }

    return out;
}

double multichoose(int* f, size_t flen) {
    double out = 1;
    size_t i, j, k;
    k = 0;
    for (i=0; i<flen; i++) {
        for (j=1; j<=f[i]; j++) {
            k++;
            out = out * k / j;
        }
    }
    return out;
}
