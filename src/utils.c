#include "utils.h"

SEXP resize_row(SEXP x, size_t m, size_t n, size_t d) {
    if (TYPEOF(x) == INTSXP) {
        SEXP y = PROTECT(Rf_allocMatrix(INTSXP, d, n));
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
        SEXP y = PROTECT(Rf_allocMatrix(REALSXP, d, n));
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
        SEXP y = PROTECT(Rf_allocMatrix(STRSXP, d, n));
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
        SEXP y = PROTECT(Rf_allocMatrix(INTSXP, n, d));
        int* yp = INTEGER(y);
        int* xp = INTEGER(x);
        size_t i;
        for (i=0; i<n*d; i++) yp[i] = xp[i];
        UNPROTECT(1);
        return y;
    } else if (TYPEOF(x) == REALSXP) {
        SEXP y = PROTECT(Rf_allocMatrix(REALSXP, n, d));
        double* yp = REAL(y);
        double* xp = REAL(x);
        size_t i;
        for (i=0; i<n*d; i++) yp[i] = xp[i];
        UNPROTECT(1);
        return y;
    } else if (TYPEOF(x) == STRSXP) {
        SEXP y = PROTECT(Rf_allocMatrix(REALSXP, n, d));
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

SEXP resize_layout(SEXP x, size_t d, char layout) {
    SEXP result;
    if (layout == 'r') {
        result = resize_row(x, Rf_nrows(x), Rf_ncols(x), d);
    } else if (layout == 'c') {
        result = resize_col(x, Rf_nrows(x), Rf_ncols(x), d);
    } else {
        result = resize_list(x, Rf_length(x), d);
    }
    return result;
}

void check_factor(SEXP result, SEXP labels) {
    SEXP resulti;
    int result_type = TYPEOF(result);
    int i, n;
    if (Rf_isFactor(labels)) {
        if (result_type == INTSXP || result_type == REALSXP) {
            Rf_setAttrib(result, R_ClassSymbol, Rf_getAttrib(labels, R_ClassSymbol));
            Rf_setAttrib(result, R_LevelsSymbol, Rf_getAttrib(labels, R_LevelsSymbol));
        } else if (result_type == VECSXP) {
            n = Rf_length(result);
            for (i = 0; i < n; i++) {
                resulti = VECTOR_ELT(result, i);
                Rf_setAttrib(resulti, R_ClassSymbol, Rf_getAttrib(labels, R_ClassSymbol));
                Rf_setAttrib(resulti, R_LevelsSymbol, Rf_getAttrib(labels, R_LevelsSymbol));
            }
        }
    }
}

char check_layout(SEXP _layout) {
    char layout;
    if (_layout == R_NilValue) {
        layout = 'r';
    } else {
        layout = CHAR(Rf_asChar(_layout))[0];
        if (layout != 'r' && layout != 'c' && layout != 'l') layout = 'r';
    }
    return layout;
}

int check_dimension(double dd, int k, char layout) {
    if (dd > INT_MAX) Rf_error("too many results");
    if (layout != 'l') {
        if (dd * k > R_XLEN_T_MAX) Rf_error("too many results");
    }
    return round(dd);
}

int variable_exist(SEXP state, char* name, int TYPE, int k, void** p) {
    SEXP v;
    int status = 0;
    if (state == R_NilValue) {
        v = R_UnboundValue;
    } else {
        v = Rf_findVarInFrame(state, Rf_install(name));
    }

    if (v == R_UnboundValue) {
        if (state == R_NilValue) {
            if (TYPE == INTSXP) {
                *p = (int*) R_alloc(k, sizeof(int));
            } else {
                Rf_error("type %d is not yet supported", TYPE);
            }
        } else {
            v = PROTECT(Rf_allocVector(TYPE, k));
            Rf_defineVar(Rf_install(name), v, state);
            UNPROTECT(1);
            if (TYPE == INTSXP) {
                *p = INTEGER(v);
            } else {
                Rf_error("type %d is not yet supported", TYPE);
            }
        }
    } else {
        if (TYPE == INTSXP) {
            *p = INTEGER(v);
        } else {
            Rf_error("type %d is not yet supported", TYPE);
        }
        status = 1;
    }
    return status;
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

double fact(size_t n) {
    double out;
    size_t i;
    out = 1;
    for(i=0; i<n; i++) {
        out = out * (n - i);
    }
    return out;
}

double fallfact(size_t n, size_t k) {
    double out;
    size_t i;
    if (n < k) {
        return 0;
    }
    out = 1;
    for(i=0; i<k; i++) {
        out = out * (n - i);
    }
    return out;
}


double choose(size_t n, size_t k) {
    double out = 1;
    size_t h, i;
    if (n < k) {
        return 0;
    }
    h = 0;
    for (i=1; i<=k; i++){
        h++;
        out = out * h / i;
    }
    for (i=1; i<=n-k; i++){
        h++;
        out = out * h / i;
    }

    return out;
}

double multichoose(int* freq, size_t flen) {
    double out = 1;
    size_t h, i, j;
    h = 0;
    for (i=0; i<flen; i++) {
        for (j=1; j<=freq[i]; j++) {
            h++;
            out = out * h / j;
        }
    }
    return out;
}
