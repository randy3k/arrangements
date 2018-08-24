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

void attach_factor_levels(SEXP result, SEXP labels) {
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

char layout_flag(SEXP _layout) {
    char layout;
    if (_layout == R_NilValue) {
        layout = 'r';
    } else {
        layout = CHAR(Rf_asChar(_layout))[0];
        if (layout != 'r' && layout != 'c' && layout != 'l') layout = 'r';
    }
    return layout;
}

int verify_dimension(double dd, int k, char layout) {
    if (dd < 0) Rf_error("d should be positive");
    if (dd >= INT_MAX) Rf_error("too many results");
    if (layout != 'l') {
        if (dd * k >= R_XLEN_T_MAX) Rf_error("too many results");
    }
    return round(dd);
}

int variable_exists(SEXP state, char* name, int TYPE, int k, void** p) {
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
    if (y != z || z < 0) Rf_error("expect integer");
    return z;
}

int* as_uint_array(SEXP x) {
    size_t i, n;
    int z;

    if (TYPEOF(x) == INTSXP) {
        int* xp = INTEGER(x);
        n = Rf_length(x);
        for (i=0; i<n; i++) {
            z = xp[i];
            if (z < 0) Rf_error("expect integer");
        }
        return xp;
    } else if (TYPEOF(x) == REALSXP) {
        int* yp;
        double* xp;
        double w;
        n = Rf_length(x);
        yp = (int*) R_alloc(n, sizeof(int));
        xp = REAL(x);
        for (i=0; i<n; i++) {
            w = xp[i];
            z = (int) w;
            if (w != z || w < 0) Rf_error("expect integer");
            yp[i] = z;
        }
        return yp;
    } else if (TYPEOF(x) == STRSXP) {
        int* yp;
        double w;
        n = Rf_length(x);
        yp = (int*) R_alloc(n, sizeof(int));
        for (i=0; i<n; i++) {
            w = atof(CHAR(STRING_ELT(x, i)));
            z = (int) w;
            if (w != z || w < 0) Rf_error("expect integer");
            yp[i] = z;
        }
        return yp;
    }
    Rf_error("expect integer");
    return NULL;
}

double fact(int n) {
    double out;
    int i;
    out = 1;
    for(i=0; i<n; i++) {
        out = out * (n - i);
    }
    return out;
}

double fallfact(int n, int k) {
    double out;
    int i;
    if (n < k) {
        return 0;
    }
    out = 1;
    for(i=0; i<k; i++) {
        out = out * (n - i);
    }
    return out;
}


double choose(int n, int k) {
    double out = 1;
    int i;
    if (n >= 0 && n < k) {
        return 0;
    } else if (k == 0) {
        return 1;
    }

    if (k > n / 2) k = n - k;

    for (i = 1; i <= k; i++) {
        out *= n - k + i;
        out /= i;
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

// from RNG.c
static SEXP GetSeedsFromVar(void)
{
    SEXP seeds = Rf_findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
    if (TYPEOF(seeds) == PROMSXP)
    seeds = Rf_eval(R_SeedsSymbol, R_GlobalEnv);
    return seeds;
}

void set_gmp_randstate(gmp_randstate_t randstate) {
    int i;
    mpz_t z;
    mpz_init(z);

    SEXP seeds = GetSeedsFromVar();
    if (seeds == R_UnboundValue) {
        PutRNGstate();
        seeds = GetSeedsFromVar();
    }

    unsigned int* seedsp = (unsigned int*) INTEGER(seeds);
    mpz_set_ui(z, round(INT_MAX * unif_rand()));
    for (i = 0; i < Rf_length(seeds); i++) {
        mpz_add_ui(z, z, seedsp[i] << 16);
    }
    gmp_randinit_mt(randstate);
    gmp_randseed(randstate, z);
    mpz_clear(z);
}

int index_length(SEXP _index) {
    if (TYPEOF(_index) == RAWSXP && Rf_inherits(_index, "bigz")) {
        return *((int* ) RAW(_index));
    } else {
        return Rf_length(_index);
    }
}
