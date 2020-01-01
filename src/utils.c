#include "utils.h"

void swap(unsigned int *ar, unsigned int first, unsigned int second) {
    unsigned int temp = ar[first];
    ar[first] = ar[second];
    ar[second] = temp;
}

void reverse(unsigned int *ar, size_t len) {
    unsigned int i, j;

    for (i = 0, j = len - 1; i < j; i++, j--) {
        swap(ar, i, j);
    }
}

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
    if (dd <= 0) Rf_error("d should be positive");
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

int* as_uint_index(SEXP x) {
    // a clone of as_uint_array for positive integers
    size_t i, n;
    int z;

    if (TYPEOF(x) == INTSXP) {
        int* xp = INTEGER(x);
        n = Rf_length(x);
        for (i=0; i<n; i++) {
            z = xp[i];
            if (z <= 0) Rf_error("invalid index");
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
            if (w != z || w <= 0) Rf_error("invalid index");
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
            if (w != z || w <= 0) Rf_error("invalid index");
            yp[i] = z;
        }
        return yp;
    }
    Rf_error("invalid index");
    return NULL;
}


static int raw_size(mpz_t z) {
    int numb = 8 * sizeof(int);
    return sizeof(int) * (2 + (mpz_sizeinbase(z, 2) + numb - 1) / numb);
}

// from biginteger::as_raw
// https://github.com/cran/gmp/blob/cb2935d6a7948e85e42b1402c6f5a3450547f0cf/src/biginteger.cc#L75
SEXP mpz_to_bigz1(mpz_t z) {
    int size = raw_size(z);
    SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size + sizeof(int)));
    unsigned char* raw = RAW(ans);
    int* r = (int*) raw;
    r[0] = 1; // scalar RAWSXP
    r[1] = size / sizeof(int) - 2;
    r[2] = (int) mpz_sgn(z);
    if (mpz_sgn(z) == 0) {
        // mpz_export writes nothing when z = 0
        r[3] = 0;
    } else {
        mpz_export(&r[3], 0, 1, sizeof(int), 0, 0, z);
    }
    Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));
    UNPROTECT(1);
    return ans;
}

// from biginteger::biginteger
// https://github.com/cran/gmp/blob/cb2935d6a7948e85e42b1402c6f5a3450547f0cf/src/biginteger.cc#L25
int as_mpz_array(mpz_t* a, size_t n, SEXP x) {
    size_t i;
    if (TYPEOF(x) == RAWSXP && Rf_inherits(x, "bigz")) {
        int* v = ((int*) RAW(x)) + 1;
        for (i = 0; i < n; i++) {
            if (v[0] > 0) {
                mpz_import(a[i], v[0], 1, sizeof(int) , 0, 0, &(v[2]));
                if(v[1] == -1) {
                    mpz_neg(a[i], a[i]);
                }
                v = v + v[0] + 2;
            } else {
                mpz_set_ui(a[i], 0);
                v = v + 1;
            }
        }
        return 0;
    } else if (TYPEOF(x) == INTSXP) {
        int* xp = INTEGER(x);
        for (i = 0; i < n; i++) {
            mpz_set_ui(a[i], abs(xp[i]));
            if (xp[i] < 0) {
                mpz_neg(a[i], a[i]);
            }
        }
        return 0;
    } else if (TYPEOF(x) == REALSXP) {
        double* xp = REAL(x);
        int w;
        for (i = 0; i < n; i++) {
            w = (int) fabs(xp[i]);
            if (w != xp[i]) return -1;
            mpz_set_ui(a[i], w);
            if (xp[i] < 0) {
                mpz_neg(a[i], a[i]);
            }
        }
        return 0;
    } else if (TYPEOF(x) == STRSXP) {
        for (i = 0; i < n; i++) {
            if (mpz_set_str(a[i], CHAR(STRING_ELT(x, i)), 10) < 0) {
                return -1;
            }
        }
        return 0;
    }
    return -1;
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
    PROTECT(seeds);

    unsigned int* seedsp = (unsigned int*) INTEGER(seeds);
    mpz_set_ui(z, round(INT_MAX * unif_rand()));
    for (i = 0; i < Rf_length(seeds); i++) {
        mpz_add_ui(z, z, seedsp[i] << 16);
    }
    gmp_randinit_mt(randstate);
    gmp_randseed(randstate, z);
    mpz_clear(z);
    UNPROTECT(1);
}

int index_length(SEXP _index) {
    if (TYPEOF(_index) == RAWSXP && Rf_inherits(_index, "bigz")) {
        return *((int* ) RAW(_index));
    } else {
        return Rf_length(_index);
    }
}
