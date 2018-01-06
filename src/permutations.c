#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "algorithms/permutation.h"
#include "utils.h"


SEXP next_permutations(SEXP _n, SEXP _d, SEXP state, SEXP labels, SEXP f, SEXP _type) {
    size_t n = Rf_asInteger(_n);
    int d = Rf_asInteger(_d);

    int ltype = TYPEOF(labels);
    int* labels_intp;
    double* labels_doublep;

    int* fp;

    char type = CHAR(Rf_asChar(_type))[0];

    size_t i, j, k;

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
        if (f == R_NilValue) {
            for(i=0; i<n; i++) ap[i] = i + 1;
        } else {
            fp = INTEGER(f);
            k = 0;
            for (i = 0; i< Rf_length(f); i++) {
                for (j = 0; j< fp[i]; j++) {
                    ap[k++] = i + 1;
                }
            }
        }

    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    SEXP result, resulti;
    int* result_intp;
    double* result_doublep;
    int nprotect = 0;

    if (type == 'r') {
        if (labels == R_NilValue) {
            result = PROTECT(Rf_allocVector(INTSXP, n*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(ltype, n*d));
            nprotect++;
            if (ltype == INTSXP) {
                result_intp = INTEGER(result);
                labels_intp = INTEGER(labels);
            } else if (ltype == REALSXP) {
                result_doublep = REAL(result);
                labels_doublep = REAL(labels);
            }
        }

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_permutation(ap, n)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                for (i=0; i<n; i++) {
                    result_intp[j + i*d] = ap[i];
                }
            } else if (ltype == INTSXP) {
                for (i=0; i<n; i++) {
                    result_intp[j + i*d] = labels_intp[ap[i] - 1];
                }
            } else if (ltype == REALSXP) {
                for (i=0; i<n; i++) {
                    result_doublep[j + i*d] = labels_doublep[ap[i] - 1];
                }
            } else if (ltype == STRSXP) {
                for (i=0; i<n; i++) {
                    SET_STRING_ELT(result, j + i*d, STRING_ELT(labels, ap[i] - 1));
                }
            }
        }
        if (status == 0) {
            result = PROTECT(resize_row(result, n, d, j));
            nprotect++;
        }
    } else if (type == 'c') {
        if (labels == R_NilValue) {
            result = PROTECT(Rf_allocVector(INTSXP, n*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(ltype, n*d));
            nprotect++;
            if (ltype == INTSXP) {
                result_intp = INTEGER(result);
                labels_intp = INTEGER(labels);
            } else if (ltype == REALSXP) {
                result_doublep = REAL(result);
                labels_doublep = REAL(labels);
            }
        }

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_permutation(ap, n)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                for (i=0; i<n; i++) {
                    result_intp[j * n + i] = ap[i];
                }
            } else if (ltype == INTSXP) {
                for (i=0; i<n; i++) {
                    result_intp[j * n + i] = labels_intp[ap[i] - 1];
                }
            } else if (ltype == REALSXP) {
                for (i=0; i<n; i++) {
                    result_doublep[j * n + i] = labels_doublep[ap[i] - 1];
                }
            } else if (ltype == STRSXP) {
                for (i=0; i<n; i++) {
                    SET_STRING_ELT(result, j * n + i, STRING_ELT(labels, ap[i] - 1));
                }
            }
        }
        if (status == 0) {
            result = PROTECT(resize_col(result, n, d, j));
            nprotect++;
        }
    } else {
        // type == "list"
        result = PROTECT(Rf_allocVector(VECSXP, d));
        if (ltype == INTSXP) {
            labels_intp = INTEGER(labels);
        } else if (ltype == REALSXP) {
            labels_doublep = REAL(labels);
        }

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_permutation(ap, n)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                resulti = Rf_allocVector(INTSXP, n);
                result_intp = INTEGER(resulti);
                for (i=0; i<n; i++) {
                    result_intp[i] = ap[i];
                }
            } else if (ltype == INTSXP) {
                resulti = Rf_allocVector(INTSXP, n);
                result_intp = INTEGER(resulti);
                for (i=0; i<n; i++) {
                    result_intp[i] = labels_intp[ap[i] - 1];
                }
            } else if (ltype == REALSXP) {
                resulti = Rf_allocVector(REALSXP, n);
                result_doublep = REAL(resulti);
                for (i=0; i<n; i++) {
                    result_doublep[i] = labels_doublep[ap[i] - 1];
                }
            } else if (ltype == STRSXP) {
                resulti = Rf_allocVector(STRSXP, n);
                for (i=0; i<n; i++) {
                    SET_STRING_ELT(resulti, i, STRING_ELT(labels, ap[i] - 1));
                }
            }
            SET_VECTOR_ELT(result, j, resulti);
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

double _multichoose(int* f, size_t flen) {
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

SEXP multichoose(SEXP f) {
    int* fp = INTEGER(f);
    size_t flen = Rf_length(f);
    return Rf_ScalarReal(_multichoose(fp, flen));
}

char* _multichoose_bigz(int* f, size_t flen) {
    mpz_t z;
    mpz_init(z);
    mpz_set_ui(z, 1);

    size_t i, j, k;
    k = 0;
    for (i=0; i<flen; i++) {
        for (j=1; j<=f[i]; j++) {
            k++;
            mpz_mul_ui(z, z, k);
            mpz_cdiv_q_ui(z, z, j);
        }
    }
    char* out = mpz_get_str(NULL, 10, z);
    mpz_clear(z);
    return out;
}

SEXP multichoose_bigz(SEXP f) {
    int* fp = INTEGER(f);
    size_t flen = Rf_length(f);
    char* c = _multichoose_bigz(fp, flen);
    SEXP out = Rf_mkString(c);
    free(c);
    return out;
}
