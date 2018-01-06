#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "algorithms/cartesian_product.h"
#include "utils.h"


SEXP next_replace_permutations(SEXP _n, SEXP _r, SEXP _d, SEXP state, SEXP labels, SEXP _type) {
    size_t n = Rf_asInteger(_n);
    size_t r = Rf_asInteger(_r);
    int d = Rf_asInteger(_d);

    int ltype = TYPEOF(labels);
    int* labels_intp;
    double* labels_doublep;

    int* fp;

    char type = CHAR(Rf_asChar(_type))[0];

    size_t i, j, k;

    size_t *sizes;
    sizes = (size_t*) malloc(r*sizeof(*sizes));
    for(i=0; i<r; i++) sizes[i] = n;

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
            ap = (unsigned int*) malloc(r * sizeof(int));
        } else {
            as = PROTECT(Rf_allocVector(INTSXP, r));
            Rf_defineVar(Rf_install("a"), as, state);
            UNPROTECT(1);
            ap = (unsigned int*) INTEGER(as);
        }
        for(i=0; i<r; i++) ap[i] = 0;

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
            result = PROTECT(Rf_allocVector(INTSXP, r*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(ltype, r*d));
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
                if (!next_cartesian_product(ap, r, sizes)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                for (i=0; i<r; i++) {
                    result_intp[j + i*d] = ap[i] + 1;
                }
            } else if (ltype == INTSXP) {
                for (i=0; i<r; i++) {
                    result_intp[j + i*d] = labels_intp[ap[i]];
                }
            } else if (ltype == REALSXP) {
                for (i=0; i<r; i++) {
                    result_doublep[j + i*d] = labels_doublep[ap[i]];
                }
            } else if (ltype == STRSXP) {
                for (i=0; i<r; i++) {
                    SET_STRING_ELT(result, j + i*d, STRING_ELT(labels, ap[i]));
                }
            }
        }
        if (status == 0) {
            result = PROTECT(resize_row(result, r, d, j));
            nprotect++;
        }
    } else if (type == 'c') {
        if (labels == R_NilValue) {
            result = PROTECT(Rf_allocVector(INTSXP, r*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(ltype, r*d));
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
                if (!next_cartesian_product(ap, r, sizes)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                for (i=0; i<r; i++) {
                    result_intp[j * n + i] = ap[i] + 1;
                }
            } else if (ltype == INTSXP) {
                for (i=0; i<r; i++) {
                    result_intp[j * n + i] = labels_intp[ap[i]];
                }
            } else if (ltype == REALSXP) {
                for (i=0; i<r; i++) {
                    result_doublep[j * n + i] = labels_doublep[ap[i]];
                }
            } else if (ltype == STRSXP) {
                for (i=0; i<r; i++) {
                    SET_STRING_ELT(result, j * n + i, STRING_ELT(labels, ap[i]));
                }
            }
        }
        if (status == 0) {
            result = PROTECT(resize_col(result, r, d, j));
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
                if (!next_cartesian_product(ap, r, sizes)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                resulti = Rf_allocVector(INTSXP, r);
                result_intp = INTEGER(resulti);
                for (i=0; i<r; i++) {
                    result_intp[i] = ap[i] + 1;
                }
            } else if (ltype == INTSXP) {
                resulti = Rf_allocVector(INTSXP, r);
                result_intp = INTEGER(resulti);
                for (i=0; i<r; i++) {
                    result_intp[i] = labels_intp[ap[i]];
                }
            } else if (ltype == REALSXP) {
                resulti = Rf_allocVector(REALSXP, r);
                result_doublep = REAL(resulti);
                for (i=0; i<r; i++) {
                    result_doublep[i] = labels_doublep[ap[i]];
                }
            } else if (ltype == STRSXP) {
                resulti = Rf_allocVector(STRSXP, r);
                for (i=0; i<r; i++) {
                    SET_STRING_ELT(resulti, i, STRING_ELT(labels, ap[i]));
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

    free(sizes);

    UNPROTECT(nprotect);
    return result;
}
