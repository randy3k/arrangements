#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/k_permutation.h"
#include "utils.h"

double npermutations_f(int* freq, size_t flen, size_t k);

SEXP next_k_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type) {
    size_t i, j;

    size_t n = as_uint(_n);
    size_t k = as_uint(_k);
    int d;
    double dd;
    if (Rf_asInteger(_d) == -1) {
        dd = fallfact(n, k);
    } else {
        dd = as_uint(_d);
    }

    int labels_type = TYPEOF(labels);
    int* labels_intp;
    double* labels_doublep;

    char type;
    if (_type == R_NilValue) {
        type = 'r';
    } else {
        type = CHAR(Rf_asChar(_type))[0];
        if (type != 'r' && type != 'c' && type != 'l') type = 'r';
    }

    if (dd > INT_MAX) Rf_error("too many results");
    if (type != 'l') {
        if (dd * k > R_XLEN_T_MAX) Rf_error("too many results");
    }
    d = round(dd);

    SEXP as, cycles;
    unsigned int* ap;
    unsigned int* cyclep;
    int nprotect = 0;

    int status = 0;

    if (state == R_NilValue) {
        as = R_UnboundValue;
        cycles = R_UnboundValue;
    } else {
        as = Rf_findVarInFrame(state, Rf_install("a"));
        cycles = Rf_findVarInFrame(state, Rf_install("cycle"));
    }

    if (as == R_UnboundValue) {
        if (state == R_NilValue) {
            ap = (unsigned int*) R_alloc(n, sizeof(int));
            cyclep = (unsigned int*) R_alloc(k, sizeof(int));
        } else {
            as = PROTECT(Rf_allocVector(INTSXP, n));
            Rf_defineVar(Rf_install("a"), as, state);
            UNPROTECT(1);
            ap = (unsigned int*) INTEGER(as);

            cycles = PROTECT(Rf_allocVector(INTSXP, k));
            Rf_defineVar(Rf_install("cycle"), cycles, state);
            UNPROTECT(1);
            cyclep = (unsigned int*) INTEGER(cycles);
        }
        for(i=0; i<n; i++) ap[i] = i;
        for(i=0; i<k; i++) cyclep[i] = n - i;

    } else {
        ap = (unsigned int*) INTEGER(as);
        cyclep = (unsigned int*) INTEGER(cycles);
        status = 1;
    }

    SEXP rdim;
    SEXP result, resulti;
    int* result_intp;
    double* result_doublep;

    if (type == 'r') {
        if (labels == R_NilValue) {
            result = PROTECT(Rf_allocVector(INTSXP, k*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(labels_type, k*d));
            nprotect++;
            if (labels_type == INTSXP) {
                result_intp = INTEGER(result);
                labels_intp = INTEGER(labels);
            } else if (labels_type == REALSXP) {
                result_doublep = REAL(result);
                labels_doublep = REAL(labels);
            }
        }

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_k_permutation(ap, cyclep, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (labels_type == NILSXP) {
                for (i=0; i<k; i++) {
                    result_intp[j + i*d] = ap[i] + 1;
                }
            } else if (labels_type == INTSXP) {
                for (i=0; i<k; i++) {
                    result_intp[j + i*d] = labels_intp[ap[i]];
                }
            } else if (labels_type == REALSXP) {
                for (i=0; i<k; i++) {
                    result_doublep[j + i*d] = labels_doublep[ap[i]];
                }
            } else if (labels_type == STRSXP) {
                for (i=0; i<k; i++) {
                    SET_STRING_ELT(result, j + i*d, STRING_ELT(labels, ap[i]));
                }
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
        if (labels == R_NilValue) {
            result = PROTECT(Rf_allocVector(INTSXP, k*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(labels_type, k*d));
            nprotect++;
            if (labels_type == INTSXP) {
                result_intp = INTEGER(result);
                labels_intp = INTEGER(labels);
            } else if (labels_type == REALSXP) {
                result_doublep = REAL(result);
                labels_doublep = REAL(labels);
            }
        }

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_k_permutation(ap, cyclep, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (labels_type == NILSXP) {
                for (i=0; i<k; i++) {
                    result_intp[j * k + i] = ap[i] + 1;
                }
            } else if (labels_type == INTSXP) {
                for (i=0; i<k; i++) {
                    result_intp[j * k + i] = labels_intp[ap[i]];
                }
            } else if (labels_type == REALSXP) {
                for (i=0; i<k; i++) {
                    result_doublep[j * k + i] = labels_doublep[ap[i]];
                }
            } else if (labels_type == STRSXP) {
                for (i=0; i<k; i++) {
                    SET_STRING_ELT(result, j * k + i, STRING_ELT(labels, ap[i]));
                }
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

    } else {
        // type == "list"
        result = PROTECT(Rf_allocVector(VECSXP, d));
        nprotect++;
        if (labels_type == INTSXP) {
            labels_intp = INTEGER(labels);
        } else if (labels_type == REALSXP) {
            labels_doublep = REAL(labels);
        }

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_k_permutation(ap, cyclep, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (labels_type == NILSXP) {
                resulti = Rf_allocVector(INTSXP, k);
                result_intp = INTEGER(resulti);
                for (i=0; i<k; i++) {
                    result_intp[i] = ap[i] + 1;
                }
            } else if (labels_type == INTSXP) {
                resulti = Rf_allocVector(INTSXP, k);
                result_intp = INTEGER(resulti);
                for (i=0; i<k; i++) {
                    result_intp[i] = labels_intp[ap[i]];
                }
            } else if (labels_type == REALSXP) {
                resulti = Rf_allocVector(REALSXP, k);
                result_doublep = REAL(resulti);
                for (i=0; i<k; i++) {
                    result_doublep[i] = labels_doublep[ap[i]];
                }
            } else if (labels_type == STRSXP) {
                resulti = Rf_allocVector(STRSXP, k);
                for (i=0; i<k; i++) {
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

    UNPROTECT(nprotect);
    return result;
}

SEXP nperm_k(SEXP _n, SEXP _k) {
    size_t n = as_uint(_n);
    size_t k = as_uint(_k);
    return Rf_ScalarReal(fallfact(n, k));
}

void npermutations_k_bigz(mpz_t p, size_t n, size_t k) {
    size_t i;
    if (n < k) {
        mpz_set_ui(p, 0);
        return;
    }
    mpz_set_ui(p, 1);
    for(i=0; i<k; i++) {
        mpz_mul_ui(p, p, n - i);
    }
}

SEXP nperm_k_bigz(SEXP _n, SEXP _k) {
    size_t n = as_uint(_n);
    size_t k = as_uint(_k);
    mpz_t p;
    mpz_init(p);
    npermutations_k_bigz(p, n, k);
    char* c = mpz_get_str(NULL, 10, p);
    SEXP out = Rf_mkString(c);
    free(c);
    mpz_clear(p);
    return out;
}


void ith_permutation_k(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    unsigned int i, j;

    for (i = 0; i < k; i++) {
        j = fallfact(n - 1 - i, k - 1 - i);
        ar[i] = index / j;
        index = index % j;
    }

    for (i = k - 1; i > 0; i--) {
        j = i;
        while (j-- > 0) {
            if (ar[j] <= ar[i]) {
                ar[i]++;
            }
        }
    }
}

void ith_permutation_k_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    unsigned int i, j;

    mpz_t q;
    mpz_init(q);
    mpz_t p;
    mpz_init(p);

    for (i = 0; i < k; i++) {
        npermutations_k_bigz(p, n - 1 - i, k - 1 - i);
        mpz_tdiv_qr(q, index, index, p);
        ar[i] = mpz_get_ui(q);
    }

    for (i = k - 1; i > 0; i--) {
        j = i;
        while (j-- > 0) {
            if (ar[j] <= ar[i]) {
                ar[i]++;
            }
        }
    }

    mpz_clear(q);
    mpz_clear(p);
}

SEXP ith_perm_k(SEXP _n, SEXP _k, SEXP _index) {
    unsigned int i;
    int n = as_uint(_n);
    int k = as_uint(_k);
    SEXP as = PROTECT(Rf_allocVector(INTSXP, k));
    unsigned int* ar = (unsigned int*) INTEGER(as);

    if (fallfact(n, k) > INT_MAX || TYPEOF(_index) == STRSXP) {
        mpz_t z;
        mpz_init(z);

        if (TYPEOF(_index) == STRSXP) {
            mpz_set_str(z, CHAR(STRING_ELT(_index, 0)), 10);
            mpz_sub_ui(z, z, 1);
        } else {
            mpz_set_ui(z, as_uint(_index) - 1);
        }
        ith_permutation_k_bigz(ar, n, k, z);
        mpz_clear(z);
    } else {
        ith_permutation_k(ar, n, k, as_uint(_index) - 1);
    }

    for (i = 0; i < k; i++) {
        ar[i]++;
    }
    UNPROTECT(1);
    return as;
}
