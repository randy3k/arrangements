#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/permutation.h"
#include "utils.h"

SEXP next_permutations(SEXP _n, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout) {
    size_t i, j, h;

    size_t n = as_uint(_n);
    int d;
    double dd;
    if (Rf_asInteger(_d) == -1) {
        if (freq == R_NilValue) {
            dd = fact(n);
        } else {
            dd = multichoose(INTEGER(freq), Rf_length(freq));
        }
    } else {
        dd = as_uint(_d);
    }

    int labels_type = TYPEOF(labels);
    int* labels_intp;
    double* labels_doublep;

    int* fp;

    char layout;
    if (_layout == R_NilValue) {
        layout = 'r';
    } else {
        layout = CHAR(Rf_asChar(_layout))[0];
        if (layout != 'r' && layout != 'c' && layout != 'l') layout = 'r';
    }

    if (dd > INT_MAX) Rf_error("too many results");
    if (layout != 'l') {
        if (dd * n > R_XLEN_T_MAX) Rf_error("too many results");
    }
    d = round(dd);

    SEXP as;
    unsigned int* ap;
    int nprotect = 0;

    int status = 0;

    if (state == R_NilValue) {
        as = R_UnboundValue;
    } else {
        as = Rf_findVarInFrame(state, Rf_install("a"));
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
        if (freq == R_NilValue) {
            for(i=0; i<n; i++) ap[i] = i;
        } else {
            fp = INTEGER(freq);
            h = 0;
            for (i = 0; i< Rf_length(freq); i++) {
                for (j = 0; j< fp[i]; j++) {
                    ap[h++] = i;
                }
            }
        }

    } else {
        ap = (unsigned int*) INTEGER(as);
        status = 1;
    }

    SEXP rdim;
    SEXP result, resulti;
    int* result_intp;
    double* result_doublep;

    if (layout == 'r') {
        if (labels == R_NilValue) {
            result = PROTECT(Rf_allocVector(INTSXP, n*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(labels_type, n*d));
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
                if (!next_permutation(ap, n)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (labels_type == NILSXP) {
                for (i=0; i<n; i++) {
                    result_intp[j + i*d] = ap[i] + 1;
                }
            } else if (labels_type == INTSXP) {
                for (i=0; i<n; i++) {
                    result_intp[j + i*d] = labels_intp[ap[i]];
                }
            } else if (labels_type == REALSXP) {
                for (i=0; i<n; i++) {
                    result_doublep[j + i*d] = labels_doublep[ap[i]];
                }
            } else if (labels_type == STRSXP) {
                for (i=0; i<n; i++) {
                    SET_STRING_ELT(result, j + i*d, STRING_ELT(labels, ap[i]));
                }
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

    } else if (layout == 'c') {
        if (labels == R_NilValue) {
            result = PROTECT(Rf_allocVector(INTSXP, n*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(labels_type, n*d));
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
                if (!next_permutation(ap, n)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (labels_type == NILSXP) {
                for (i=0; i<n; i++) {
                    result_intp[j * n + i] = ap[i] + 1;
                }
            } else if (labels_type == INTSXP) {
                for (i=0; i<n; i++) {
                    result_intp[j * n + i] = labels_intp[ap[i]];
                }
            } else if (labels_type == REALSXP) {
                for (i=0; i<n; i++) {
                    result_doublep[j * n + i] = labels_doublep[ap[i]];
                }
            } else if (labels_type == STRSXP) {
                for (i=0; i<n; i++) {
                    SET_STRING_ELT(result, j * n + i, STRING_ELT(labels, ap[i]));
                }
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

    } else {
        // layout == 'l'
        result = PROTECT(Rf_allocVector(VECSXP, d));
        nprotect++;
        if (labels_type == INTSXP) {
            labels_intp = INTEGER(labels);
        } else if (labels_type == REALSXP) {
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
            if (labels_type == NILSXP) {
                resulti = Rf_allocVector(INTSXP, n);
                result_intp = INTEGER(resulti);
                for (i=0; i<n; i++) {
                    result_intp[i] = ap[i] + 1;
                }
            } else if (labels_type == INTSXP) {
                resulti = Rf_allocVector(INTSXP, n);
                result_intp = INTEGER(resulti);
                for (i=0; i<n; i++) {
                    result_intp[i] = labels_intp[ap[i]];
                }
            } else if (labels_type == REALSXP) {
                resulti = Rf_allocVector(REALSXP, n);
                result_doublep = REAL(resulti);
                for (i=0; i<n; i++) {
                    result_doublep[i] = labels_doublep[ap[i]];
                }
            } else if (labels_type == STRSXP) {
                resulti = Rf_allocVector(STRSXP, n);
                for (i=0; i<n; i++) {
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

SEXP num_multiset_n_permutations(SEXP freq) {
    int* fp = INTEGER(as_uint_array(freq));
    size_t flen = Rf_length(freq);
    return Rf_ScalarReal(multichoose(fp, flen));
}

void n_multiset_n_permutations_bigz(mpz_t z, int* freq, size_t flen) {
    mpz_set_ui(z, 1);
    size_t i, j, h;
    h = 0;
    for (i=0; i<flen; i++) {
        for (j=1; j<=freq[i]; j++) {
            h++;
            mpz_mul_ui(z, z, h);
            mpz_cdiv_q_ui(z, z, j);
        }
    }
}

SEXP num_multiset_n_permutations_bigz(SEXP freq) {
    int* fp = INTEGER(as_uint_array(freq));
    size_t flen = Rf_length(freq);
    mpz_t z;
    mpz_init(z);
    n_multiset_n_permutations_bigz(z, fp, flen);
    char* c = mpz_get_str(NULL, 10, z);
    SEXP out = Rf_mkString(c);
    mpz_clear(z);
    free(c);
    return out;
}

void ith_permutation(unsigned int* ar, unsigned int n, unsigned int index) {
    unsigned int i, j;
    unsigned int* fact = (unsigned int*) malloc(n * sizeof(unsigned int));

    fact[0] = 1;
    for (i = 1; i < n; i++) {
        fact[i] = fact[i-1] * i;
    }
    for (i = 0; i < n; i++) {
        ar[i] = index / fact[n - 1 - i];
        index = index % fact[n - 1 - i];
    }

    for (i = n - 1; i > 0; i--) {
        j = i;
        while (j-- > 0) {
            if (ar[j] <= ar[i]) {
                ar[i]++;
            }
        }
    }

    free(fact);
}

void ith_permutation_bigz(unsigned int* ar, unsigned int n, mpz_t index) {
    unsigned int i, j;

    mpz_t q;
    mpz_init(q);

    mpz_t* fact = (mpz_t*) malloc(n * sizeof(mpz_t));
    for (i=0; i< n; i++) mpz_init(fact[i]);

    mpz_set_ui(fact[0], 1);
    for (i=1; i< n; i++) mpz_mul_ui(fact[i], fact[i-1], i);

    for (i = 0; i < n; i++) {
        mpz_tdiv_qr(q, index, index, fact[n - 1 - i]);
        ar[i] = mpz_get_ui(q);
    }

    for (i = n - 1; i > 0; i--) {
        j = i;
        while (j-- > 0) {
            if (ar[j] <= ar[i]) {
                ar[i]++;
            }
        }
    }

    mpz_clear(q);
    for (i=0; i< n; i++) mpz_clear(fact[i]);
    free(fact);
}

SEXP get_ith_permutation(SEXP _n, SEXP _index) {
    unsigned int i;
    int n = as_uint(_n);
    SEXP as = PROTECT(Rf_allocVector(INTSXP, n));
    unsigned int* ar = (unsigned int*) INTEGER(as);

    if (TYPEOF(_index) == STRSXP || n > 12) {
        mpz_t z;
        mpz_init(z);

        if (TYPEOF(_index) == STRSXP) {
            mpz_set_str(z, CHAR(STRING_ELT(_index, 0)), 10);
            mpz_sub_ui(z, z, 1);
        } else {
            mpz_set_ui(z, as_uint(_index) - 1);
        }
        ith_permutation_bigz(ar, n, z);
        mpz_clear(z);
    } else {
        ith_permutation(ar, n, as_uint(_index) - 1);
    }

    for (i = 0; i < n; i++) {
        ar[i]++;
    }
    UNPROTECT(1);
    return as;
}
