#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/multiset_permutation.h"
#include "utils.h"

double npermutations_f(int* freq, size_t flen, size_t k);

SEXP next_multiset_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _type) {
    size_t i, j, h;

    size_t n = as_uint(_n);
    size_t k = as_uint(_k);
    int d;
    double dd;
    if (Rf_asInteger(_d) == -1) {
        if (freq == R_NilValue) {
            dd = fallfact(n, k);
        } else {
            dd = npermutations_f(INTEGER(freq), Rf_length(freq), k);
        }
    } else {
        dd = as_uint(_d);
    }

    int ltype = TYPEOF(labels);
    int* labels_intp;
    double* labels_doublep;

    int* fp;

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
        if (dd * k > INT_MAX) Rf_error("too many results");
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

    if (type == 'r') {
        if (labels == R_NilValue) {
            result = PROTECT(Rf_allocVector(INTSXP, k*d));
            nprotect++;
            result_intp = INTEGER(result);
        } else {
            result = PROTECT(Rf_allocVector(ltype, k*d));
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
                if (!next_multiset_permutation(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                for (i=0; i<k; i++) {
                    result_intp[j + i*d] = ap[i] + 1;
                }
            } else if (ltype == INTSXP) {
                for (i=0; i<k; i++) {
                    result_intp[j + i*d] = labels_intp[ap[i]];
                }
            } else if (ltype == REALSXP) {
                for (i=0; i<k; i++) {
                    result_doublep[j + i*d] = labels_doublep[ap[i]];
                }
            } else if (ltype == STRSXP) {
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
            result = PROTECT(Rf_allocVector(ltype, k*d));
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
                if (!next_multiset_permutation(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                for (i=0; i<k; i++) {
                    result_intp[j * k + i] = ap[i] + 1;
                }
            } else if (ltype == INTSXP) {
                for (i=0; i<k; i++) {
                    result_intp[j * k + i] = labels_intp[ap[i]];
                }
            } else if (ltype == REALSXP) {
                for (i=0; i<k; i++) {
                    result_doublep[j * k + i] = labels_doublep[ap[i]];
                }
            } else if (ltype == STRSXP) {
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
        if (ltype == INTSXP) {
            labels_intp = INTEGER(labels);
        } else if (ltype == REALSXP) {
            labels_doublep = REAL(labels);
        }

        for (j=0; j<d; j++) {
            if (status) {
                if (!next_multiset_permutation(ap, n, k)) {
                    status = 0;
                    break;
                }
            } else {
                status = 1;
            }
            if (ltype == NILSXP) {
                resulti = Rf_allocVector(INTSXP, k);
                result_intp = INTEGER(resulti);
                for (i=0; i<k; i++) {
                    result_intp[i] = ap[i] + 1;
                }
            } else if (ltype == INTSXP) {
                resulti = Rf_allocVector(INTSXP, k);
                result_intp = INTEGER(resulti);
                for (i=0; i<k; i++) {
                    result_intp[i] = labels_intp[ap[i]];
                }
            } else if (ltype == REALSXP) {
                resulti = Rf_allocVector(REALSXP, k);
                result_doublep = REAL(resulti);
                for (i=0; i<k; i++) {
                    result_doublep[i] = labels_doublep[ap[i]];
                }
            } else if (ltype == STRSXP) {
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


double npermutations_f(int* freq, size_t flen, size_t k) {
    int n = 0;
    int i, j, h;
    for (i=0; i<flen; i++) n += freq[i];
    if (k > n) {
        return 0;
    }

    int maxf;
    maxf = 0;
    for (i=0; i<flen; i++) {
        if (freq[i] > maxf) maxf = freq[i];
    }

    double rfact = 1;
    for (j=2; j<=k; j++) {
        rfact = rfact*j;
    }
    size_t factlen = (k < maxf ? k : maxf) + 1;
    double* fact = (double*) malloc(factlen * sizeof(double));
    fact[0] = 1;
    for (j=1; j< factlen; j++) fact[j] = j * fact[j-1];

    double* p = (double*) malloc((k+1) * sizeof(double));
    for (j=0; j<=k; j++) p[j] = 0;

    double ptemp;

    for (i=0; i<flen; i++) {
        if (i == 0) {
            for (j=0; j<=k && j<=freq[i]; j++) {
                p[j] = rfact / fact[j];
            }
            ptemp = p[k];
        } else if (i < flen - 1){
            for (j=k; j>0; j--) {
                ptemp = 0;
                for(h=0; h<=freq[i] && h<=j; h++) {
                    ptemp += p[j-h] / fact[h];
                }
                p[j] = ptemp;
            }
        } else {
            ptemp = 0;
            for(h=0; h<=freq[i] && h<=k; h++) {
                ptemp += p[k-h] / fact[h];
            }
        }
    }

    free(fact);
    free(p);
    return ptemp;
}

SEXP nperm_f(SEXP freq, SEXP _k) {
    int* fp = INTEGER(freq);
    size_t flen = Rf_length(freq);
    size_t k = as_uint(_k);
    return Rf_ScalarReal(npermutations_f(fp, flen, k));
}

void npermutations_f_bigz(mpz_t z, int* freq, size_t flen, size_t k) {
    int n = 0;
    int i, j, h;
    for (i=0; i<flen; i++) n += freq[i];
    if (k > n) {
        mpz_set(z, 0);
        return;
    }

    int maxf;
    maxf = 0;
    for (i=0; i<flen; i++) {
        if (freq[i] > maxf) maxf = freq[i];
    }

    mpz_t rfact;
    mpz_init(rfact);
    mpz_fac_ui(rfact, k);

    size_t factlen = (k < maxf ? k : maxf) + 1;
    mpz_t* fact = (mpz_t*) malloc(factlen * sizeof(mpz_t));
    for (j=0; j< factlen; j++) mpz_init(fact[j]);

    mpz_set_ui(fact[0], 1);
    for (j=1; j< factlen; j++) mpz_mul_ui(fact[j], fact[j-1], j);

    mpz_t* p = (mpz_t*) malloc((k+1) * sizeof(mpz_t));
    for (j=0; j<=k; j++) mpz_init(p[j]);

    mpz_t ptemp;
    mpz_init(ptemp);

    for (i=0; i<flen; i++) {
        if (i == 0) {
            for (j=0; j<=k && j<=freq[i]; j++) {
                mpz_cdiv_q(p[j], rfact, fact[j]);
            }
            mpz_set(z, p[k]);
        } else if (i < flen - 1){
            for (j=k; j>0; j--) {
                mpz_set_ui(z, 0);
                for(h=0; h<=freq[i] && h<=j; h++) {
                    mpz_cdiv_q(ptemp, p[j-h], fact[h]);
                    mpz_add(z, z, ptemp);
                }
                mpz_set(p[j], z);
            }
        } else {
            mpz_set_ui(z, 0);
            for(h=0; h<=freq[i] && h<=k; h++) {
                mpz_cdiv_q(ptemp, p[k-h], fact[h]);
                mpz_add(z, z, ptemp);
            }
        }
    }

    for (j=0; j< factlen; j++) mpz_clear(fact[j]);
    for (j=0; j<=k; j++) mpz_clear(p[j]);
    free(fact);
    free(p);
    mpz_clear(rfact);
    mpz_clear(ptemp);
}

SEXP nperm_f_bigz(SEXP freq, SEXP _k) {
    int* fp = INTEGER(freq);
    size_t flen = Rf_length(freq);
    size_t k = as_uint(_k);
    mpz_t z;
    mpz_init(z);
    npermutations_f_bigz(z, fp, flen, k);
    char* c = mpz_get_str(NULL, 10, z);
    SEXP out = Rf_mkString(c);
    mpz_clear(z);
    free(c);
    return out;
}
