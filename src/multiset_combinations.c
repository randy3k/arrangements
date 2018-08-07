#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/multiset_combination.h"
#include "utils.h"

double n_multiset_combinations(int* freq, size_t flen, size_t k);

SEXP next_multiset_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout) {
    size_t i, j, h;

    size_t n = as_uint(_n);
    size_t k = as_uint(_k);
    int d;
    double dd;
    if (Rf_asInteger(_d) == -1) {
        dd = n_multiset_combinations(INTEGER(freq), Rf_length(freq), k);
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
        if (dd * k > R_XLEN_T_MAX) Rf_error("too many results");
    }
    d = round(dd);

    SEXP ms;
    unsigned int* mp;
    SEXP as;
    unsigned int* ap;
    int nprotect = 0;

    int status = 0;

    if (state == R_NilValue) {
        ms = R_UnboundValue;
        as = R_UnboundValue;
    } else {
        ms = PROTECT(Rf_findVarInFrame(state, Rf_install("m")));
        as = PROTECT(Rf_findVarInFrame(state, Rf_install("a")));
        nprotect += 2;
    }

    if (ms == R_UnboundValue) {
        if (state == R_NilValue) {
            mp = (unsigned int*) R_alloc(n, sizeof(int));
        } else {
            ms = PROTECT(Rf_allocVector(INTSXP, n));
            Rf_defineVar(Rf_install("m"), ms, state);
            UNPROTECT(1);
            mp = (unsigned int*) INTEGER(ms);
        }
        fp = INTEGER(freq);
        h = 0;
        for (i = 0; i< Rf_length(freq); i++) {
            for (j = 0; j< fp[i]; j++) {
                mp[h++] = i;
            }
        }

    } else {
        mp = (unsigned int*) INTEGER(ms);
    }

    if (as == R_UnboundValue) {
        if (state == R_NilValue) {
            ap = (unsigned int*) R_alloc(k, sizeof(int));
        } else {
            as = PROTECT(Rf_allocVector(INTSXP, k));
            Rf_defineVar(Rf_install("a"), as, state);
            UNPROTECT(1);
            ap = (unsigned int*) INTEGER(as);
        }
        fp = INTEGER(freq);
        h = 0;
        for (i = 0; i< k; i++) {
            ap[i] = mp[i];
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
                if (!next_multiset_combination(mp, ap, n, k)) {
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

    } else if (layout == 'c') {
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
                if (!next_multiset_combination(mp, ap, n, k)) {
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
                if (!next_multiset_combination(mp, ap, n, k)) {
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


double n_multiset_combinations(int* freq, size_t flen, size_t k) {
    int n = 0;
    int i, j, h;
    for (i=0; i<flen; i++) n += freq[i];
    if (k > n) {
        return 0;
    }

    double* p = (double*) malloc((k+1) * sizeof(double));
    for (j=0; j<=k; j++) p[j] = 0;

    double ptemp;

    for (i=0; i<flen; i++) {
        if (i == 0) {
            for (j=0; j<=k && j<=freq[i]; j++) {
                p[j] = 1;
            }
            ptemp = p[k];
        } else if (i < flen - 1){
            for (j=k; j>0; j--) {
                ptemp = 0;
                for(h=0; h<=freq[i] && h<=j; h++) {
                    ptemp += p[j-h];
                }
                p[j] = ptemp;
            }
        } else {
            ptemp = 0;
            for(h=0; h<=freq[i] && h<=k; h++) {
                ptemp += p[k-h];
            }
        }
    }

    free(p);
    return ptemp;
}

SEXP num_multiset_combinations(SEXP freq, SEXP _k) {
    int* fp = INTEGER(freq);
    size_t flen = Rf_length(freq);
    size_t k = as_uint(_k);
    return Rf_ScalarReal(n_multiset_combinations(fp, flen, k));
}

void n_multiset_combinations_bigz(mpz_t z, int* freq, size_t flen, size_t k) {
    int n = 0;
    int i, j, h;
    for (i=0; i<flen; i++) n += freq[i];
    if (k > n) {
        mpz_set_ui(z, 0);
        return;
    }

    mpz_t* p = (mpz_t*) malloc((k+1) * sizeof(mpz_t));
    for (j=0; j<=k; j++) mpz_init(p[j]);

    for (i=0; i<flen; i++) {
        if (i == 0) {
            for (j=0; j<=k && j<=freq[i]; j++) {
                mpz_set_ui(p[j], 1);
            }
            mpz_set(z, p[k]);
        } else if (i < flen - 1){
            for (j=k; j>0; j--) {
                mpz_set_ui(z, 0);
                for(h=0; h<=freq[i] && h<=j; h++) {
                    mpz_add(z, z, p[j-h]);
                }
                mpz_set(p[j], z);
            }
        } else {
            mpz_set_ui(z, 0);
            for(h=0; h<=freq[i] && h<=k; h++) {
                mpz_add(z, z, p[k-h]);
            }
        }
    }
}

SEXP num_multiset_combinations_bigz(SEXP freq, SEXP _k) {
    int* fp = INTEGER(freq);
    size_t flen = Rf_length(freq);
    size_t k = as_uint(_k);
    mpz_t z;
    mpz_init(z);
    n_multiset_combinations_bigz(z, fp, flen, k);
    char* c = mpz_get_str(NULL, 10, z);
    SEXP out = Rf_mkString(c);
    mpz_clear(z);
    free(c);
    return out;
}


void ith_multiset_combination(unsigned int* ar, int* freq, size_t flen, size_t k, unsigned int index) {
    unsigned int i, j;
    unsigned int start = 0;
    unsigned int count, this_count;
    int* subfreq = (int*) malloc(flen * sizeof(int));

    for (i = 0; i < flen; i++) subfreq[i] = freq[i];

    for (i = 0; i < k; i++) {
        count = 0;
        for (j = start; j < flen; j++) {
            if (subfreq[j] == 0) continue;
            subfreq[j]--;
            this_count = count + n_multiset_combinations(subfreq, flen, k - i - 1);
            if (this_count > index) {
                ar[i] = j;
                start = j;
                index -= count;
                break;
            }
            count = this_count;
            subfreq[j] = 0;
        }
    }

    free(subfreq);
}

void ith_multiset_combination_bigz(unsigned int* ar, int* freq, size_t flen, size_t k, mpz_t index) {
    unsigned int i, j;
    unsigned int start = 0;
    mpz_t count;
    mpz_init(count);
    mpz_t this_count;
    mpz_init(this_count);

    int* subfreq = (int*) malloc(flen * sizeof(int));

    for (i = 0; i < flen; i++) subfreq[i] = freq[i];

    for (i = 0; i < k; i++) {
        mpz_set_ui(count, 0);
        for (j = start; j < flen; j++) {
            if (subfreq[j] == 0) continue;
            subfreq[j]--;
            n_multiset_combinations_bigz(this_count, subfreq, flen, k - i - 1);
            mpz_add(this_count, this_count, count);
            if (mpz_cmp(this_count, index) > 0) {
                ar[i] = j;
                start = j;
                mpz_sub(index, index, count);
                break;
            }
            mpz_set(count, this_count);
            subfreq[j] = 0;
        }
    }

    free(subfreq);
    mpz_clear(count);
    mpz_clear(this_count);
}

SEXP get_ith_multiset_combination(SEXP freq, SEXP _k, SEXP _index) {
    unsigned int i;
    int* fp = INTEGER(freq);
    size_t flen = Rf_length(freq);
    int k = as_uint(_k);
    SEXP as = PROTECT(Rf_allocVector(INTSXP, k));
    unsigned int* ar = (unsigned int*) INTEGER(as);

    if (TYPEOF(_index) == STRSXP || n_multiset_combinations(fp, flen, k) > INT_MAX) {
        mpz_t z;
        mpz_init(z);

        if (TYPEOF(_index) == STRSXP) {
            mpz_set_str(z, CHAR(STRING_ELT(_index, 0)), 10);
            mpz_sub_ui(z, z, 1);
        } else {
            mpz_set_ui(z, as_uint(_index) - 1);
        }
        ith_multiset_combination_bigz(ar, fp, flen, k, z);
        mpz_clear(z);
    } else {
        ith_multiset_combination(ar, fp, flen, k, as_uint(_index) - 1);
    }

    for (i = 0; i < k; i++) {
        ar[i]++;
    }
    UNPROTECT(1);
    return as;
}
