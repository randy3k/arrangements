#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "arrangements.h"
#include "next/multiset_permutation.h"
#include "utils.h"
#include "macros.h"

double n_multiset_permutations(int* freq, size_t flen, size_t k);

SEXP next_multiset_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout) {
    int i, j, h;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    int n = as_uint(_n);
    int k = as_uint(_k);
    int* fp = INTEGER(freq);
    char layout = check_layout(_layout);

    double dd = Rf_asInteger(_d) == -1 ? n_multiset_permutations(fp, Rf_length(freq), k) : as_uint(_d);
    int d = check_dimension(dd, k, layout);

    unsigned int* ap;

    if (!variable_exist(state, "a", INTSXP, n, (void**) &ap)) {
        if (freq == R_NilValue) {
            for(i=0; i<n; i++) ap[i] = i;
        } else {
            h = 0;
            for (i = 0; i< Rf_length(freq); i++) {
                for (j = 0; j< fp[i]; j++) {
                    ap[h++] = i;
                }
            }
        }
        status = 0;
    }

    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_multiset_permutation(ap, n, k)) { \
            status = 0; \
            break; \
        }

    int labels_type = TYPEOF(labels);
    if (labels_type == NILSXP) {
        RESULT_NILSXP(k);
    } else if (labels_type == INTSXP) {
        RESULT_INTSXP(k);
    } else if (labels_type == REALSXP) {
        RESULT_REALSXP(k);
    } else if (labels_type == STRSXP) {
        RESULT_STRSXP(k);
    }

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    check_factor(result, labels);
    UNPROTECT(nprotect);
    return result;
}


double n_multiset_permutations(int* freq, size_t flen, size_t k) {
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

SEXP num_multiset_permutations(SEXP freq, SEXP _k) {
    int* fp = INTEGER(freq);
    size_t flen = Rf_length(freq);
    size_t k = as_uint(_k);
    return Rf_ScalarReal(n_multiset_permutations(fp, flen, k));
}

void n_multiset_permutations_bigz(mpz_t z, int* freq, size_t flen, size_t k) {
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

SEXP num_multiset_permutations_bigz(SEXP freq, SEXP _k) {
    int* fp = INTEGER(freq);
    size_t flen = Rf_length(freq);
    size_t k = as_uint(_k);
    mpz_t z;
    mpz_init(z);
    n_multiset_permutations_bigz(z, fp, flen, k);
    char* c = mpz_get_str(NULL, 10, z);
    SEXP out = Rf_mkString(c);
    mpz_clear(z);
    free(c);
    return out;
}


void ith_multiset_permutation(unsigned int* ar, int* freq, size_t flen, size_t k, unsigned int index) {
    unsigned int i, j;
    unsigned int count, this_count;
    int* subfreq = (int*) malloc(flen * sizeof(int));

    for (i = 0; i < flen; i++) subfreq[i] = freq[i];

    for (i = 0; i < k; i++) {
        count = 0;
        for (j = 0; j < flen; j++) {
            if (subfreq[j] == 0) continue;
            subfreq[j]--;
            this_count = count + n_multiset_permutations(subfreq, flen, k - i - 1);
            if (this_count > index) {
                ar[i] = j;
                index -= count;
                break;
            }
            count = this_count;
            subfreq[j]++;
        }
    }

    free(subfreq);
}

void ith_multiset_permutation_bigz(unsigned int* ar, int* freq, size_t flen, size_t k, mpz_t index) {
    unsigned int i, j;
    mpz_t count;
    mpz_init(count);
    mpz_t this_count;
    mpz_init(this_count);

    int* subfreq = (int*) malloc(flen * sizeof(int));

    for (i = 0; i < flen; i++) subfreq[i] = freq[i];

    for (i = 0; i < k; i++) {
        mpz_set_ui(count, 0);
        for (j = 0; j < flen; j++) {
            if (subfreq[j] == 0) continue;
            subfreq[j]--;
            n_multiset_permutations_bigz(this_count, subfreq, flen, k - i - 1);
            mpz_add(this_count, this_count, count);
            if (mpz_cmp(this_count, index) > 0) {
                ar[i] = j;
                mpz_sub(index, index, count);
                break;
            }
            mpz_set(count, this_count);
            subfreq[j]++;
        }
    }

    free(subfreq);
    mpz_clear(count);
    mpz_clear(this_count);
}

SEXP get_multiset_permutation(SEXP freq, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample) {
    int i, j;
    int nprotect = 0;
    SEXP result = R_NilValue;

    int* fp = INTEGER(freq);
    size_t flen = Rf_length(freq);
    int k = as_uint(_k);
    char layout = check_layout(_layout);

    double dd = _index == R_NilValue ? as_uint(_nsample) : Rf_length(_index);
    int d = check_dimension(dd, k, layout);

    unsigned int* ap;
    ap = (unsigned int*) R_alloc(k, sizeof(int));

    double max = n_multiset_permutations(fp, flen, k);
    int sampling = _index == R_NilValue;
    int bigz = TYPEOF(_index) == STRSXP || max > INT_MAX;

    if (bigz) {
        gmp_randstate_t randstate;
        mpz_t z;
        mpz_t maxz;
        mpz_init(z);

        if (sampling) {
            GetRNGstate();
            set_gmp_state(randstate);
            mpz_init(maxz);
            n_multiset_permutations_bigz(maxz, fp, flen, k);
        } else {
            if (TYPEOF(_index) != STRSXP) {
                _index = Rf_coerceVector(_index, STRSXP);
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                mpz_urandomm(z, randstate, maxz); \
            } else { \
                mpz_set_str(z, CHAR(STRING_ELT(_index, j)), 10); \
                mpz_sub_ui(z, z, 1); \
            } \
            ith_multiset_permutation_bigz(ap, fp, flen, k, z);

        int labels_type = TYPEOF(labels);
        if (labels_type == NILSXP) {
            RESULT_NILSXP(k);
        } else if (labels_type == INTSXP) {
            RESULT_INTSXP(k);
        } else if (labels_type == REALSXP) {
            RESULT_REALSXP(k);
        } else if (labels_type == STRSXP) {
            RESULT_STRSXP(k);
        }

        mpz_clear(z);
        if (sampling){
            mpz_clear(maxz);
            gmp_randclear(randstate);
            PutRNGstate();
        }

    } else {
        int* index;
        if (sampling) {
            GetRNGstate();
        } else {
            index = INTEGER(_index);
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                ith_multiset_permutation(ap, fp, flen, k, floor(max * unif_rand())); \
            } else { \
                ith_multiset_permutation(ap, fp, flen, k, index[j] - 1); \
            }

        int labels_type = TYPEOF(labels);
        if (labels_type == NILSXP) {
            RESULT_NILSXP(k);
        } else if (labels_type == INTSXP) {
            RESULT_INTSXP(k);
        } else if (labels_type == REALSXP) {
            RESULT_REALSXP(k);
        } else if (labels_type == STRSXP) {
            RESULT_STRSXP(k);
        }

        if (sampling){
            PutRNGstate();
        }
    }

    check_factor(result, labels);
    UNPROTECT(nprotect);
    return result;
}
