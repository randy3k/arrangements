#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "next.h"
#include "macros.h"
#include "utils.h"


double n_multiset_combinations(int* fp, size_t flen, size_t k) {
    int n = 0;
    int i, j, h;
    for (i=0; i<flen; i++) n += fp[i];
    if (k > n) {
        return 0;
    }

    double* p = (double*) malloc((k+1) * sizeof(double));
    for (j=0; j<=k; j++) p[j] = 0;

    double ptemp;

    for (i=0; i<flen; i++) {
        if (i == 0) {
            for (j=0; j<=k && j<=fp[i]; j++) {
                p[j] = 1;
            }
            ptemp = p[k];
        } else if (i < flen - 1){
            for (j=k; j>0; j--) {
                ptemp = 0;
                for(h=0; h<=fp[i] && h<=j; h++) {
                    ptemp += p[j-h];
                }
                p[j] = ptemp;
            }
        } else {
            ptemp = 0;
            for(h=0; h<=fp[i] && h<=k; h++) {
                ptemp += p[k-h];
            }
        }
    }

    free(p);
    return ptemp;
}


void n_multiset_combinations_bigz(mpz_t z, int* fp, size_t flen, size_t k) {
    int n = 0;
    int i, j, h;
    for (i=0; i<flen; i++) n += fp[i];
    if (k > n) {
        mpz_set_ui(z, 0);
        return;
    }

    mpz_t* p = (mpz_t*) malloc((k+1) * sizeof(mpz_t));
    for (j=0; j<=k; j++) mpz_init(p[j]);

    for (i=0; i<flen; i++) {
        if (i == 0) {
            for (j=0; j<=k && j<=fp[i]; j++) {
                mpz_set_ui(p[j], 1);
            }
            mpz_set(z, p[k]);
        } else if (i < flen - 1){
            for (j=k; j>0; j--) {
                mpz_set_ui(z, 0);
                for(h=0; h<=fp[i] && h<=j; h++) {
                    mpz_add(z, z, p[j-h]);
                }
                mpz_set(p[j], z);
            }
        } else {
            mpz_set_ui(z, 0);
            for(h=0; h<=fp[i] && h<=k; h++) {
                mpz_add(z, z, p[k-h]);
            }
        }
    }
    for (j=0; j<=k; j++) mpz_clear(p[j]);
    free(p);
}


void nth_multiset_combination(unsigned int* ar, int* fp, size_t flen, size_t k, unsigned int index) {
    unsigned int i, j;
    unsigned int start = 0;
    unsigned int count, this_count;
    int* subfreq = (int*) malloc(flen * sizeof(int));

    for (i = 0; i < flen; i++) subfreq[i] = fp[i];

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


void nth_multiset_combination_bigz(unsigned int* ar, int* fp, size_t flen, size_t k, mpz_t index) {
    unsigned int i, j;
    unsigned int start = 0;
    mpz_t count;
    mpz_init(count);
    mpz_t this_count;
    mpz_init(this_count);

    int* subfreq = (int*) malloc(flen * sizeof(int));

    for (i = 0; i < flen; i++) subfreq[i] = fp[i];

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


SEXP next_multiset_combinations(int* fp, size_t flen, int k, SEXP labels, char layout, int d, SEXP _skip, SEXP state) {
    int i, j, h;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    int n = 0;
    for (i=0; i<flen; i++) n += fp[i];

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_multiset_combinations(fp, flen, k);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* mp;
    unsigned int* ap;

    if (!variable_exists(state, "m", INTSXP, n, (void**) &mp)) {
        h = 0;
        for (i = 0; i < flen; i++) {
            for (j = 0; j < fp[i]; j++) {
                mp[h++] = i;
            }
        }
        status = 0;
    }

    if (!variable_exists(state, "a", INTSXP, n, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            for (i = 0; i < n; i++) {
                ap[i] = mp[i];
            }
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_multiset_combinations_bigz(maxz, fp, flen, k);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set(skipz, 0);
                }
                mpz_clear(maxz);
                nth_multiset_combination_bigz(ap, fp, flen, k, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_multiset_combination(ap, fp, flen, k, skip);
            }

            int* subfreq = (int*) malloc(flen * sizeof(int));

            for (i = 0; i < flen; i++) subfreq[i] = fp[i];
            for (i = 0; i < k; i++) subfreq[ap[i]]--;
            h = k;
            for (i = 0; i< flen; i++) {
                for (j = 0; j< subfreq[i]; j++) {
                    ap[h++] = i;
                }
            }
            free(subfreq);

        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_multiset_combination(mp, ap, n, k)) { \
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
    UNPROTECT(nprotect);
    return result;
}


SEXP catch_multiset_combinations(int* fp, size_t flen, int k, SEXP labels, char layout, SEXP _index, SEXP _nsample) {
    int i, j;
    int nprotect = 0;
    int bigz = 0;
    int sampling = _index == R_NilValue;
    SEXP result = R_NilValue;

    double dd;
    if (sampling) {
        dd = as_uint(_nsample);
    } else if (TYPEOF(_index) == RAWSXP && Rf_inherits(_index, "bigz")) {
        dd = *((int* ) RAW(_index));
        bigz = 1;
    } else {
        dd = Rf_length(_index);
    }
    int d = verify_dimension(dd, k, layout);

    double maxd;
    if (!bigz) {
        maxd = n_multiset_combinations(fp, flen, k);
        bigz = maxd > INT_MAX;
    }

    unsigned int* ap;
    ap = (unsigned int*) R_alloc(k, sizeof(int));

    if (bigz) {
        mpz_t* index;
        gmp_randstate_t randstate;
        mpz_t z;
        mpz_t maxz;
        mpz_init(z);
        mpz_init(maxz);
        n_multiset_combinations_bigz(maxz, fp, flen , k);

        if (sampling) {
            GetRNGstate();
            set_gmp_randstate(randstate);
        } else {
            index = (mpz_t*) R_alloc(d, sizeof(mpz_t));
            for (i = 0; i < d; i++) mpz_init(index[i]);
            int status = as_mpz_array(index, d, _index);
            for(i = 0; i < d; i++) {
                if (status < 0 || mpz_sgn(index[i]) <= 0 || mpz_cmp(index[i], maxz) > 0) {
                    for (i = 0; i < d; i++) mpz_clear(index[i]);
                    mpz_clear(maxz);
                    mpz_clear(z);
                    Rf_error("invalid index");
                }
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                mpz_urandomm(z, randstate, maxz); \
            } else { \
                mpz_sub_ui(z, index[j], 1); \
            } \
            nth_multiset_combination_bigz(ap, fp, flen, k, z);

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
        mpz_clear(maxz);
        if (sampling){
            gmp_randclear(randstate);
            PutRNGstate();
        } else {
            for (i = 0; i < d; i++) mpz_clear(index[i]);
        }

    } else {
        int* index;
        if (sampling) {
            GetRNGstate();
        } else {
            index = as_uint_index(_index);
            for (i = 0; i < d; i++) {
                if (index[i] <= 0 || index[i] > maxd) {
                    Rf_error("invalid index");
                }
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                nth_multiset_combination(ap, fp, flen, k, floor(maxd * unif_rand())); \
            } else { \
                nth_multiset_combination(ap, fp, flen, k, index[j] - 1); \
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

    UNPROTECT(nprotect);
    return result;
}
