#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "stdlib.h"
#include "../macros.h"
#include "../math.h"
#include "../utils.h"


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


void n_multiset_permutations_bigz(mpz_t z, int* freq, size_t flen, size_t k) {
    int n = 0;
    int i, j, h;
    for (i=0; i<flen; i++) n += freq[i];
    if (k > n) {
        mpz_set_ui(z, 0);
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


// A Simple, Efficient P(n,k) Algorithm by Alistair Israel
// http://alistairisrael.wordpress.com/2009/09/22/simple-efficient-pnk-algorithm/
// c implementation by 2017 Randy Lai
// http://randycity.github.io

unsigned int next_multiset_permutation(unsigned int *ar, size_t n, size_t k)
{
    unsigned int i;
    unsigned int j;
    unsigned int edge = k-1;

    if(k<n){
        j = k;
        // search for largest j such that nth_j > nth_edge (a is increasing for j>=k)
        while(j<n && ar[edge]>=ar[j]) j++;
    }
    if(k<n && j<n){
        swap(ar, edge, j);
    }else{
        if (k<n){
            reverse(ar+k, n-k);
        }

        // find rightmost ascent to left of edge
        for (i = edge -1; ; i--) {
            if (ar[i] < ar[i+1]) {
                break;
            } else if (i == 0) {
                return 0;
            }
        }

        // find smallest j>=i+1 where nth_j>nth_i (a is decreasing for j>=i+1)
        j = n-1;
        while(j>i && ar[i] >= ar[j]) j--;

        swap(ar, i, j);

        reverse(ar+i+1, n-i-1);
    }

    return 1;
}


void nth_multiset_permutation(unsigned int* ar, int* freq, size_t flen, size_t k, unsigned int index) {
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


void nth_multiset_permutation_bigz(unsigned int* ar, int* freq, size_t flen, size_t k, mpz_t index) {
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


SEXP next_multiset_permutations(int* fp, size_t flen, int k, SEXP labels, char layout, int d, SEXP _skip, SEXP state) {
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
        maxd = n_multiset_permutations(fp, flen, k);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;

    if (!variable_exists(state, (char*)"a", INTSXP, n, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            h = 0;
            for (i = 0; i< flen; i++) {
                for (j = 0; j< fp[i]; j++) {
                    ap[h++] = i;
                }
            }
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_multiset_permutations_bigz(maxz, fp, flen, k);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set_ui(skipz, 0);
                }
                mpz_clear(maxz);
                nth_multiset_permutation_bigz(ap, fp, flen, k, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                nth_multiset_permutation(ap, fp, flen, k, skip);
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
    } else if (labels_type == LGLSXP) {
        RESULT_LGLSXP(k);
    } else {
        Rf_error("label type not supported");
    }

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    UNPROTECT(nprotect);
    return result;
}


SEXP draw_multiset_permutations(int* fp, size_t flen, int k, SEXP labels, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = n_multiset_permutations(fp, flen, k);
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
        n_multiset_permutations_bigz(maxz, fp, flen , k);

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
            nth_multiset_permutation_bigz(ap, fp, flen, k, z);

        int labels_type = TYPEOF(labels);
        if (labels_type == NILSXP) {
            RESULT_NILSXP(k);
        } else if (labels_type == INTSXP) {
            RESULT_INTSXP(k);
        } else if (labels_type == REALSXP) {
            RESULT_REALSXP(k);
        } else if (labels_type == STRSXP) {
            RESULT_STRSXP(k);
        } else if (labels_type == LGLSXP) {
            RESULT_LGLSXP(k);
        } else {
            Rf_error("label type not supported");
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
                nth_multiset_permutation(ap, fp, flen, k, floor(maxd * unif_rand())); \
            } else { \
                nth_multiset_permutation(ap, fp, flen, k, index[j] - 1); \
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
        } else if (labels_type == LGLSXP) {
            RESULT_LGLSXP(k);
        } else {
            Rf_error("label type not supported");
        }

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}
