#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>
#include "../combinatorics.h"
#include "../utils.h"
#include "../macros.h"


double n_k_partitions(int n, int k) {
    if (n < k) {
        return 0;
    } else if (n == 0 && k == 0) {
        return 1;
    } else if (k == 0) {
        return 0;
    }
    int n1 = n-k+1;
    double* p = (double*) malloc(n1*k * sizeof(double));
    int i, j, h;

    for (j=0; j<k; j++) {
        p[j] = 1;
    }
    for (i=1; i<n1; i++) {
        p[i*k] = 1;
        for (j=1; j<k; j++) {
            h = i*k + j;
            if (i > j) {
                p[h] =  p[h - 1] + p[h - (j + 1)*k];
            } else {
                p[h] =  p[h - 1];
            }
        }
    }
    double out = p[n1*k - 1];
    free(p);
    return out;
}


void n_k_partitions_bigz(mpz_t z, int n, int k) {
    if (n < k) {
        mpz_set_ui(z, 0);
        return;
    } else if (n == 0 && k == 0) {
        mpz_set_ui(z, 1);
        return;
    } else if (k == 0) {
        mpz_set_ui(z, 0);
        return;
    }

    int n1 = n-k+1;
    int i, j, h;

    mpz_t* p = (mpz_t*) malloc(n1*k * sizeof(mpz_t));
    for (i=0; i<n1*k; i++) mpz_init(p[i]);

    for (j=0; j<k; j++) {
        mpz_set_ui(p[j], 1);
    }
    for (i=1; i<n1; i++) {
        mpz_set_ui(p[i*k], 1);
        for (j=1; j<k; j++) {
            h = i*k + j;
            if (i > j) {
                mpz_add(p[h], p[h - 1], p[h - (j + 1)*k]);
            } else {
                mpz_set(p[h], p[h - 1]);
            }
        }
    }
    mpz_set(z, p[n1*k - 1]);
    for (i=0; i<n1*k; i++) mpz_clear(p[i]);
    free(p);
}

double nkm(int n, int k, int m) {
    // number of partitions of n into at most k parts of sizes <= m
    // note that number of partitions of n into exactly k parts
    // is p(n, k, m) - p(n, k-1, m) = p(n-k, k, m-1)

    if (n > m*k) {
        return 0;
    } else if (n == 0) {
        return 1;
    } else if (k == 0) {
        return 0;
    }

    int i, j, h;
    double* p = (double*) malloc((n + 1) * sizeof(double));
    for (j = 1; j <= n; j++) {
        p[j] = 0;
    }
    p[0] = 1;
    for (i = 1; i <= m; i++) {
        for (j = n; j >= k + i; j--) {
            p[j] -= p[j - k - i];
        }
        for (j = n; j >= 0; j--) {
            for (h = i; h <= j; h += i) {
                p[j] += p[j - h];
            }
        }
    }
    double pn = p[n];
    free(p);
    return pn;
}

double n_k_m_partitions(int n, int k, int m) {
    return nkm(n-k, k, m-1);
}


void nkm_bigz(mpz_t z, int n, int k, int m) {
    // number of partitions of n into at most k parts of sizes <= m
    // note that number of partitions of n into exactly k parts
    // is p(n, k, m) - p(n, k-1, m) = p(n-k, k, m-1)
    if (n > m*k) {
        mpz_set_ui(z, 0);
        return;
    } else if (n == 0) {
        mpz_set_ui(z, 1);
        return;
    } else if (k == 0) {
        mpz_set_ui(z, 0);
        return;
    }

    int i, j, h;
    mpz_t* p = (mpz_t*) malloc((n+1) * sizeof(mpz_t));
    for (j = 0; j <= n; j++) mpz_init(p[j]);
    for (j = 1; j <= n; j++) {
        mpz_set_ui(p[j], 0);
    }
    mpz_set_ui(p[0], 1);
    for (i = 1; i <= m; i++) {
        for (j = n; j >= k + i; j--) {
            mpz_sub(p[j], p[j], p[j - k - i]);
        }
        for (j = n; j >= 0; j--) {
            for (h = i; h <= j; h += i) {
                mpz_add(p[j], p[j], p[j - h]);
            }
        }
    }
    mpz_set(z, p[n]);
    for (j = 0; j <= n; j++) mpz_clear(p[j]);
    free(p);
}

void n_k_m_partitions_bigz(mpz_t z, int n, int k, int m) {
    nkm_bigz(z, n-k, k, m-1);
}


void identify_asc_k_partition(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    unsigned int i, j;
    unsigned int start = 1;
    unsigned int count, this_count;

    for (i = 0; i < k; i++) {
        count = 0;
        for (j = start; j <= n; j++) {
            this_count = count + n_k_partitions(n - j - (j - 1) * (k - i - 1), k - i - 1);
            if (this_count > index) {
                ar[i] = j;
                start = j;
                n -= j;
                index -= count;
                break;
            }
            count = this_count;
        }
    }
}


void identify_asc_k_partition_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    unsigned int i, j;
    unsigned int start = 1;
    mpz_t count, this_count;
    mpz_init(count);
    mpz_init(this_count);

    for (i = 0; i < k; i++) {
        mpz_set_ui(count, 0);
        for (j = start; j <= n; j++) {
            n_k_partitions_bigz(this_count, n - j - (j - 1) * (k - i - 1), k - i - 1);
            mpz_add(this_count, this_count, count);
            if (mpz_cmp(this_count, index) > 0) {
                ar[i] = j;
                start = j;
                n -= j;
                mpz_sub(index, index, count);
                break;
            }
            mpz_set(count, this_count);
        }
    }

    mpz_clear(count);
    mpz_clear(this_count);
}


SEXP next_asc_k_partitions(int n, int k, char layout, int d, SEXP _skip, SEXP state) {
    int i, j;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_k_partitions(n, k);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;

    if (!variable_exists(state, "a", INTSXP, k, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            for(i=0; i<k-1; i++) ap[i] = 1;
            ap[k-1] = n - k + 1;
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_k_partitions_bigz(maxz, n, k);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set(skipz, 0);
                }
                mpz_clear(maxz);
                identify_asc_k_partition_bigz(ap, n, k, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                identify_asc_k_partition(ap, n, k, skip);
            }
        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_asc_k_partition(ap, n, k)) { \
            status = 0; \
            break; \
        }

    RESULT_K_PART();

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    UNPROTECT(nprotect);
    return result;
}


SEXP obtain_asc_k_partitions(int n, int k, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = n_k_partitions(n, k);
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
        n_k_partitions_bigz(maxz, n, k);

        if (sampling) {
            GetRNGstate();
            set_gmp_randstate(randstate);
        } else {
            index = (mpz_t*) R_alloc(d, sizeof(mpz_t));
            for (i = 0; i < d; i++) mpz_init(index[i]);
            int status = as_mpz_array(index, d, _index);
            for(i = 0; i < d; i++) {
                if (status < 0 || mpz_sgn(index[i]) <= 0) {
                    for (i = 0; i < d; i++) mpz_clear(index[i]);
                    mpz_clear(maxz);
                    mpz_clear(z);
                    Rf_error("expect integer");
                } else if (mpz_cmp(index[i], maxz) > 0) {
                    mpz_set(index[i], maxz);
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
            identify_asc_k_partition_bigz(ap, n, k, z);

        RESULT_K_PART();

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
            index = as_uint_array(_index);
            for (i = 0; i < d; i++) {
                if (index[i] <= 0) {
                    Rf_error("expect integer");
                } else if (index[i] > maxd) {
                    index[i] = maxd;
                }
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                identify_asc_k_partition(ap, n, k, floor(maxd * unif_rand())); \
            } else { \
                identify_asc_k_partition(ap, n, k, index[j] - 1); \
            }

        RESULT_K_PART();

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}


void identify_desc_k_partition(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index) {
    unsigned int i, j;
    unsigned int start = n - k + 1;
    unsigned int count, this_count;

    for (i = 0; i < k; i++) {
        count = 0;
        for (j = start; j >= 1; j--) {
            this_count = count + n_k_m_partitions(n - j, k - i - 1, j);
            if (this_count > index) {
                ar[i] = j;
                n -= j;
                start = n - (k - i - 2);
                start = start > j? j : start;
                index -= count;
                break;
            }
            count = this_count;
        }
    }
}


void identify_desc_k_partition_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index) {
    unsigned int i, j;
    unsigned int start = n - k + 1;
    mpz_t count, this_count;
    mpz_init(count);
    mpz_init(this_count);

    for (i = 0; i < k; i++) {
        mpz_set_ui(count, 0);
        for (j = start; j >= 1; j--) {
            n_k_m_partitions_bigz(this_count, n - j, k - i - 1, j);
            mpz_add(this_count, this_count, count);
            if (mpz_cmp(this_count, index) > 0) {
                ar[i] = j;
                n -= j;
                start = n - (k - i - 2);
                start = start > j? j : start;
                mpz_sub(index, index, count);
                break;
            }
            mpz_set(count, this_count);
        }
    }

    mpz_clear(count);
    mpz_clear(this_count);
}


SEXP next_desc_k_partitions(int n, int k, char layout, int d, SEXP _skip, SEXP state) {
    int i, j;
    int nprotect = 0;
    int status = 1;
    SEXP result;

    double dd;
    double maxd;
    int bigz = TYPEOF(_skip) == RAWSXP && Rf_inherits(_skip, "bigz");
    if (d == -1 || !Rf_isNull(_skip)) {
        maxd = n_k_partitions(n, k);
        bigz = bigz || maxd >= INT_MAX;
    }
    dd = d == -1 ? maxd : d;
    d = verify_dimension(dd, n, layout);

    unsigned int* ap;

    if (!variable_exists(state, "a", INTSXP, k, (void**) &ap)) {
        mpz_t maxz;
        int skip;
        mpz_t skipz;
        if (Rf_isNull(_skip)) {
            for(i=1; i<k; i++) ap[i] = 1;
            ap[0] = n - k + 1;
        } else {
            if (bigz) {
                mpz_init(maxz);
                mpz_init(skipz);
                n_k_partitions_bigz(maxz, n, k);
                if (as_mpz_array(&skipz, 1, _skip) < 0 || mpz_sgn(skipz) < 0) {
                    mpz_clear(skipz);
                    mpz_clear(maxz);
                    Rf_error("expect integer");
                } else if (mpz_cmp(skipz, maxz) >= 0) {
                    mpz_set(skipz, 0);
                }
                mpz_clear(maxz);
                identify_desc_k_partition_bigz(ap, n, k, skipz);
                mpz_clear(skipz);
            } else {
                skip = as_uint(_skip);
                if (skip >= (int) maxd) {
                    skip = 0;
                }
                identify_desc_k_partition(ap, n, k, skip);
            }
        }
        status = 0;
    }

    #undef NEXT
    #define NEXT() \
        if (status == 0) { \
            status = 1; \
        } else if (!next_desc_k_partition(ap, n, k)) { \
            status = 0; \
            break; \
        }

    RESULT_K_PART();

    if (status == 0) {
        result = PROTECT(resize_layout(result, j, layout));
        nprotect++;
    }
    UNPROTECT(nprotect);
    return result;
}


SEXP obtain_desc_k_partitions(int n, int k, char layout, SEXP _index, SEXP _nsample) {
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
        maxd = n_k_partitions(n, k);
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
        n_k_partitions_bigz(maxz, n, k);

        if (sampling) {
            GetRNGstate();
            set_gmp_randstate(randstate);
        } else {
            index = (mpz_t*) R_alloc(d, sizeof(mpz_t));
            for (i = 0; i < d; i++) mpz_init(index[i]);
            int status = as_mpz_array(index, d, _index);
            for(i = 0; i < d; i++) {
                if (status < 0 || mpz_sgn(index[i]) <= 0) {
                    for (i = 0; i < d; i++) mpz_clear(index[i]);
                    mpz_clear(maxz);
                    mpz_clear(z);
                    Rf_error("expect integer");
                } else if (mpz_cmp(index[i], maxz) > 0) {
                    mpz_set(index[i], maxz);
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
            identify_desc_k_partition_bigz(ap, n, k, z);

        RESULT_K_PART();

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
            index = as_uint_array(_index);
            for (i = 0; i < d; i++) {
                if (index[i] <= 0) {
                    Rf_error("expect integer");
                } else if (index[i] > maxd) {
                    index[i] = maxd;
                }
            }
        }

        #undef NEXT
        #define NEXT() \
            if (sampling) { \
                identify_desc_k_partition(ap, n, k, floor(maxd * unif_rand())); \
            } else { \
                identify_desc_k_partition(ap, n, k, index[j] - 1); \
            }

        RESULT_K_PART();

        if (sampling){
            PutRNGstate();
        }
    }

    UNPROTECT(nprotect);
    return result;
}
