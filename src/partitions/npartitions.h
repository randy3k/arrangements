#ifndef NPARTITIONS_H__
#define NPARTITIONS_H__


double n_partitions(int n) {
    if (n == 0) return 1;
    // find P(1),...,P(n) sequentially
    int i, j, k, s;
    double out;
    double* p = (double*) malloc((n+1) * sizeof(double));
    p[0] = p[1] = 1;
    for(i=2 ; i<=n ; i++){
        p[i] = 0;
        for (j=1, k=1, s=1; i-j>=0; k+=3, j+=k, s=-s) {
            p[i] += s*p[i-j];
        }
        for (j=2, k=2, s=1; i-j>=0; k+=3, j+=k, s=-s) {
            p[i] += s*p[i-j];
        }
    }
    out = p[n];
    free(p);
    return out;
}


void n_partitions_bigz(mpz_t z, int n) {
    // find P(1),...,P(n) sequentially
    if (n == 0) {
        mpz_set_ui(z, 1);
        return;
    }
    int i, j, h, s;
    mpz_t* p = (mpz_t*) malloc((n+1) * sizeof(mpz_t));
    for (i=0; i<n+1; i++) mpz_init(p[i]);

    mpz_set_ui(p[0], 1);
    mpz_set_ui(p[1], 1);
    for(i=2 ; i<=n ; i++){
        for (j=1, h=1, s=1; i-j>=0; h+=3, j+=h, s=-s) {
            if (s > 0){
                mpz_add(p[i], p[i], p[i-j]);
            } else {
                mpz_sub(p[i], p[i], p[i-j]);
            }
        }
        for (j=2, h=2, s=1; i-j>=0; h+=3, j+=h, s=-s) {
            if (s > 0){
                mpz_add(p[i], p[i], p[i-j]);
            } else {
                mpz_sub(p[i], p[i], p[i-j]);
            }
        }
    }
    mpz_set(z, p[n]);
    for (i=0; i<n+1; i++) mpz_clear(p[i]);
    free(p);
}


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

double n_k_max_partitions(int n, int k, int m) {
    return nkm(n-k, k, m-1);
}

void n_k_max_partitions_bigz(mpz_t z, int n, int k, int m) {
    nkm_bigz(z, n-k, k, m-1);
}

double n_min_partitions(int n, int m) {
    if (n == 0) {
        return 1;
    }

    int i, j, h;
    double* p = (double*) malloc((n + 1) * sizeof(double));
    for (j = 1; j <= n; j++) {
        p[j] = 0;
    }
    p[0] = 1;
    for (i = m; i <= n; i++) {
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

void n_min_partitions_bigz(mpz_t z, int n, int m) {
    if (n == 0) {
        mpz_set_ui(z, 1);
    }

    int i, j, h;
    mpz_t* p = (mpz_t*) malloc((n+1) * sizeof(mpz_t));
    for (j = 0; j <= n; j++) mpz_init(p[j]);
    for (j = 1; j <= n; j++) {
        mpz_set_ui(p[j], 0);
    }
    mpz_set_ui(p[0], 1);
    for (i = m; i <= n; i++) {
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

#endif /* end of include guard: NPARTITIONS_H__ */
