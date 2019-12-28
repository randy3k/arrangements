#ifndef COMPOSITIONS_UTILS_H
#define COMPOSITIONS_UTILS_H

#include <gmp.h>
#include <math.h>
#include <stdlib.h>

double n_compositions(int n);

void n_compositions_bigz(mpz_t z, int n);

double n_k_compositions(int n, int k);

void n_k_compositions_bigz(mpz_t z, int n, int k);


#endif /* end of include guard: COMPOSITIONS_UTILS_H */
