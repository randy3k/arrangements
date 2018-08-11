#ifndef GMP_UTILS_H__
#define GMP_UTILS_H__

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>

SEXP mpz_to_bigz1(mpz_t z);

int as_mpz_array(mpz_t* a, size_t n, SEXP x);

#endif /* end of include guard: GMP_UTILS_H__ */
