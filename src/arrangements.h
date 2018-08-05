#ifndef _ARRANGEMENTS_H_
#define _ARRANGEMENTS_H_ 1

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP next_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type);

SEXP next_asc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _type);

SEXP next_desc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _type);

SEXP npart_k(SEXP _n, SEXP _k);

SEXP npart_k_bigz(SEXP _n, SEXP _k);

SEXP next_k_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type);

SEXP nperm_k(SEXP _n, SEXP _k);

SEXP nperm_k_bigz(SEXP _n, SEXP _k);

SEXP next_multiset_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _type);

SEXP nperm_f(SEXP freq, SEXP _k);

SEXP nperm_f_bigz(SEXP freq, SEXP _k);

SEXP next_multiset_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _type);

SEXP ncomb_f(SEXP freq, SEXP _k);

SEXP ncomb_f_bigz(SEXP freq, SEXP _k);

SEXP next_asc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _type);

SEXP next_desc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _type);

SEXP npart(SEXP _n);

SEXP npart_bigz(SEXP _n);

SEXP next_permutations(SEXP _n, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _type);

SEXP nperm_n(SEXP freq);

SEXP nperm_n_bigz(SEXP freq);

SEXP next_replace_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type);

SEXP next_replace_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type);

#endif
