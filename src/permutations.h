#ifndef PERMUTATIONS_H__
#define PERMUTATIONS_H__

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP next_k_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout);
SEXP num_k_permutations(SEXP _n, SEXP _k);
SEXP num_k_permutations_bigz(SEXP _n, SEXP _k);
SEXP get_k_permutations(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_multiset_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout);
SEXP num_multiset_permutations(SEXP freq, SEXP _k);
SEXP num_multiset_permutations_bigz(SEXP freq, SEXP _k);
SEXP get_multiset_permutation(SEXP freq, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_ordinary_permutations(SEXP _n, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout);
SEXP num_multiset_n_permutations(SEXP freq);
SEXP num_multiset_n_permutations_bigz(SEXP freq);
SEXP get_ordinary_permutations(SEXP _n, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_replacement_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout);
SEXP get_replacement_permutations(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);


#endif /* end of include guard: PERMUTATIONS_H__ */
