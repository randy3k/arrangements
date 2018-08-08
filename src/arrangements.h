#ifndef _ARRANGEMENTS_H_
#define _ARRANGEMENTS_H_ 1

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP next_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout);

SEXP get_combinations(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_asc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _layout);

SEXP next_desc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _layout);

SEXP num_k_partitions(SEXP _n, SEXP _k);

SEXP num_k_partitions_bigz(SEXP _n, SEXP _k);

SEXP next_k_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout);

SEXP num_k_permutations(SEXP _n, SEXP _k);

SEXP num_k_permutations_bigz(SEXP _n, SEXP _k);

SEXP get_k_permutations(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_multiset_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout);

SEXP num_multiset_permutations(SEXP freq, SEXP _k);

SEXP num_multiset_permutations_bigz(SEXP freq, SEXP _k);

SEXP get_multiset_permutation(SEXP freq, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_multiset_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout);

SEXP num_multiset_combinations(SEXP freq, SEXP _k);

SEXP num_multiset_combinations_bigz(SEXP freq, SEXP _k);

SEXP get_multiset_combination(SEXP freq, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_asc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _layout);

SEXP next_desc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _layout);

SEXP num_partitions(SEXP _n);

SEXP num_partitions_bigz(SEXP _n);

SEXP next_permutations(SEXP _n, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout);

SEXP num_multiset_n_permutations(SEXP freq);

SEXP num_multiset_n_permutations_bigz(SEXP freq);

SEXP get_permutation(SEXP _n, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_replacement_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout);

SEXP get_replacement_combination(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_replacement_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout);

SEXP get_replacement_permutation(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

#endif
