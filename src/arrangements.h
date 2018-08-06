#ifndef _ARRANGEMENTS_H_
#define _ARRANGEMENTS_H_ 1

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP next_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type);

SEXP get_ith_combination(SEXP _n, SEXP _k, SEXP _index);

SEXP next_asc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _type);

SEXP next_desc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _type);

SEXP num_k_partitions(SEXP _n, SEXP _k);

SEXP num_k_partitions_bigz(SEXP _n, SEXP _k);

SEXP next_k_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type);

SEXP num_k_permutations(SEXP _n, SEXP _k);

SEXP num_k_permutations_bigz(SEXP _n, SEXP _k);

SEXP get_ith_k_permutation(SEXP _n, SEXP _k, SEXP _index);

SEXP next_multiset_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _type);

SEXP num_multiset_permutations(SEXP freq, SEXP _k);

SEXP num_multiset_permutations_bigz(SEXP freq, SEXP _k);

SEXP get_ith_multiset_permutation(SEXP freq, SEXP _k, SEXP _index);

SEXP next_multiset_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _type);

SEXP num_multiset_combinations(SEXP freq, SEXP _k);

SEXP num_multiset_combinations_bigz(SEXP freq, SEXP _k);

SEXP next_asc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _type);

SEXP next_desc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _type);

SEXP num_partitions(SEXP _n);

SEXP num_partitions_bigz(SEXP _n);

SEXP next_permutations(SEXP _n, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _type);

SEXP num_multiset_n_permutations(SEXP freq);

SEXP num_multiset_n_permutations_bigz(SEXP freq);

SEXP get_ith_permutation(SEXP _n, SEXP _index);

SEXP next_replacement_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type);

SEXP next_replacement_permutations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _type);

SEXP get_ith_replacement_permutation(SEXP _n, SEXP _k, SEXP _index);

#endif
