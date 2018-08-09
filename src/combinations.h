#ifndef COMBINATIONS_H__
#define COMBINATIONS_H__

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP next_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout);
SEXP get_combinations(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_multiset_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP freq, SEXP _layout);
SEXP num_multiset_combinations(SEXP freq, SEXP _k);
SEXP num_multiset_combinations_bigz(SEXP freq, SEXP _k);
SEXP get_multiset_combination(SEXP freq, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

SEXP next_replacement_combinations(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP labels, SEXP _layout);
SEXP get_replacement_combination(SEXP _n, SEXP _k, SEXP labels, SEXP _layout, SEXP _index, SEXP _nsample);

#endif /* end of include guard: COMBINATIONS_H__ */
