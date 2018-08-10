#ifndef COMBINATIONS_H__
#define COMBINATIONS_H__

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP ncombinations(SEXP _x, SEXP _k, SEXP _n, SEXP _v, SEXP _freq, SEXP _replace, SEXP _bigz);
SEXP get_combinations(SEXP _x, SEXP _k, SEXP _n, SEXP _v, SEXP _freq, SEXP _replace,
                      SEXP _layout, SEXP _d, SEXP _index, SEXP _nsample, SEXP state, SEXP _skip);

#endif /* end of include guard: COMBINATIONS_H__ */
