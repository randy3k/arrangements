#ifndef PERMUTATIONS_H__
#define PERMUTATIONS_H__

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP npermutations(SEXP _x, SEXP _k, SEXP _n, SEXP _v, SEXP _freq, SEXP _replace, SEXP _bigz);
SEXP get_permutations(SEXP _x, SEXP _k, SEXP _n, SEXP _v, SEXP _freq, SEXP _replace,
                      SEXP _layout, SEXP _d, SEXP _index, SEXP _nsample, SEXP state, SEXP _skip);

#endif /* end of include guard: PERMUTATIONS_H__ */
