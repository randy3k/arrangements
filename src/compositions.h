#ifndef COMPOSITIONS_H__
#define COMPOSITIONS_H__

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP ncompositions(SEXP _n, SEXP _k, SEXP _bigz);
SEXP get_compositions(SEXP _n, SEXP _k, SEXP _descending, SEXP _layout, SEXP _d, SEXP _index, SEXP _nsample, SEXP state, SEXP _skip, SEXP _drop);

#endif /* end of include guard: COMPOSITIONS_H__ */
