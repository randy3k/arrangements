#ifndef PARTITIONS_H__
#define PARTITIONS_H__

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP npartitions(SEXP _n, SEXP _k, SEXP _bigz);
SEXP collect_partitions(SEXP _n, SEXP _k, SEXP _descending, SEXP _layout, SEXP _d, SEXP _index, SEXP _nsample, SEXP state, SEXP _skip, SEXP _drop);

#endif /* end of include guard: PARTITIONS_H__ */
