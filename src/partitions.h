#ifndef PARTITIONS_H__
#define PARTITIONS_H__

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP next_asc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _layout);
SEXP next_desc_k_partitions(SEXP _n, SEXP _k, SEXP _d, SEXP state, SEXP _layout);
SEXP num_k_partitions(SEXP _n, SEXP _k);
SEXP num_k_partitions_bigz(SEXP _n, SEXP _k);

SEXP next_asc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _layout);
SEXP next_desc_partitions(SEXP _n, SEXP _d, SEXP state, SEXP _layout);
SEXP num_partitions(SEXP _n);
SEXP num_partitions_bigz(SEXP _n);


#endif /* end of include guard: PARTITIONS_H__ */
