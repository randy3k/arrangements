#ifndef _UTILS_H_
#define _UTILS_H_ 1

#include <R.h>
#include <Rinternals.h>

SEXP resize_row(SEXP x, size_t n, size_t k, size_t d);

SEXP resize_col(SEXP x, size_t n, size_t k, size_t d);

SEXP resize_list(SEXP x, size_t k, size_t d);

int as_uint(SEXP x);

SEXP as_uint_array(SEXP x);

double fallfact(size_t n, size_t k);

double choose(size_t n, size_t k);

double multichoose(int* f, size_t flen);

#endif
