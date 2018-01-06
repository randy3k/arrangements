#ifndef _UTILS_H_
#define _UTILS_H_ 1

#include <R.h>
#include <Rinternals.h>

SEXP resize_row(SEXP x, size_t n, size_t m, size_t d);

SEXP resize_col(SEXP x, size_t n, size_t m, size_t d);

SEXP resize_list(SEXP x, size_t m, size_t d);

#endif
