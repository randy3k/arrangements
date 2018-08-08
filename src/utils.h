#ifndef _UTILS_H_
#define _UTILS_H_ 1

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <gmp.h>

SEXP resize_row(SEXP x, size_t n, size_t k, size_t d);

SEXP resize_col(SEXP x, size_t n, size_t k, size_t d);

SEXP resize_list(SEXP x, size_t k, size_t d);

SEXP resize_layout(SEXP x, size_t d, char layout);

void attach_factor_levels(SEXP result, SEXP lables);

char layout_flag(SEXP _layout);

int verify_dimension(double dd, int k, char layout);

int variable_exist(SEXP state, char* name, int TYPE, int k, void** p);

int as_uint(SEXP x);

int* as_uint_array(SEXP x);

double fact(size_t n);

double fallfact(size_t n, size_t k);

double choose(size_t n, size_t k);

double multichoose(int* f, size_t flen);

SEXP verify_parameters(SEXP _n, SEXP _k, SEXP _x, SEXP _freq);

void set_gmp_randstate(gmp_randstate_t randstate);

#endif
