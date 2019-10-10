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

int variable_exists(SEXP state, char* name, int TYPE, int k, void** p);

int as_uint(SEXP x);

int* as_uint_array(SEXP x);

int* as_uint_index(SEXP x);

double fact(int n);

double fallfact(int n, int k);

double choose(int n, int k);

double multichoose(int* f, size_t flen);

SEXP mpz_to_bigz1(mpz_t z);

int as_mpz_array(mpz_t* a, size_t n, SEXP x);

void set_gmp_randstate(gmp_randstate_t randstate);

int index_length(SEXP _index);

#endif
