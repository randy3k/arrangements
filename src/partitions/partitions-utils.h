#ifndef PARTITIONS_UTILS_H__
#define PARTITIONS_UTILS_H__
#include <gmp.h>
#include <stdlib.h>

void make_partition_table(double* table, int n);

void make_partition_table_bigz(mpz_t* table, int n);

void make_partition_max_table(double* table, int n);

void make_partition_max_table_bigz(mpz_t* table, int n);

void make_distinct_partition_table(double* table, int n);

void make_distinct_partition_table_bigz(mpz_t* table, int n);

void make_distinct_partition_max_table(double* table, int n);

void make_distinct_partition_max_table_bigz(mpz_t* table, int n);

void make_k_partition_table(double* table, int n, int k);

void make_k_partition_table_bigz(mpz_t* table, int n, int k);

void make_nkm_table(double* table, int n_max);

void make_nkm_table_bigz(mpz_t* table, int n_max);

double n_partitions(int n);

void n_partitions_bigz(mpz_t z, int n);

double n_min_partitions(int n, int m);

void n_min_partitions_bigz(mpz_t z, int n, int m);

double n_max_partitions(int n, int m);

void n_max_partitions_bigz(mpz_t z, int n, int m);

double n_k_partitions(int n, int k);

void n_k_partitions_bigz(mpz_t z, int n, int k);

double n_k_min_partitions(int n, int k, int m);

void n_k_min_partitions_bigz(mpz_t z, int n, int k, int m);

double n_k_max_partitions(int n, int k, int m);

void n_k_max_partitions_bigz(mpz_t z, int n, int k, int m);

double n_distinct_partitions(int n);

void n_distinct_partitions_bigz(mpz_t z, int n);

double n_min_distinct_partitions(int n, int m);

void n_min_distinct_partitions_bigz(mpz_t z, int n, int m);

double n_max_distinct_partitions(int n, int m);

void n_max_distinct_partitions_bigz(mpz_t z, int n, int m);

double n_k_distinct_partitions(int n, int k);

void n_k_distinct_partitions_bigz(mpz_t z, int n, int k);


#endif /* end of include guard: PARTITIONS_UTILS_H__ */
