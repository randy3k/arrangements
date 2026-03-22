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

void nth_asc_partition_table(unsigned int* ar, unsigned int n, unsigned int index, double* table, int table_n);
void nth_asc_partition_table_bigz(unsigned int* ar, unsigned int n, mpz_t index, mpz_t* table, int table_n);

void nth_desc_partition_table(unsigned int* ar, unsigned int n, unsigned int index, double* table, int table_n);
void nth_desc_partition_table_bigz(unsigned int* ar, unsigned int n, mpz_t index, mpz_t* table, int table_n);

void nth_asc_k_partition_table(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index, double* table, int table_k);
void nth_asc_k_partition_table_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index, mpz_t* table, int table_k);

void nth_desc_k_partition_table(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index, double* table, int table_n);
void nth_desc_k_partition_table_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index, mpz_t* table, int table_n);

void nth_asc_distinct_partition_table(unsigned int* ar, unsigned int m, unsigned int n, unsigned int index, double* table, int table_n);
void nth_asc_distinct_partition_table_bigz(unsigned int* ar, unsigned int m, unsigned int n, mpz_t index, mpz_t* table, int table_n);

void nth_desc_distinct_partition_table(unsigned int* ar, unsigned int m, unsigned int n, unsigned int index, double* table, int table_n);
void nth_desc_distinct_partition_table_bigz(unsigned int* ar, unsigned int m, unsigned int n, mpz_t index, mpz_t* table, int table_n);

void nth_asc_k_distinct_partition_table(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index, double* table, int table_k);
void nth_asc_k_distinct_partition_table_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index, mpz_t* table, int table_k);

void nth_desc_k_distinct_partition_table(unsigned int* ar, unsigned int n, unsigned int k, unsigned int index, double* table, int table_n);
void nth_desc_k_distinct_partition_table_bigz(unsigned int* ar, unsigned int n, unsigned int k, mpz_t index, mpz_t* table, int table_n);

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
