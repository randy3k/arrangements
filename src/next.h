#ifndef _NEXT_H_
#define _NEXT_H_ 1

unsigned int next_cartesian_product(unsigned int *ar, size_t len, const size_t *sizes);

unsigned int next_combination(unsigned int *ar, size_t n, unsigned int k);

unsigned int next_asc_k_partition(unsigned int *ar, size_t n, unsigned int k);
unsigned int next_desc_k_partition(unsigned int *ar, size_t n, unsigned int k);

unsigned int next_k_permutation(unsigned int *ar, unsigned int *cycle, size_t n, size_t k);

unsigned int next_multicombination(unsigned int *ar, size_t n, unsigned int k);

unsigned int next_multiset_combination(
    const unsigned int *multiset, unsigned int *ar, size_t n, unsigned int k);

unsigned int next_multiset_permutation(unsigned int *ar, size_t n, size_t k);

unsigned int next_asc_partition(unsigned int *ar, int* k);

unsigned int next_desc_partition(unsigned int *ar, int* hp, int* kp);

unsigned int next_permutation(unsigned int *ar, size_t len);

#endif
