#ifndef _PART_H
#define _PART_H 1

#include <stddef.h>

unsigned int next_asc_partition(unsigned int *ar, size_t* k);

unsigned int next_desc_partition(unsigned int *ar, size_t* hp, size_t* kp);

unsigned int next_asc_k_partition(unsigned int *ar, size_t n, unsigned int k);

unsigned int next_desc_k_partition(unsigned int *ar, size_t n, unsigned int k);

#endif

#ifndef _K_PART_H
#define _K_PART_H 1

#include <stddef.h>

unsigned int next_asc_k_partition(unsigned int *ar, size_t n, unsigned int k);

unsigned int next_desc_k_partition(unsigned int *ar, size_t n, unsigned int k);

#endif
