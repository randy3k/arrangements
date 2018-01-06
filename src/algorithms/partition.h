#ifndef _PART_H
#define _PART_H 1

#include <stddef.h>

unsigned int next_asc_partition(unsigned int *ar, size_t* k);

unsigned int next_desc_partition(unsigned int *ar, size_t* hp, size_t* kp);

#endif
