#ifndef _M_PERM_H
#define _M_PERM_H 1

#include <stddef.h>
#include "utils.h"

// A Simple, Efficient P(n,k) Algorithm by Alistair Israel
// http://alistairisrael.wordpress.com/2009/09/22/simple-efficient-pnk-algorithm/
// c implementation by 2017 Randy Lai
// http://randycity.github.io

unsigned int next_multiset_permutation(unsigned int *ar, size_t n, size_t k)
{
    unsigned int i;
    unsigned int j;
    unsigned int edge = k-1;

    if(k<n){
        j = k;
        // search for largest j such that nth_j > nth_edge (a is increasing for j>=k)
        while(j<n && ar[edge]>=ar[j]) j++;
    }
    if(k<n && j<n){
        swap(ar, edge, j);
    }else{
        if (k<n){
            reverse(ar+k, n-k);
        }

        // find rightmost ascent to left of edge
        for (i = edge -1; ; i--) {
            if (ar[i] < ar[i+1]) {
                break;
            } else if (i == 0) {
                return 0;
            }
        }

        // find smallest j>=i+1 where nth_j>nth_i (a is decreasing for j>=i+1)
        j = n-1;
        while(j>i && ar[i] >= ar[j]) j--;

        swap(ar, i, j);

        reverse(ar+i+1, n-i-1);
    }

    return 1;
}


#endif
