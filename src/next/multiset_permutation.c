#include "multiset_permutation.h"

// A Simple, Efficient P(n,k) Algorithm by Alistair Israel
// http://alistairisrael.wordpress.com/2009/09/22/simple-efficient-pnk-algorithm/
// c implementation by 2017 Randy Lai
// http://randycity.github.io

static void swap(unsigned int *ar, unsigned int first, unsigned int second)
{
    unsigned int temp = ar[first];
    ar[first] = ar[second];
    ar[second] = temp;
}

static void reverse(unsigned int *ar, size_t len)
{
    unsigned int i, j;

    for (i = 0, j = len - 1; i < j; i++, j--) {
        swap(ar, i, j);
    }
}

unsigned int next_multiset_permutation(unsigned int *ar, size_t n, size_t k)
{
    long i;
    long j;
    long edge = k-1;

    if(k<n){
        j = k;
        // search for largest j such that a_j > a_edge (a is increasing for j>=k)
        while(j<n && ar[edge]>=ar[j]) j++;
    }
    if(k<n && j<n){
        swap(ar, edge, j);
    }else{
        if (k<n){
            reverse(ar+k, n-k);
        }

        // find rightmost ascent to left of edge
        i = edge -1;
        while(i>=0 && ar[i]>=ar[i+1]) i--;

        if (i<0) return 0;

        // find smallest j>=i+1 where a_j>a_i (a is decreasing for j>=i+1)
        j = n-1;
        while(j>i && ar[i] >= ar[j]) j--;

        swap(ar, i, j);

        reverse(ar+i+1, n-i-1);
    }

    return 1;
}
