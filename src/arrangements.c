#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "arrangements.h"

static const R_CallMethodDef CallEntries[] = {
    {"next_combinations", (DL_FUNC) &next_combinations, 6},
    {"next_asc_k_partitions", (DL_FUNC) &next_asc_k_partitions, 5},
    {"next_desc_k_partitions", (DL_FUNC) &next_desc_k_partitions, 5},
    {"npart_k", (DL_FUNC) &npart_k, 2},
    {"npart_k_bigz", (DL_FUNC) &npart_k_bigz, 2},
    {"next_k_permutations", (DL_FUNC) &next_k_permutations, 7},
    {"nperm_k", (DL_FUNC) &nperm_k, 2},
    {"nperm_k_bigz", (DL_FUNC) &nperm_k_bigz, 2},
    {"nperm_f", (DL_FUNC) &nperm_f, 2},
    {"nperm_f_bigz", (DL_FUNC) &nperm_f_bigz, 2},
    {"next_multiset_combinations", (DL_FUNC) &next_multiset_combinations, 7},
    {"ncomb_f", (DL_FUNC) &ncomb_f, 2},
    {"ncomb_f_bigz", (DL_FUNC) &ncomb_f_bigz, 2},
    {"next_asc_partitions", (DL_FUNC) &next_asc_partitions, 4},
    {"next_desc_partitions", (DL_FUNC) &next_multiset_combinations, 4},
    {"npart", (DL_FUNC) &npart, 1},
    {"npart_bigz", (DL_FUNC) &npart_bigz, 1},
    {"nperm_n", (DL_FUNC) &nperm_n, 1},
    {"nperm_n_bigz", (DL_FUNC) &nperm_n_bigz, 1},
    {"next_replace_combinations", (DL_FUNC) &next_replace_combinations, 6},
    {"next_replace_permutations", (DL_FUNC) &next_replace_permutations, 6},
    {NULL, NULL, 0}
};

void R_init_xptr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
