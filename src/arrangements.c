#include "arrangements.h"
#include "utils.h"

static const R_CallMethodDef CallEntries[] = {
    {"next_combinations", (DL_FUNC) &next_combinations, 6},
    {"get_ith_combination", (DL_FUNC) &get_ith_combination, 3},
    {"next_asc_k_partitions", (DL_FUNC) &next_asc_k_partitions, 5},
    {"next_desc_k_partitions", (DL_FUNC) &next_desc_k_partitions, 5},
    {"num_k_partitions", (DL_FUNC) &num_k_partitions, 2},
    {"num_k_partitions_bigz", (DL_FUNC) &num_k_partitions_bigz, 2},
    {"next_k_permutations", (DL_FUNC) &next_k_permutations, 6},
    {"num_k_permutations", (DL_FUNC) &num_k_permutations, 2},
    {"num_k_permutations_bigz", (DL_FUNC) &num_k_permutations_bigz, 2},
    {"get_ith_k_permutation", (DL_FUNC) &get_ith_k_permutation, 3},
    {"next_multiset_permutations", (DL_FUNC) &next_multiset_permutations, 7},
    {"num_multiset_permutations", (DL_FUNC) &num_multiset_permutations, 2},
    {"num_multiset_permutations_bigz", (DL_FUNC) &num_multiset_permutations_bigz, 2},
    {"get_ith_multiset_permutation", (DL_FUNC) &get_ith_multiset_permutation, 3},
    {"next_multiset_combinations", (DL_FUNC) &next_multiset_combinations, 7},
    {"num_multiset_combinations", (DL_FUNC) &num_multiset_combinations, 2},
    {"num_multiset_combinations_bigz", (DL_FUNC) &num_multiset_combinations_bigz, 2},
    {"next_asc_partitions", (DL_FUNC) &next_asc_partitions, 4},
    {"next_desc_partitions", (DL_FUNC) &next_desc_partitions, 4},
    {"num_partitions", (DL_FUNC) &num_partitions, 1},
    {"num_partitions_bigz", (DL_FUNC) &num_partitions_bigz, 1},
    {"next_permutations", (DL_FUNC) &next_permutations, 6},
    {"num_multiset_n_permutations", (DL_FUNC) &num_multiset_n_permutations, 1},
    {"num_multiset_n_permutations_bigz", (DL_FUNC) &num_multiset_n_permutations_bigz, 1},
    {"get_ith_permutation", (DL_FUNC) &get_ith_permutation, 2},
    {"next_replacement_combinations", (DL_FUNC) &next_replacement_combinations, 6},
    {"next_replacement_permutations", (DL_FUNC) &next_replacement_permutations, 6},
    {"get_ith_replacement_permutation", (DL_FUNC) &get_ith_replacement_permutation, 3},
    {"as_uint_array", (DL_FUNC) &as_uint_array, 1},
    {NULL, NULL, 0}
};

void R_init_arrangements(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
