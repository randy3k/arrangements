#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "combinations.h"
#include "permutations.h"
#include "partitions.h"
#include "utils.h"

static const R_CallMethodDef CallEntries[] = {
    {"ncombinations", (DL_FUNC) &ncombinations, 7},
    {"get_combinations", (DL_FUNC) &get_combinations, 13},

    {"npermutations", (DL_FUNC) &npermutations, 7},
    {"get_permutations", (DL_FUNC) &get_permutations, 13},

    {"next_asc_k_partitions", (DL_FUNC) &next_asc_k_partitions, 5},
    {"next_desc_k_partitions", (DL_FUNC) &next_desc_k_partitions, 5},
    {"num_k_partitions", (DL_FUNC) &num_k_partitions, 2},
    {"num_k_partitions_bigz", (DL_FUNC) &num_k_partitions_bigz, 2},

    {"next_asc_partitions", (DL_FUNC) &next_asc_partitions, 4},
    {"next_desc_partitions", (DL_FUNC) &next_desc_partitions, 4},
    {"num_partitions", (DL_FUNC) &num_partitions, 1},
    {"num_partitions_bigz", (DL_FUNC) &num_partitions_bigz, 1},

    {NULL, NULL, 0}
};

void R_init_arrangements(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
