#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "combinations.h"
#include "permutations.h"
#include "partitions.h"


static const R_CallMethodDef CallEntries[] = {
    {"ncombinations", (DL_FUNC) &ncombinations, 7},
    {"collect_combinations", (DL_FUNC) &collect_combinations, 13},

    {"npermutations", (DL_FUNC) &npermutations, 7},
    {"collect_permutations", (DL_FUNC) &collect_permutations, 13},

    {"npartitions", (DL_FUNC) &npartitions, 3},
    {"collect_partitions", (DL_FUNC) &collect_partitions, 10},

    {NULL, NULL, 0}
};

void R_init_arrangements(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
