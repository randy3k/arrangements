#define RESULT_NILSXP_PART() \
    if (layout == 'r') { \
        result = PROTECT(Rf_allocMatrix(INTSXP, d, n)); \
        nprotect++; \
        int* resultp = INTEGER(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j + i*d] = ap[i]; \
            } \
            for (i=k; i<n; i++) { \
                resultp[j + i*d] = 0; \
            } \
        } \
    } else if (layout == 'c') { \
        result = PROTECT(Rf_allocMatrix(INTSXP, n, d)); \
        nprotect++; \
        int * resultp = INTEGER(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j * n + i] = ap[i]; \
            } \
            for (i=k; i<n; i++) { \
                resultp[j * n + i] = 0; \
            } \
        } \
    } else if (layout == 'l') { \
        result = PROTECT(Rf_allocVector(VECSXP, d)); \
        nprotect++; \
        SEXP resulti; \
        int* resultp; \
        for (j=0; j<d; j++) { \
            NEXT(); \
            resulti = Rf_allocVector(INTSXP, k); \
            resultp = INTEGER(resulti); \
            for (i=0; i<k; i++) { \
                resultp[i] = ap[i]; \
            } \
            SET_VECTOR_ELT(result, j, resulti); \
        } \
    }

#define RESULT_NILSXP_K_PART() \
    if (layout == 'r') { \
        result = PROTECT(Rf_allocMatrix(INTSXP, d, k)); \
        nprotect++; \
        int* resultp = INTEGER(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j + i*d] = ap[i]; \
            } \
        } \
    } else if (layout == 'c') { \
        result = PROTECT(Rf_allocMatrix(INTSXP, k, d)); \
        nprotect++; \
        int * resultp = INTEGER(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j * k + i] = ap[i]; \
            } \
        } \
    } else if (layout == 'l') { \
        result = PROTECT(Rf_allocVector(VECSXP, d)); \
        nprotect++; \
        SEXP resulti; \
        int* resultp; \
        for (j=0; j<d; j++) { \
            NEXT(); \
            resulti = Rf_allocVector(INTSXP, k); \
            resultp = INTEGER(resulti); \
            for (i=0; i<k; i++) { \
                resultp[i] = ap[i]; \
            } \
            SET_VECTOR_ELT(result, j, resulti); \
        } \
    }


#define RESULT_NILSXP(k) \
    if (layout == 'r') { \
        result = PROTECT(Rf_allocMatrix(INTSXP, d, k)); \
        nprotect++; \
        int* resultp = INTEGER(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j + i*d] = ap[i] + 1; \
            } \
        } \
    } else if (layout == 'c') { \
        result = PROTECT(Rf_allocMatrix(INTSXP, k, d)); \
        nprotect++; \
        int * resultp = INTEGER(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j * k + i] = ap[i] + 1; \
            } \
        } \
    } else if (layout == 'l') { \
        result = PROTECT(Rf_allocVector(VECSXP, d)); \
        nprotect++; \
        SEXP resulti; \
        int* resultp; \
        for (j=0; j<d; j++) { \
            NEXT(); \
            resulti = Rf_allocVector(INTSXP, k); \
            resultp = INTEGER(resulti); \
            for (i=0; i<k; i++) { \
                resultp[i] = ap[i] + 1; \
            } \
            SET_VECTOR_ELT(result, j, resulti); \
        } \
    }

#define RESULT_INTSXP(k) \
    int* labelsp = INTEGER(labels); \
    if (layout == 'r') { \
        result = PROTECT(Rf_allocMatrix(INTSXP, d, k)); \
        nprotect++; \
        int* resultp = INTEGER(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j + i*d] = labelsp[ap[i]]; \
            } \
        } \
    } else if (layout == 'c') { \
        result = PROTECT(Rf_allocMatrix(INTSXP, k, d)); \
        nprotect++; \
        int * resultp = INTEGER(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j * k + i] = labelsp[ap[i]]; \
            } \
        } \
    } else if (layout == 'l') { \
        result = PROTECT(Rf_allocVector(VECSXP, d)); \
        nprotect++; \
        SEXP resulti; \
        int* resultp; \
        for (j=0; j<d; j++) { \
            NEXT(); \
            resulti = Rf_allocVector(INTSXP, k); \
            resultp = INTEGER(resulti); \
            for (i=0; i<k; i++) { \
                resultp[i] = labelsp[ap[i]]; \
            } \
            SET_VECTOR_ELT(result, j, resulti); \
        } \
    }

#define RESULT_REALSXP(k) \
    double* labelsp = REAL(labels); \
    if (layout == 'r') { \
        result = PROTECT(Rf_allocMatrix(REALSXP, d, k)); \
        nprotect++; \
        double* resultp = REAL(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j + i*d] = labelsp[ap[i]]; \
            } \
        } \
    } else if (layout == 'c') { \
        result = PROTECT(Rf_allocMatrix(REALSXP, k, d)); \
        nprotect++; \
        double * resultp = REAL(result); \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                resultp[j * k + i] = labelsp[ap[i]]; \
            } \
        } \
    } else if (layout == 'l') { \
        result = PROTECT(Rf_allocVector(VECSXP, d)); \
        nprotect++; \
        SEXP resulti; \
        double* resultp; \
        for (j=0; j<d; j++) { \
            NEXT(); \
            resulti = Rf_allocVector(REALSXP, k); \
            resultp = REAL(resulti); \
            for (i=0; i<k; i++) { \
                resultp[i] = labelsp[ap[i]]; \
            } \
            SET_VECTOR_ELT(result, j, resulti); \
        } \
    }

#define RESULT_STRSXP(k) \
    if (layout == 'r') { \
        result = PROTECT(Rf_allocMatrix(STRSXP, d, k)); \
        nprotect++; \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                SET_STRING_ELT(result, j + i*d, STRING_ELT(labels, ap[i])); \
            } \
        } \
    } else if (layout == 'c') { \
        result = PROTECT(Rf_allocMatrix(STRSXP, k, d)); \
        nprotect++; \
        for (j=0; j<d; j++) { \
            NEXT(); \
            for (i=0; i<k; i++) { \
                SET_STRING_ELT(result, j * k + i, STRING_ELT(labels, ap[i])); \
            } \
        } \
    } else if (layout == 'l') { \
        result = PROTECT(Rf_allocVector(VECSXP, d)); \
        nprotect++; \
        SEXP resulti; \
        for (j=0; j<d; j++) { \
            NEXT(); \
            resulti = Rf_allocVector(STRSXP, k); \
            for (i=0; i<k; i++) { \
                SET_STRING_ELT(resulti, i, STRING_ELT(labels, ap[i])); \
            } \
            SET_VECTOR_ELT(result, j, resulti); \
        } \
    }
