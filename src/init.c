#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _skpr_AOptimality(SEXP);
extern SEXP _skpr_calcAliasTrace(SEXP, SEXP);
extern SEXP _skpr_calculateAOptimalityPseudo(SEXP);
extern SEXP _skpr_calculateDEfficiency(SEXP);
extern SEXP _skpr_covarianceMatrix(SEXP);
extern SEXP _skpr_covarianceMatrixPseudo(SEXP);
extern SEXP _skpr_DOptimality(SEXP);
extern SEXP _skpr_DOptimalityBlocked(SEXP, SEXP);
extern SEXP _skpr_genBlockedOptimalDesign(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _skpr_genOptimalDesign(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _skpr_getPseudoInverse(SEXP);
extern SEXP _skpr_IOptimality(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_skpr_AOptimality",                (DL_FUNC) &_skpr_AOptimality,                 1},
    {"_skpr_calcAliasTrace",             (DL_FUNC) &_skpr_calcAliasTrace,              2},
    {"_skpr_calculateAOptimalityPseudo", (DL_FUNC) &_skpr_calculateAOptimalityPseudo,  1},
    {"_skpr_calculateDEfficiency",       (DL_FUNC) &_skpr_calculateDEfficiency,        1},
    {"_skpr_covarianceMatrix",           (DL_FUNC) &_skpr_covarianceMatrix,            1},
    {"_skpr_covarianceMatrixPseudo",     (DL_FUNC) &_skpr_covarianceMatrixPseudo,      1},
    {"_skpr_DOptimality",                (DL_FUNC) &_skpr_DOptimality,                 1},
    {"_skpr_DOptimalityBlocked",         (DL_FUNC) &_skpr_DOptimalityBlocked,          2},
    {"_skpr_genBlockedOptimalDesign",    (DL_FUNC) &_skpr_genBlockedOptimalDesign,    13},
    {"_skpr_genOptimalDesign",           (DL_FUNC) &_skpr_genOptimalDesign,            8},
    {"_skpr_getPseudoInverse",           (DL_FUNC) &_skpr_getPseudoInverse,            1},
    {"_skpr_IOptimality",                (DL_FUNC) &_skpr_IOptimality,                 3},
    {NULL, NULL, 0}
};

void R_init_skpr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
