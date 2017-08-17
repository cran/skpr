// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// AOptimality
double AOptimality(const arma::mat& currentDesign);
RcppExport SEXP _skpr_AOptimality(SEXP currentDesignSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type currentDesign(currentDesignSEXP);
    rcpp_result_gen = Rcpp::wrap(AOptimality(currentDesign));
    return rcpp_result_gen;
END_RCPP
}
// calcAliasTrace
double calcAliasTrace(const arma::mat& currentDesign, const arma::mat& aliasMatrix);
RcppExport SEXP _skpr_calcAliasTrace(SEXP currentDesignSEXP, SEXP aliasMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type currentDesign(currentDesignSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aliasMatrix(aliasMatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(calcAliasTrace(currentDesign, aliasMatrix));
    return rcpp_result_gen;
END_RCPP
}
// calculateAOptimalityPseudo
double calculateAOptimalityPseudo(const arma::mat& currentDesign);
RcppExport SEXP _skpr_calculateAOptimalityPseudo(SEXP currentDesignSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type currentDesign(currentDesignSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateAOptimalityPseudo(currentDesign));
    return rcpp_result_gen;
END_RCPP
}
// calculateDEfficiency
double calculateDEfficiency(const arma::mat& currentDesign);
RcppExport SEXP _skpr_calculateDEfficiency(SEXP currentDesignSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type currentDesign(currentDesignSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateDEfficiency(currentDesign));
    return rcpp_result_gen;
END_RCPP
}
// covarianceMatrix
arma::mat covarianceMatrix(const arma::mat& design);
RcppExport SEXP _skpr_covarianceMatrix(SEXP designSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type design(designSEXP);
    rcpp_result_gen = Rcpp::wrap(covarianceMatrix(design));
    return rcpp_result_gen;
END_RCPP
}
// covarianceMatrixPseudo
arma::mat covarianceMatrixPseudo(const arma::mat& design);
RcppExport SEXP _skpr_covarianceMatrixPseudo(SEXP designSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type design(designSEXP);
    rcpp_result_gen = Rcpp::wrap(covarianceMatrixPseudo(design));
    return rcpp_result_gen;
END_RCPP
}
// DOptimality
double DOptimality(const arma::mat& currentDesign);
RcppExport SEXP _skpr_DOptimality(SEXP currentDesignSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type currentDesign(currentDesignSEXP);
    rcpp_result_gen = Rcpp::wrap(DOptimality(currentDesign));
    return rcpp_result_gen;
END_RCPP
}
// DOptimalityBlocked
double DOptimalityBlocked(const arma::mat& currentDesign, const arma::mat& blockedVar);
RcppExport SEXP _skpr_DOptimalityBlocked(SEXP currentDesignSEXP, SEXP blockedVarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type currentDesign(currentDesignSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type blockedVar(blockedVarSEXP);
    rcpp_result_gen = Rcpp::wrap(DOptimalityBlocked(currentDesign, blockedVar));
    return rcpp_result_gen;
END_RCPP
}
// genOptimalDesign
List genOptimalDesign(arma::mat initialdesign, const arma::mat& candidatelist, const std::string condition, const arma::mat& momentsmatrix, NumericVector initialRows, arma::mat aliasdesign, const arma::mat& aliascandidatelist, double minDopt);
RcppExport SEXP _skpr_genOptimalDesign(SEXP initialdesignSEXP, SEXP candidatelistSEXP, SEXP conditionSEXP, SEXP momentsmatrixSEXP, SEXP initialRowsSEXP, SEXP aliasdesignSEXP, SEXP aliascandidatelistSEXP, SEXP minDoptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type initialdesign(initialdesignSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type candidatelist(candidatelistSEXP);
    Rcpp::traits::input_parameter< const std::string >::type condition(conditionSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type momentsmatrix(momentsmatrixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initialRows(initialRowsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aliasdesign(aliasdesignSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aliascandidatelist(aliascandidatelistSEXP);
    Rcpp::traits::input_parameter< double >::type minDopt(minDoptSEXP);
    rcpp_result_gen = Rcpp::wrap(genOptimalDesign(initialdesign, candidatelist, condition, momentsmatrix, initialRows, aliasdesign, aliascandidatelist, minDopt));
    return rcpp_result_gen;
END_RCPP
}
// genBlockedOptimalDesign
List genBlockedOptimalDesign(arma::mat initialdesign, arma::mat candidatelist, const arma::mat& blockeddesign, const std::string condition, const arma::mat& momentsmatrix, IntegerVector initialRows, const arma::mat& blockedVar, const arma::mat& aliasdesign, const arma::mat& aliascandidatelist, double minDopt, List interactions, const arma::mat& disallowed, const bool anydisallowed);
RcppExport SEXP _skpr_genBlockedOptimalDesign(SEXP initialdesignSEXP, SEXP candidatelistSEXP, SEXP blockeddesignSEXP, SEXP conditionSEXP, SEXP momentsmatrixSEXP, SEXP initialRowsSEXP, SEXP blockedVarSEXP, SEXP aliasdesignSEXP, SEXP aliascandidatelistSEXP, SEXP minDoptSEXP, SEXP interactionsSEXP, SEXP disallowedSEXP, SEXP anydisallowedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type initialdesign(initialdesignSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type candidatelist(candidatelistSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type blockeddesign(blockeddesignSEXP);
    Rcpp::traits::input_parameter< const std::string >::type condition(conditionSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type momentsmatrix(momentsmatrixSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type initialRows(initialRowsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type blockedVar(blockedVarSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aliasdesign(aliasdesignSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aliascandidatelist(aliascandidatelistSEXP);
    Rcpp::traits::input_parameter< double >::type minDopt(minDoptSEXP);
    Rcpp::traits::input_parameter< List >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type disallowed(disallowedSEXP);
    Rcpp::traits::input_parameter< const bool >::type anydisallowed(anydisallowedSEXP);
    rcpp_result_gen = Rcpp::wrap(genBlockedOptimalDesign(initialdesign, candidatelist, blockeddesign, condition, momentsmatrix, initialRows, blockedVar, aliasdesign, aliascandidatelist, minDopt, interactions, disallowed, anydisallowed));
    return rcpp_result_gen;
END_RCPP
}
// getPseudoInverse
arma::mat getPseudoInverse(const arma::mat& currentDesign);
RcppExport SEXP _skpr_getPseudoInverse(SEXP currentDesignSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type currentDesign(currentDesignSEXP);
    rcpp_result_gen = Rcpp::wrap(getPseudoInverse(currentDesign));
    return rcpp_result_gen;
END_RCPP
}
// IOptimality
double IOptimality(const arma::mat& currentDesign, const arma::mat& momentsMatrix, const arma::mat& blockedVar);
RcppExport SEXP _skpr_IOptimality(SEXP currentDesignSEXP, SEXP momentsMatrixSEXP, SEXP blockedVarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type currentDesign(currentDesignSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type momentsMatrix(momentsMatrixSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type blockedVar(blockedVarSEXP);
    rcpp_result_gen = Rcpp::wrap(IOptimality(currentDesign, momentsMatrix, blockedVar));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_skpr_AOptimality", (DL_FUNC) &_skpr_AOptimality, 1},
    {"_skpr_calcAliasTrace", (DL_FUNC) &_skpr_calcAliasTrace, 2},
    {"_skpr_calculateAOptimalityPseudo", (DL_FUNC) &_skpr_calculateAOptimalityPseudo, 1},
    {"_skpr_calculateDEfficiency", (DL_FUNC) &_skpr_calculateDEfficiency, 1},
    {"_skpr_covarianceMatrix", (DL_FUNC) &_skpr_covarianceMatrix, 1},
    {"_skpr_covarianceMatrixPseudo", (DL_FUNC) &_skpr_covarianceMatrixPseudo, 1},
    {"_skpr_DOptimality", (DL_FUNC) &_skpr_DOptimality, 1},
    {"_skpr_DOptimalityBlocked", (DL_FUNC) &_skpr_DOptimalityBlocked, 2},
    {"_skpr_genOptimalDesign", (DL_FUNC) &_skpr_genOptimalDesign, 8},
    {"_skpr_genBlockedOptimalDesign", (DL_FUNC) &_skpr_genBlockedOptimalDesign, 13},
    {"_skpr_getPseudoInverse", (DL_FUNC) &_skpr_getPseudoInverse, 1},
    {"_skpr_IOptimality", (DL_FUNC) &_skpr_IOptimality, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_skpr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
