// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// LRTq
double LRTq(NumericVector expr, IntegerMatrix geno, NumericVector causal_ratio, int perm);
RcppExport SEXP _LRTq_LRTq(SEXP exprSEXP, SEXP genoSEXP, SEXP causal_ratioSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type expr(exprSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type geno(genoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type causal_ratio(causal_ratioSEXP);
    Rcpp::traits::input_parameter< int >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(LRTq(expr, geno, causal_ratio, perm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LRTq_LRTq", (DL_FUNC) &_LRTq_LRTq, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_LRTq(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
