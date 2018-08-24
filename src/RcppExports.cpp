// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dm_model_estimator
Rcpp::List dm_model_estimator(IntegerMatrix Y, IntegerVector z, int iter, IntegerMatrix S, bool aggregate, bool fold_change, NumericMatrix b, NumericVector b_0);
RcppExport SEXP _mbTest_dm_model_estimator(SEXP YSEXP, SEXP zSEXP, SEXP iterSEXP, SEXP SSEXP, SEXP aggregateSEXP, SEXP fold_changeSEXP, SEXP bSEXP, SEXP b_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< bool >::type aggregate(aggregateSEXP);
    Rcpp::traits::input_parameter< bool >::type fold_change(fold_changeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b_0(b_0SEXP);
    rcpp_result_gen = Rcpp::wrap(dm_model_estimator(Y, z, iter, S, aggregate, fold_change, b, b_0));
    return rcpp_result_gen;
END_RCPP
}
// zinb_model_estimator
Rcpp::List zinb_model_estimator(IntegerMatrix Y, IntegerVector z, NumericVector s, int iter, bool DPP, double c_s, IntegerMatrix S, bool aggregate, bool fold_change, NumericMatrix b, NumericVector b_0);
RcppExport SEXP _mbTest_zinb_model_estimator(SEXP YSEXP, SEXP zSEXP, SEXP sSEXP, SEXP iterSEXP, SEXP DPPSEXP, SEXP c_sSEXP, SEXP SSEXP, SEXP aggregateSEXP, SEXP fold_changeSEXP, SEXP bSEXP, SEXP b_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type DPP(DPPSEXP);
    Rcpp::traits::input_parameter< double >::type c_s(c_sSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< bool >::type aggregate(aggregateSEXP);
    Rcpp::traits::input_parameter< bool >::type fold_change(fold_changeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b_0(b_0SEXP);
    rcpp_result_gen = Rcpp::wrap(zinb_model_estimator(Y, z, s, iter, DPP, c_s, S, aggregate, fold_change, b, b_0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mbTest_dm_model_estimator", (DL_FUNC) &_mbTest_dm_model_estimator, 8},
    {"_mbTest_zinb_model_estimator", (DL_FUNC) &_mbTest_zinb_model_estimator, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_mbTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
