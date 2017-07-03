// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// selinv2julia
int selinv2julia(int nnodes, int nnz, IntegerVector colptr_, IntegerVector rowind_, NumericVector nzval_);
RcppExport SEXP sInverse_selinv2julia(SEXP nnodesSEXP, SEXP nnzSEXP, SEXP colptr_SEXP, SEXP rowind_SEXP, SEXP nzval_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nnodes(nnodesSEXP);
    Rcpp::traits::input_parameter< int >::type nnz(nnzSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type colptr_(colptr_SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rowind_(rowind_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nzval_(nzval_SEXP);
    rcpp_result_gen = Rcpp::wrap(selinv2julia(nnodes, nnz, colptr_, rowind_, nzval_));
    return rcpp_result_gen;
END_RCPP
}
