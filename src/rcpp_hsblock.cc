#include "rcpp_hsblock.hh"

RcppExport SEXP rcpp_hsblock(SEXP options_sexp) {
  // SEXP Adj, SEXP effect_se_sexp, SEXP x_sexp,
  //                              SEXP c_sexp, SEXP c_delta_sexp,

  BEGIN_RCPP

  // Rcpp::traits::input_parameter<const SpMat>::type effect(effect_sexp);
  // Rcpp::traits::input_parameter<const Mat>::type effect_se(effect_se_sexp);
  // Rcpp::traits::input_parameter<const Mat>::type C(c_sexp);
  // Rcpp::traits::input_parameter<const Mat>::type Cdelta(c_delta_sexp);
  // Rcpp::traits::input_parameter<const Mat>::type X(x_sexp);
  Rcpp::List options_list(options_sexp);

  options_t opt;
  // set_options_from_list(options_list, opt);
  // return Rcpp::wrap(impl_fit_zqtl(effect, effect_se, X, C, Cdelta, opt));
  END_RCPP
}
