// #define DEBUG 1

#include "rcpp_hsblock.hh"

// [[Rcpp::export]]
RcppExport SEXP vem_hsb_bern(SEXP adj_sexp, SEXP z_sexp, SEXP options_sexp) {
  BEGIN_RCPP

  Rcpp::traits::input_parameter<const SpMat>::type adj(adj_sexp);
  Rcpp::traits::input_parameter<const SpMat>::type latent_init(z_sexp);
  Rcpp::List options_list(options_sexp);
  options_t opt;
  set_options_from_list(options_list, opt);

  using TD = hsb_bern_t<Scalar>;
  using Tree = btree_t<TD>;
  using UD = hsb_update_data_t<Tree, NON_DEGREE_CORRECTED>;
  if(!valid_bern_data(adj)) {
    ELOG("Invalid adjacency matrix");
    return Rcpp::List::create();
  }

  auto ret = var_em<Tree, UD>(adj, latent_init, opt);
  return Rcpp::wrap(ret);
  END_RCPP
}

// [[Rcpp::export]]
RcppExport SEXP vem_hsb_pois(SEXP adj_sexp, SEXP z_sexp, SEXP options_sexp) {
  BEGIN_RCPP

  Rcpp::traits::input_parameter<const SpMat>::type adj(adj_sexp);
  Rcpp::traits::input_parameter<const SpMat>::type latent_init(z_sexp);
  Rcpp::List options_list(options_sexp);
  options_t opt;
  set_options_from_list(options_list, opt);

  using TD = hsb_pois_t<Scalar>;
  using Tree = btree_t<TD>;
  using UD = hsb_update_data_t<Tree, NON_DEGREE_CORRECTED>;

  if(!valid_pois_data(adj)) {
    ELOG("Invalid adjacency matrix");
    return Rcpp::List::create();
  }

  auto ret = var_em<Tree, UD>(adj, latent_init, opt);
  return Rcpp::wrap(ret);
  END_RCPP
}

// [[Rcpp::export]]
RcppExport SEXP vem_dhsb_pois(SEXP adj_sexp, SEXP z_sexp, SEXP options_sexp) {
  BEGIN_RCPP

  Rcpp::traits::input_parameter<const SpMat>::type adj(adj_sexp);
  Rcpp::traits::input_parameter<const SpMat>::type latent_init(z_sexp);
  Rcpp::List options_list(options_sexp);
  options_t opt;
  set_options_from_list(options_list, opt);

  using TD = hsb_pois_t<Scalar>;
  using Tree = btree_t<TD>;
  using UD = hsb_update_data_t<Tree, DEGREE_CORRECTED>;

  if(!valid_pois_data(adj)) {
    ELOG("Invalid adjacency matrix");
    return Rcpp::List::create();
  }

  auto ret = var_em<Tree, UD>(adj, latent_init, opt);
  return Rcpp::wrap(ret);
  END_RCPP
}
