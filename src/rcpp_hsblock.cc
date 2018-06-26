#include "rcpp_hsblock.hh"

// adj_sexp = Matrix::sparseMatrix object

// [[Rcpp::export]]
RcppExport SEXP rcpp_hsblock(SEXP adj_sexp, SEXP z_sexp, SEXP depth_sexp) {
  // SEXP Adj, SEXP effect_se_sexp, SEXP x_sexp,
  //                              SEXP c_sexp, SEXP c_delta_sexp,

  BEGIN_RCPP

  Rcpp::traits::input_parameter<const SpMat>::type adj(adj_sexp);
  Rcpp::traits::input_parameter<Mat>::type latent_init(z_sexp);
  Rcpp::traits::input_parameter<Index>::type depth(depth_sexp);

  SpMat A = adj;
  Mat Z = latent_init;

  TLOG("depth = " << depth);

  if(A.rows() != A.cols()) {
    ELOG("Invalid adjacency matrix");
    return Rcpp::List::create();
  }

  if(A.rows() != Z.rows()) {
    ELOG("Adj and Latent should match dimensionality (in rows)");
    return Rcpp::List::create();
  }

  if(Z.cols() != (2 << (depth - 1))) {
    ELOG("Number clusters should be 2^(depth - 1)");
    return Rcpp::List::create();
  }

  // check

  // Rcpp::traits::input_parameter<const Mat>::type C(c_sexp);
  // Rcpp::traits::input_parameter<const Mat>::type Cdelta(c_delta_sexp);
  // Rcpp::traits::input_parameter<const Mat>::type X(x_sexp);

  // Rcpp::List options_list(options_sexp);
  // options_t opt;

  // set_options_from_list(options_list, opt);
  // return Rcpp::wrap(impl_fit_zqtl(effect, effect_se, X, C, Cdelta, opt));

  using tree_data_t = hsb_bern_t<Scalar>;
  using model_tree_t = btree_t<tree_data_t>;
  using update_t = hsb_update_func_t<model_tree_t>;

  // For each newtok, we pre-calculate cluster-specific degree matrix
  //

  Mat C = A * Z;
  model_tree_t tree(depth);
  update_t update_func(C, Z);

  auto root = tree.root_node_obj();

  // 1. initialize the tree model
  update_func.clear_tree_data(root);

  Index ii = 0;

  update_func.collect_delta_stat(root, ii);

  // 2. initialize the latent variable matrix
  // 3. calculate cluster-specific deree matrix
  // 4. for each vertex: update cluster membership
  //   a. remove current assignment
  //   b. delta-update cluster-specific degree, bottom
  //   c. calculate assignment score
  //   d. make a change
  //   e. delta-update cluster-specific degree, bottom
  // 5. update the tree model
  //   a. take gradient
  //   b. update using ADAM?

  END_RCPP

  return Rcpp::List::create();
}
