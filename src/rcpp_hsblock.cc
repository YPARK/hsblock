#define DEBUG 1

#include "rcpp_hsblock.hh"





// [[Rcpp::export]]
RcppExport SEXP rcpp_hsblock(SEXP adj_sexp, SEXP z_sexp, SEXP options_sexp) {
  BEGIN_RCPP

  Rcpp::traits::input_parameter<const SpMat>::type adj(adj_sexp);
  Rcpp::traits::input_parameter<Mat>::type latent_init(z_sexp);
  Rcpp::List options_list(options_sexp);

  options_t opt;
  set_options_from_list(options_list, opt);

  const Index depth = opt.tree_depth();
  SpMat A = adj;
  Mat Z = latent_init;

  TLOG("depth = " << depth);

  if (A.rows() != A.cols()) {
    ELOG("Invalid adjacency matrix");
    return Rcpp::List::create();
  }

  if (A.rows() != Z.cols()) {
    ELOG("Adj and Latent should match dimensionality (in rows)");
    return Rcpp::List::create();
  }

  if (Z.rows() != (1 << (depth - 1))) {
    ELOG("Number clusters should be 2^(depth - 1)");
    return Rcpp::List::create();
  }

  std::mt19937 rng(opt.rseed());

  Vec llik(opt.vbiter());

  //  using tree_data_t = hsb_bern_t<Scalar>;
  using tree_data_t = hsb_pois_t<Scalar>;
  using model_tree_t = btree_t<tree_data_t>;
  using update_t = hsb_update_func_t<model_tree_t>;

  Mat C = Z * A;
  Vec N = Z * Mat::Ones(Z.cols(), 1);
  model_tree_t tree(depth - 1);
  if (tree.num_leaves() != Z.rows()) {
    ELOG("Tree and latent matrix do match");
    return Rcpp::List::create();
  }

  update_t update_func(tree, A, C, Z, N);

  // 1. Initialize the tree model
  update_func.clear_tree_data();
  update_func.init_tree_var();

  const Index n = Z.cols();

  discrete_sampler_t<Vec> randK(Z.rows());

  hsblock_latent_inference(n, randK, rng, opt, std::make_tuple(update_func));

  update_func.dump_tree_data();

  // // For each newtok, we pre-calculate cluster-specific degree matrix
  // //



  // // // Compute log-likelihood
  // // Scalar score = update_func.eval_tree_score(root);
  // // llik(iter) = score;

  // // if (opt.verbose()) {
  // //   TLOG("Iter = " << iter << ", Score = " << score);
  // // }

  // // This should happen out of latent variable inference

  // update_func.clear_tree_data(root);
  // for (Index ii = 0; ii < n; ++ii) {
  //   update_func.increase_tree_stat(root, ii);
  // }
  // // update_func.update_tree_var(root, 0.01);
  // // update_func.dump_tree_data(root);

  // Mat Zmean = update_func.Zstat.mean();

  return Rcpp::List::create(Rcpp::_["Z"] = Z, Rcpp::_["llik"] = llik);
  END_RCPP
  // return Rcpp::wrap(impl_fit_zqtl(effect, effect_se, X, C, Cdelta, opt));
}
