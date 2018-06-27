#define DEBUG 1

#include "rcpp_hsblock.hh"

// adj_sexp = Matrix::sparseMatrix object

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

  // return Rcpp::wrap(impl_fit_zqtl(effect, effect_se, X, C, Cdelta, opt));

  using tree_data_t = hsb_bern_t<Scalar>;
  using model_tree_t = btree_t<tree_data_t>;
  using update_t = hsb_update_func_t<model_tree_t>;

  // For each newtok, we pre-calculate cluster-specific degree matrix
  //

  Index n = A.rows();
  Mat C = Z * A;
  Vec N = Z * Mat::Ones(n, 1);
  model_tree_t tree(depth - 1);
  update_t update_func(C, Z, N);

  if (tree.num_leaves() != Z.rows()) {
    ELOG("Tree and latent matrix do match");
    return Rcpp::List::create();
  }

  auto root = tree.root_node_obj();

  std::mt19937 rng(opt.rseed());
  const Index K = Z.rows();
  discrete_sampler_t<Vec> discrete(K);

  // 1. Initialize the tree model
  update_func.clear_tree_data(root);

  for (Index ii = 0; ii < n; ++ii) {
    update_func.increase_tree_stat(root, ii);
  }

  update_func.dump_tree_data(root);

  // 2. For each vertex: update cluster membership
  running_stat_t<Mat> Zstat(K, n);

  // Must keep this progress obj; otherwise segfault will occur
  Progress prog(opt.vbiter(), !opt.verbose());
  Vec llik(opt.vbiter());

  for (Index iter = 0; iter < opt.vbiter(); ++iter) {

    if (Progress::check_abort()) {
      break;
    }
    prog.increment();

    for (Index ii = 0; ii < n; ++ii) {
      // a. Remove ii's contribution on the tree
      update_func.decrease_tree_stat(root, ii);

      // don't forget the stat on the neighbors

      for (SpMat::InnerIterator jt(A, ii); jt; ++jt) {
        Index jj = jt.index();
        if (jj == ii) continue;
        update_func.decrease_tree_stat(root, jj);
      }

      // b. Remove ii's contribution on C, Z, N

      C -= Z.col(ii) * A.row(ii);
      N -= Z.col(ii);
      Z.col(ii).setZero();

      update_func.increase_tree_stat(root, ii);

      // don't forget the stat on the neighbors

      for (SpMat::InnerIterator jt(A, ii); jt; ++jt) {
        Index jj = jt.index();
        if (jj == ii) continue;
        update_func.increase_tree_stat(root, jj);
      }

      // c. calculate delta score
      const Vec& log_score = update_func.eval_tree_delta_score(root, ii);

      Index kk = discrete(log_score, rng);

      // remove temporary

      for (SpMat::InnerIterator jt(A, ii); jt; ++jt) {
        Index jj = jt.index();
        if (jj == ii) continue;
        update_func.decrease_tree_stat(root, jj);
      }

      // d. Add ii's contribution on C, Z, N

      update_func.decrease_tree_stat(root, ii);

      Z(kk, ii) = 1.0;
      C += Z.col(ii) * A.row(ii);
      N += Z.col(ii);

      for (SpMat::InnerIterator jt(A, ii); jt; ++jt) {
        Index jj = jt.index();
        if (jj == ii) continue;
        update_func.increase_tree_stat(root, jj);
      }

      // f. Add ii's contribution on the tree
      update_func.increase_tree_stat(root, ii);
    }

    // Compute log-likelihood
    Scalar score = update_func.eval_tree_score(root);
    llik(iter) = score;

    TLOG("Iter = " << iter << ", Score = " << score);

    // if (opt.verbose()) conv.print(Rcpp::Rcerr);
    Zstat(Z);
  }

  update_func.dump_tree_data(root);

  Mat Zmean = Zstat.mean();

  return Rcpp::List::create(Rcpp::_["Z"] = Zmean,
			    Rcpp::_["llik"] = llik);
  END_RCPP
  // return Rcpp::wrap(impl_fit_zqtl(effect, effect_se, X, C, Cdelta, opt));
}
