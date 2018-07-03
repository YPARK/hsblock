// #define DEBUG 1

#include "rcpp_hsblock.hh"

// [[Rcpp::export]]
RcppExport SEXP rcpp_hsblock(const SpMat adj, SpMat latent_init,
                             SEXP options_sexp) {
  BEGIN_RCPP

  // Rcpp::traits::input_parameter<const SpMat>::type adj(adj_sexp);
  // Rcpp::traits::input_parameter<const Mat>::type latent_init(z_sexp);

  Rcpp::List options_list(options_sexp);
  const Scalar ZERO = 0.0;
  const Scalar TOL = 1e-4;

  options_t opt;
  set_options_from_list(options_list, opt);

  const Index depth = opt.tree_depth();

  SpMat A(adj.rows(), adj.cols());
  A = adj;
  SpMat Z = latent_init;

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

  // using tree_data_t = hsb_bern_t<Scalar>;
  using tree_data_t = hsb_pois_t<Scalar>;
  using model_tree_t = btree_t<tree_data_t>;
  using update_data_t = hsb_update_data_t<model_tree_t, DEGREE_CORRECTED>;
  // using null_update_t = hsb_update_data_null_t<model_tree_t>;

  model_tree_t tree(depth - 1);
  // model_tree_t null_tree(depth - 1);
  if (tree.num_leaves() != Z.rows()) {
    ELOG("Tree and latent matrix do match");
    return Rcpp::List::create();
  }

  hsb_func_dummy_t dummy_func;

  update_data_t model_update_data(tree, A, Z);
  // null_update_t null_update_func(null_tree, A, Z);

  // Initialize the tree model
  clear_tree_data(model_update_data);
  init_tree_var(model_update_data);
  // null_update_func.clear_tree_data();
  // null_update_func.init_tree_var();

  const Index n = Z.cols();

  discrete_sampler_t<Vec> randK(Z.rows());
  running_stat_t<SpMat> Zstat(Z.rows(), Z.cols());

  const Scalar rate0 = opt.rate0();
  const Scalar decay = opt.decay();
  const Scalar delay = opt.delay();
  Progress prog(opt.vbiter(), !opt.verbose());

  for (Index iter = 0; iter < opt.vbiter(); ++iter) {
    // Latent variable inference
    Zstat.reset();

    hsblock_latent_inference(n, Zstat, randK, rng, opt,
                             std::make_tuple(model_update_data),
                             std::make_tuple(dummy_func));

    Z = Zstat.mean().pruned(ZERO, TOL);

    if (Progress::check_abort()) {
      break;
    }
    prog.increment();

    const Scalar rate = rate0 * std::pow(iter + delay, decay);

    Scalar score =
        hsblock_param_inference(opt, rate, std::make_tuple(model_update_data),
                                std::make_tuple(dummy_func));

    if (opt.verbose()) {
      TLOG("Iter [" << iter << "] [" << rate << "] Score = " << score);
    }

    llik(iter) = score;
  }

  return Rcpp::List::create(Rcpp::_["Z"] = Z, Rcpp::_["llik"] = llik);
  END_RCPP
  // return Rcpp::wrap(impl_fit_zqtl(effect, effect_se, X, C, Cdelta, opt));
}
