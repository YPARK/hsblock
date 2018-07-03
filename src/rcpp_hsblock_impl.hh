#ifndef RCPP_HSBLOCK_IMPL_HH_
#define RCPP_HSBLOCK_IMPL_HH_

template <typename Tree, typename UpdateData>
Rcpp::List var_em(const SpMat adj, const SpMat latent_init,
                  const options_t& opt) {
  const Scalar ZERO = 0.0;
  const Scalar TOL = 1e-4;
  const Index depth = opt.tree_depth();

  SpMat A = adj;
  SpMat Z = latent_init;

  if (opt.verbose()) {
    TLOG("Tree depth = " << depth);
  }

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

  Tree tree(depth - 1);
  if (tree.num_leaves() != Z.rows()) {
    ELOG("Tree and latent matrix do match");
    return Rcpp::List::create();
  }

  hsb_func_dummy_t empty;
  UpdateData update_data(tree, A, Z);

  // Initialize the tree model
  clear_tree_data(update_data);
  init_tree_var(update_data);

  const Index n = Z.cols();
  discrete_sampler_t<Vec> randK(Z.rows());
  running_stat_t<SpMat> Zstat(Z.rows(), Z.cols());

  const Scalar rate0 = opt.rate0();
  const Scalar decay = opt.decay();
  const Scalar delay = opt.delay();
  Progress prog(opt.vbiter(), !opt.verbose());

  for (Index iter = 0; iter < opt.vbiter(); ++iter) {
    Zstat.reset();
    hsblock_latent_inference(n, Zstat, randK, rng, opt,
                             std::make_tuple(update_data),
                             std::make_tuple(empty));
    Z = Zstat.mean().pruned(ZERO, TOL);

    if (Progress::check_abort()) {
      break;
    }
    prog.increment();

    const Scalar rate = rate0 * std::pow(iter + delay, decay);
    const Scalar score = hsblock_param_inference(
        opt, rate, std::make_tuple(update_data), std::make_tuple(empty));

    if (opt.verbose()) {
      TLOG("Iter [" << iter << "] [" << rate << "] Score = " << score);
    }

    llik(iter) = score;
  }

  return Rcpp::List::create(Rcpp::_["Z"] = Z, Rcpp::_["llik"] = llik);
}

#endif
