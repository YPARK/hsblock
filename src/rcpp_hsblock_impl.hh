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

  const Scalar rate0 = opt.rate0();
  const Scalar decay = opt.decay();
  const Scalar delay = opt.delay();

  if (rate0 <= 0 || rate0 > 1) {
    ELOG("rate0 (0, 1] parameter is wrong: " << rate0);
    return Rcpp::List::create();
  }

  if (decay >= 0) {
    ELOG("decay parameter must be negative: " << decay);
    return Rcpp::List::create();
  }

  if (delay < 1) {
    ELOG("delay parameter must be at least 1: " << delay);
    return Rcpp::List::create();
  }

  std::mt19937 rng(opt.rseed());
  Vec llik_v(opt.vbiter() + 1);
  Vec llik_e(opt.vbiter() + 1);

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

  Progress prog(opt.vbiter(), !opt.verbose());
  Scalar escore = hsblock_empirical_score(std::make_tuple(update_data));
  Scalar vscore = hsblock_var_score(std::make_tuple(update_data));

  llik_e(0) = escore;
  llik_v(0) = vscore;

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

    vscore = hsblock_param_inference(rate, std::make_tuple(update_data),
                                     std::make_tuple(empty));

    escore = hsblock_empirical_score(std::make_tuple(update_data));

    if (opt.verbose()) {
      TLOG("Iter [" << iter << "] [" << rate << "] vScore [" << std::setw(10)
                    << vscore << "] eScore [" << std::setw(10) << escore
                    << "]");
    }

    llik_e(iter + 1) = escore;
    llik_v(iter + 1) = vscore;
  }

  return Rcpp::List::create(Rcpp::_["Z"] = Z,
                            Rcpp::_["llik.variational"] = llik_v,
                            Rcpp::_["llik.empirical"] = llik_e);
}

bool valid_bern_data(const SpMat adj) {
  calc_stat_t<Scalar> calc_stat;
  visit(adj, calc_stat);
  const Scalar max_val = calc_stat.max();
  const Scalar min_val = calc_stat.min();
  if (max_val > 1.0) return false;
  if (min_val < 0.0) return false;
  return true;
}

bool valid_pois_data(const SpMat adj) {
  calc_stat_t<Scalar> calc_stat;
  visit(adj, calc_stat);
  const Scalar min_val = calc_stat.min();
  if (min_val < 0.0) return false;
  return true;
}

#endif
