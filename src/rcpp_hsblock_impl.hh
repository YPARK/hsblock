#ifndef RCPP_HSBLOCK_IMPL_HH_
#define RCPP_HSBLOCK_IMPL_HH_

#include <algorithm>
#include <vector>

////////////////////////////////////////////////////////////////
// A collection of functors needed for inference in btree
template <typename TT>
struct hsb_update_func_t {
  using Tree = TT;
  using Node = typename TT::node_ptr_t;
  using Data = typename TT::node_data_t;
  using Scalar = typename Data::Unit;

  explicit hsb_update_func_t(TT& tt, SpMat& adj, SpMat& zz)
      : tree(tt),
        n(zz.cols()),
        K(zz.rows()),
        A(adj),
        At(A.transpose()),
        Z(zz),
        C(K, n),
        ClustSize(K),
        delta_score(K) {
    calibrate_cnz();
  }

  // Initialize tree model with respect to Z
  void clear_tree_data() { clear_tree_data(tree.root_node_obj()); }

  void clear_tree_data(Node r) {
    if (r->is_leaf()) {
      clear(r->data);
    } else {
      clear_tree_data(r->get_left());
      clear_tree_data(r->get_right());
      clear(r->data);
    }
  }

  void calibrate_cnz() {
    C = Z * At;
    ClustSize = Z * Mat::Ones(n, 1);
  }

  void init_tree_var() { init_tree_var(tree.root_node_obj()); }

  void init_tree_var(Node r) {
    recursive_tree_delta_stat(r, 0);
    recursive_init_tree_var(r);
  }

  void update_tree_var(Scalar rate) {
    update_tree_var(tree.root_node_obj(), rate);
  }

  void update_tree_var(Node r, Scalar rate) {
    if (r->is_leaf()) {
      update(r->data, rate);
    } else {
      update_tree_var(r->get_left(), rate);
      update_tree_var(r->get_right(), rate);
      update(r->data, rate);
    }
  }

  const Vec& eval_tree_delta_score(Index ii) {
    return eval_tree_delta_score(tree.root_node_obj(), ii);
  }

  const Vec& eval_tree_delta_score(Node r, Index ii) {
    Scalar accum = 0.0;
    recursive_tree_delta_stat(r, ii);
    recursive_eval_delta_score(r, ii);
    recursive_accum_delta_score(r, ii, accum);
    return delta_score;
  }

  Scalar eval_tree_score() { return eval_tree_score(tree.root_node_obj()); }

  Scalar eval_tree_score(Node r) {
    if (r->is_leaf()) {
      eval_score(r->data);
      return r->data.score;
    } else {
      Scalar sc_left = eval_tree_score(r->get_left());
      Scalar sc_right = eval_tree_score(r->get_right());
      eval_score(r->data);
      return sc_left + sc_right;
    }
  }

  void dump_tree_data(Node r) {
    if (r->is_leaf()) {
      dump(r->data);
    } else {
      dump_tree_data(r->get_left());
      dump_tree_data(r->get_right());
      dump(r->data);
    }
  }

  void dump_tree_data() { dump_tree_data(tree.root_node_obj()); }

  TT& tree;  // binary tree object
  const Index n;
  const Index K;
  const SpMat& A;  // n x n adjaency matrix
  const SpMat At;  // n x n adjaency transpose matrix
  SpMat& Z;        // K x n latent membership matrix
  SpMat C;         // K x n cluster degree matrix
  Vec ClustSize;   // K x 1 size vector

  /////////////////////////////////////////////
  // Note: this doesn't work in delta-update //
  /////////////////////////////////////////////

  void increase_tree_stat() {
    for (Index ii = 0; ii < n; ++ii)
      increase_tree_stat(tree.root_node_obj(), ii);
  }

  void increase_tree_stat(Index ii) {
    increase_tree_stat(tree.root_node_obj(), ii);
  }

  void increase_tree_stat(Node r, Index ii) {
    auto inc_op = [](auto prev, auto delt) { return prev + delt; };
    recursive_tree_delta_stat(r, ii);
    recursive_tree_stat(r, ii, inc_op);
  }

  /////////////////////////////////////////////
  // Note: this doesn't work in delta-update //
  /////////////////////////////////////////////

  void decrease_tree_stat(Index ii) {
    decrease_tree_stat(tree.root_node_obj(), ii);
  }

  void decrease_tree_stat(Node r, Index ii) {
    auto dec_op = [](auto prev, auto delt) { return prev - delt; };
    recursive_tree_delta_stat(r, ii);
    recursive_tree_stat(r, ii, dec_op);
  }

 private:
  Vec delta_score;  // K x 1 delta_score vector

  // initialize
  Scalar recursive_init_tree_var(Node r) {
    if (r->is_leaf()) {
      const Scalar level = 0.0;
      init_param(r->data, level);
      return level;
    } else {
      Node left = r->get_left();
      Node right = r->get_right();
      const Scalar one_val = 1.0;
      const Scalar lv_left = recursive_init_tree_var(left);
      const Scalar lv_right = recursive_init_tree_var(right);
      const Scalar level = std::max(lv_left, lv_right) + one_val;
      init_param_intern(r->data, left->data, right->data, level);
      return level;
    }
  }

  // accumulate delta scores
  inline void recursive_accum_delta_score(Node r, Index ii, Scalar accum) {
    if (r->is_leaf()) {
      const Index kk = r->leaf_idx();
      const Scalar accum_tot = accum + r->data.delta_score;
      delta_score(kk) = accum_tot;
    } else {
      Scalar accum_left = accum + r->data.delta_score_left;
      Scalar accum_right = accum + r->data.delta_score_right;
      Scalar denom = std::max(accum_left, accum_right);
      accum_left -= denom;
      accum_right -= denom;
      recursive_accum_delta_score(r->get_left(), ii, accum_left);
      recursive_accum_delta_score(r->get_right(), ii, accum_right);
    }
  }

  // evaluate delta score given delta statistics
  inline void recursive_eval_delta_score(Node r, Index ii) {
    if (r->is_leaf()) {
      eval_delta_score(r->data);
#ifdef DEBUG
      // dump(r->data);
#endif
    } else {
      Node left = r->get_left();
      Node right = r->get_right();

      recursive_eval_delta_score(left, ii);
      recursive_eval_delta_score(right, ii);

      // Choosing the left means modeling on the right
      // using the delta_stat = (diR, nR)
      eval_delta_left(r->data, right->data);

      // Choosing the right means modeling on the left
      // using the delta_stat = (diL, nL)
      eval_delta_right(r->data, left->data);

      eval_delta_score(r->data);

#ifdef DEBUG
      // dump(r->data);
#endif
    }
  }

  // collect partial statistics wrt ii
  inline void recursive_tree_delta_stat(Node r, Index ii) {
    if (r->is_leaf()) {
      const Index kk = r->leaf_idx();

      // statistics to push up
      r->data.delta_stat_dik = C.coeff(kk, ii);
      r->data.delta_stat_nik = Z.coeff(kk, ii);
      r->data.delta_stat_nk = ClustSize(kk);

    } else {
      Node left = r->get_left();
      Node right = r->get_right();
      recursive_tree_delta_stat(left, ii);
      recursive_tree_delta_stat(right, ii);
      merge_left_right_delta(left->data, right->data, r->data);
    }
  }

  /////////////////////////////////////////////
  // Note: this doesn't work in delta-update //
  /////////////////////////////////////////////

  // collect statistics in tree given delta statistics
  template <typename Op>
  inline void recursive_tree_stat(Node r, Index ii, Op op) {
    const Scalar half = 0.5;

    if (r->is_leaf()) {
      const Index kk = r->leaf_idx();

      const Scalar d_ik = r->data.delta_stat_dik;
      const Scalar z_ik = r->data.delta_stat_nik;
      const Scalar n_k = r->data.delta_stat_nk;

      // E = sum_ij A_ij z_ik z_jk / 2 + sum_i z_ik e_i
      //     sum_i [ d_ik * z_ik / 2 + e_i * z_ik]
      const Scalar dE = d_ik * z_ik * half;

      // T = sum_ij z_ik z_jk I[i!=j] / 2
      //   = sum_i [ z_ik * (n_k-z_ik) /2 ]
      const Scalar dT = z_ik * (n_k - z_ik) * half;

      // increase stat edge and stat total
      r->data.stat_edge = op(r->data.stat_edge, dE);
      r->data.stat_total = op(r->data.stat_total, dT);

    } else {
      Node left = r->get_left();
      Node right = r->get_right();

      recursive_tree_stat(left, ii, op);
      recursive_tree_stat(right, ii, op);

      const Scalar niL = left->data.delta_stat_nik;
      const Scalar niR = right->data.delta_stat_nik;

      const Scalar diL = left->data.delta_stat_dik;
      const Scalar diR = right->data.delta_stat_dik;

      const Scalar nL = left->data.delta_stat_nk;
      const Scalar nR = right->data.delta_stat_nk;

      // E = sum_ij A_ij sum_{l in L} sum_{r in R} z_il z_jr
      //     sum_i sum_{l in L} z_il sum_{r in R} sum_j A_ij z_jr
      //     sum_i ( n_iL * d_iR )
      const Scalar dE = niL * diR;  //  + niR * diL;

      // T = sum_{i!=j} sum_{l in L} sum_{r in R} z_il z_jr
      //   = sum_i sum_{l in L} z_il sum_{r in R} sum_{j!=i} z_jr
      //   = sum_i n_iL sum_{r in R} [ n_r - z_ir ]
      //   = sum_i [ n_iL (n_R - n_iR) ]
      const Scalar dT = niL * (nR - niR);  //  + niR * (nL - niL);

      // increase stat edge and stat total
      r->data.stat_edge = op(r->data.stat_edge, dE);
      r->data.stat_total = op(r->data.stat_total, dT);
      // dump(r->data);
    }
  }
};

template <typename RandDisc, typename Stat, typename RNG,
          typename... UpdateFunc>
inline void hsblock_latent_inference(
    const Index numVertex, Stat& zstat, RandDisc& randK, RNG& rng,
    const options_t& opt, std::tuple<UpdateFunc...>&& update_func_tup) {
  const Index K = randK.K;

  Index vertex_ii = 0;   // vertex index
  Index cluster_kk = 0;  // cluster index

  Vec logScore(K);

  auto init_func = [](auto&& func) { func.calibrate_cnz(); };

  auto accum_eval = [&logScore, &vertex_ii](auto&& func) {
    logScore += func.eval_tree_delta_score(vertex_ii);
  };

  auto discount_matrix = [&vertex_ii](auto&& func) {

    const SpMat& At = func.At;
    SpMat& C = func.C;
    SpMat& Z = func.Z;
    Vec& ClustSize = func.ClustSize;

    // this: C -= Z.col(vertex_ii) * A.row(vertex_ii);
    for (SpMat::InnerIterator jt(At, vertex_ii); jt; ++jt) {
      const Index jj = jt.index();
      const Scalar Aij = jt.value();
      for (SpMat::InnerIterator kt(Z, vertex_ii); kt; ++kt) {
        const Index kk = kt.index();
        const Scalar Zki = kt.value();
        C.coeffRef(kk, jj) -= Zki * Aij;
      }
    }
    ClustSize -= Z.col(vertex_ii);
    Z.col(vertex_ii) = Z.col(vertex_ii) * 0.0;
  };

  auto update_matrix = [&vertex_ii, &cluster_kk](auto&& func) {

    const SpMat& At = func.At;
    SpMat& C = func.C;
    SpMat& Z = func.Z;
    Vec& ClustSize = func.ClustSize;

    Z.coeffRef(cluster_kk, vertex_ii) = 1.0;
    // this: C += Z.col(vertex_ii) * A.row(vertex_ii);
    for (SpMat::InnerIterator jt(At, vertex_ii); jt; ++jt) {
      const Index jj = jt.index();
      const Scalar Aij = jt.value();
      C.coeffRef(cluster_kk, jj) += Aij;
    }
    ClustSize += Z.col(vertex_ii);
  };

  auto prune_matrix = [](auto&& func) {
    const Scalar ZERO = 0.0;
    const Scalar TOL = 1e-4;
    func.C.prune(ZERO, TOL);
    func.Z.prune(ZERO, TOL);
  };

  /////////////////////////////////////////////////////////////////////
  // This could save some time if the graphs are totally assortative //
  /////////////////////////////////////////////////////////////////////

  Index kMin, kMax;

  auto narrow_clusters = [&kMin, &kMax, &vertex_ii](auto&& func) {
    for (SpMat::InnerIterator kt(func.C, vertex_ii); kt; ++kt) {
      if (kMin > kt.index()) kMin = kt.index();
      if (kMax < kt.index()) kMax = kt.index();
    }
  };

  auto accum_eval_econ = [&logScore, &vertex_ii, &kMin, &kMax](auto&& func) {
    auto lca = func.tree.get_lca_node_obj(kMin, kMax);
    logScore += func.eval_tree_delta_score(lca, vertex_ii);
  };

  auto update_stat = [](auto&& func) { func.update_z_stat(); };

  func_apply(init_func, std::move(update_func_tup));

  Progress prog(opt.inner_iter() + opt.burnin_iter(), false);

  const Index max_iter = opt.inner_iter() + opt.burnin_iter();

  for (Index inner_iter = 0; inner_iter < max_iter; ++inner_iter) {
    if (Progress::check_abort()) {
      break;
    }
    prog.increment();

    for (Index ii = 0; ii < numVertex; ++ii) {
      logScore.setZero();
      vertex_ii = ii;

      if (opt.economy()) {
        kMin = K;
        kMax = -1;
        func_apply(narrow_clusters, std::move(update_func_tup));
        if (kMax >= 0 && kMin != kMax) {
          func_apply(accum_eval_econ, std::move(update_func_tup));
          func_apply(discount_matrix, std::move(update_func_tup));
          cluster_kk = randK(logScore, kMin, kMax + 1, rng);
          func_apply(update_matrix, std::move(update_func_tup));
        }
      } else {
        func_apply(accum_eval, std::move(update_func_tup));
        func_apply(discount_matrix, std::move(update_func_tup));
        cluster_kk = randK(logScore, rng);
        func_apply(update_matrix, std::move(update_func_tup));
      }
    }

    if (inner_iter > opt.burnin_iter() &&
        ((inner_iter % opt.record_interval()) == 0)) {
      func_apply(prune_matrix, std::move(update_func_tup));
      auto& func = std::get<0>(update_func_tup);  // Z is shared
      zstat.add(func.Z);
    }
  }
}

template <typename... UpdateFunc, typename Scalar>
inline Scalar hsblock_param_inference(
    const options_t& opt, const Scalar rate,
    std::tuple<UpdateFunc...>&& update_func_tup) {
  Scalar score = 0.0;

  auto clear_increase_update = [&rate](auto&& func) {
    func.clear_tree_data();
    func.increase_tree_stat();
    func.update_tree_var(rate);
  };

  auto eval_tree = [&score](auto&& func) { score += func.eval_tree_score(); };

  func_apply(clear_increase_update, std::move(update_func_tup));
  func_apply(eval_tree, std::move(update_func_tup));

  return score;
}

#endif
