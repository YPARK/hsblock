#ifndef RCPP_HSBLOCK_IMPL_HH_
#define RCPP_HSBLOCK_IMPL_HH_

////////////////////////////////////////////////////////////////
// A collection of functors needed for inference in btree
template <typename TT>
struct hsb_update_func_t {
  using Tree = TT;
  using Node = typename TT::node_ptr_t;
  using Data = typename TT::node_data_t;
  using Scalar = typename Data::Unit;

  explicit hsb_update_func_t(TT& tt, SpMat& adj, Mat& cc, Mat& zz, Vec& nn)
      : tree(tt),
        A(adj),
        C(cc),
        Z(zz),
        N(nn),
        n(cc.cols()),
        K(cc.rows()),
        Zstat(cc.rows(), cc.cols()),
        delta_score(zz.rows()) {
#ifdef DEBUG
    ASSERT(C.rows() == K && C.cols() == n,    // match dim
           "C [" << K << " x " << n << "]");  //
    ASSERT(Z.rows() == K && Z.cols() == n,    // match dim
           "Z [" << K << " x " << n << "]");  //
    ASSERT(N.rows() == K && N.cols() == 1,    // match dim
           "N [" << K << " x " << 1 << "]");  //
#endif
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
    C = Z * A;
    N = Z * Mat::Ones(n, 1);
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
  SpMat& A;  // n x n adjaency matrix
  Mat& C;    // K x n cluster degree matrix
  Mat& Z;    // K x n latent membership matrix
  Vec& N;    // K x 1 size vector

  const Index n;
  const Index K;

  running_stat_t<Mat> Zstat;

  void reset_z_stat() { Zstat.reset(); }

  void update_z_stat() { Zstat(Z); }

  /////////////////////////////////////////////
  // Note: this doesn't work in delta-update //
  /////////////////////////////////////////////

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
      r->data.delta_stat_dik = C(kk, ii);
      r->data.delta_stat_nik = Z(kk, ii);
      r->data.delta_stat_nk = N(kk);

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

template <typename RandDisc, typename RNG, typename... UpdateFunc>
void hsblock_latent_inference(const Index numVertex, RandDisc& randK, RNG& rng,
                              const options_t& opt,
                              std::tuple<UpdateFunc...>&& update_func_tup) {
  const Index K = randK.K;

  Index vertex_ii = 0;   // vertex index
  Index cluster_kk = 0;  // cluster index

  Vec logScore(K);

  auto init_func = [](auto&& func) {
    func.calibrate_cnz();
    func.reset_z_stat();
  };

  auto accum_eval = [&logScore, &vertex_ii](auto&& func) {
    logScore += func.eval_tree_delta_score(vertex_ii);
  };

  auto discount_matrix = [&vertex_ii](auto&& func) {

    const SpMat& A = func.A;
    Mat& C = func.C;
    Mat& Z = func.Z;
    Vec& N = func.N;

    C -= Z.col(vertex_ii) * A.row(vertex_ii);
    N -= Z.col(vertex_ii);
    Z.col(vertex_ii).setZero();
  };

  auto update_matrix = [&vertex_ii, &cluster_kk](auto&& func) {

    const SpMat& A = func.A;
    Mat& C = func.C;
    Mat& Z = func.Z;
    Vec& N = func.N;

    Z(cluster_kk, vertex_ii) = 1.0;
    C += Z.col(vertex_ii) * A.row(vertex_ii);
    N += Z.col(vertex_ii);
  };

  auto update_stat = [&](auto&& func) { func.update_z_stat(); };

  func_apply(init_func, std::move(update_func_tup));

  Progress prog(opt.inner_iter(), !opt.verbose());

  for (Index inner_iter = 0; inner_iter < opt.inner_iter(); ++inner_iter) {
    if (Progress::check_abort()) {
      break;
    }
    prog.increment();

    for (Index ii = 0; ii < numVertex; ++ii) {
      logScore.setZero();
      vertex_ii = ii;
      func_apply(accum_eval, std::move(update_func_tup));
      func_apply(discount_matrix, std::move(update_func_tup));
      cluster_kk = randK(logScore, rng);
      func_apply(update_matrix, std::move(update_func_tup));
    }

    func_apply(update_stat, std::move(update_func_tup));
  }
}

#endif
