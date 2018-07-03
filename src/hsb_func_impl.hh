#ifndef HSB_FUNC_IMPL_HH_
#define HSB_FUNC_IMPL_HH_

////////////////////////////////////////////////////////////////
// delta update of udata
template <typename Udata>
void _remove_vertex(Udata& udata, const Index vertex_ii,
                    const HSB_NETWORK_DATA) {
  // this: C -= Z.col(vertex_ii) * A.row(vertex_ii);
  const SpMat& At = udata.At;
  SpMat& C = udata.C;
  SpMat& Z = udata.Z;
  Vec& ClustSize = udata.ClustSize;
  Vec& Volume = udata.Volume;
  const Vec& Deg = udata.Deg;

  for (SpMat::InnerIterator jt(At, vertex_ii); jt; ++jt) {
    const Index jj = jt.index();
    if (vertex_ii == jj) continue;
    const Scalar Aij = jt.value();
    for (SpMat::InnerIterator kt(Z, vertex_ii); kt; ++kt) {
      const Index kk = kt.index();
      const Scalar Zki = kt.value();
      C.coeffRef(kk, jj) -= Zki * Aij;
    }
  }
  ClustSize -= Z.col(vertex_ii);
  Volume -= Deg(vertex_ii) * Z.col(vertex_ii);
  Z.col(vertex_ii) = Z.col(vertex_ii) * 0.0;
}

template <typename Udata>
void remove_vertex(Udata& udata, const Index vertex_ii) {
  _remove_vertex(udata, vertex_ii, typename Udata::DATA_TAG());
}

template <typename Udata>
void _assign_vertex(Udata& udata, const Index vertex_ii, const Index cluster_kk,
                    const HSB_NETWORK_DATA) {
  const SpMat& At = udata.At;
  SpMat& C = udata.C;
  SpMat& Z = udata.Z;
  Vec& ClustSize = udata.ClustSize;
  Vec& Volume = udata.Volume;
  const Vec& Deg = udata.Deg;

  Z.coeffRef(cluster_kk, vertex_ii) = 1.0;
  // this: C += Z.col(vertex_ii) * A.row(vertex_ii);
  for (SpMat::InnerIterator jt(At, vertex_ii); jt; ++jt) {
    const Index jj = jt.index();
    if (vertex_ii == jj) continue;
    const Scalar Aij = jt.value();
    C.coeffRef(cluster_kk, jj) += Aij;
  }
  ClustSize += Z.col(vertex_ii);
  Volume += Deg(vertex_ii) * Z.col(vertex_ii);
}

template <typename Udata>
void assign_vertex(Udata& udata, const Index vertex_ii,
                   const Index cluster_kk) {
  _assign_vertex(udata, vertex_ii, cluster_kk, typename Udata::DATA_TAG());
}

template <typename Udata>
void _clear_tree_data(Udata& udata, const HSB_NETWORK_DATA) {
  recursive_clear_tree_data(udata.tree.root_node_obj());
}

template <typename Udata>
void clear_tree_data(Udata& udata) {
  _clear_tree_data(udata, typename Udata::DATA_TAG());
}

template <typename Node>
void recursive_clear_tree_data(Node r) {
  if (!r->is_leaf()) {
    recursive_clear_tree_data(r->get_left());
    recursive_clear_tree_data(r->get_right());
  }
  clear(r->data);
}

template <typename Udata, typename Node>
void _init_tree_var(Udata& udata, Node r, const HSB_NETWORK_DATA) {
  recursive_tree_delta_stat(udata, r, 0, typename Udata::DC_TAG());
  recursive_init_tree_var(r);
}

template <typename Udata, typename Node>
void init_tree_var(Udata& udata, Node r) {
  _init_tree_var(udata, r, typename Udata::DATA_TAG());
}

template <typename Udata>
void _init_tree_var(Udata& udata, const HSB_NETWORK_DATA) {
  _init_tree_var(udata, udata.tree.root_node_obj(), typename Udata::DATA_TAG());
}

template <typename Udata>
void init_tree_var(Udata& udata) {
  _init_tree_var(udata, typename Udata::DATA_TAG());
}

template <typename Node>
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

////////////////////////////////////////////////////////////////
// collect partial statistics wrt ii
template <typename Udata, typename Node>
inline void recursive_tree_delta_stat(Udata& udata, Node r, const Index ii,
                                      const NON_DEGREE_CORRECTED) {
  SpMat& C = udata.C;
  SpMat& Z = udata.Z;
  Vec& ClustSize = udata.ClustSize;
  Vec& Volume = udata.Volume;

  if (r->is_leaf()) {
    const Index kk = r->leaf_idx();

    // statistics to push up
    r->data.delta_stat_dik = C.coeff(kk, ii);
    r->data.delta_stat_nik = Z.coeff(kk, ii);
    r->data.delta_stat_nk = ClustSize(kk);

  } else {
    Node left = r->get_left();
    Node right = r->get_right();
    recursive_tree_delta_stat(udata, left, ii, NON_DEGREE_CORRECTED());
    recursive_tree_delta_stat(udata, right, ii, NON_DEGREE_CORRECTED());
    merge_left_right_delta(left->data, right->data, r->data);
  }
}

//*******************************//
//* Degree-corrected statistics *//
//*******************************//

// collect partial statistics wrt ii with degree correction
template <typename Udata, typename Node>
inline void recursive_tree_delta_stat(Udata& udata, Node r, const Index ii,
                                      const DEGREE_CORRECTED) {
  SpMat& C = udata.C;
  SpMat& Z = udata.Z;
  Vec& ClustSize = udata.ClustSize;
  Vec& Volume = udata.Volume;
  const Scalar DegSum = udata.DegSum;
  const Vec& Deg = udata.Deg;

  if (r->is_leaf()) {
    const Index kk = r->leaf_idx();

    // statistics to push up
    r->data.delta_stat_dik = C.coeff(kk, ii);
    r->data.delta_stat_nik = Z.coeff(kk, ii);
    r->data.delta_stat_nk = Deg(ii) * Volume(kk) / DegSum;

  } else {
    Node left = r->get_left();
    Node right = r->get_right();
    recursive_tree_delta_stat(udata, left, ii, DEGREE_CORRECTED());
    recursive_tree_delta_stat(udata, right, ii, DEGREE_CORRECTED());
    merge_left_right_delta(left->data, right->data, r->data);
  }
}

// evaluate delta score given delta statistics
template <typename Node>
inline void recursive_eval_delta_score(Node r, const Index ii) {
  if (r->is_leaf()) {
    eval_delta_score(r->data);
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

    // Unnecessary to compute the following
    // eval_delta_score(r->data);
  }
}

// accumulate delta scores
template <typename Udata, typename Node>
inline void recursive_accum_delta_score(Udata& udata, Node r, const Index ii,
                                        const Scalar accum) {
  if (r->is_leaf()) {
    const Index kk = r->leaf_idx();
    udata.delta_score(kk) = accum + r->data.delta_score;
  } else {
    // NOTE: local adjustment can make results totally different
    Scalar accum_left = accum + r->data.delta_score_left;
    Scalar accum_right = accum + r->data.delta_score_right;
    recursive_accum_delta_score(udata, r->get_left(), ii, accum_left);
    recursive_accum_delta_score(udata, r->get_right(), ii, accum_right);
  }
}

template <typename Udata, typename Node>
const auto& _eval_tree_delta_score(Udata& udata, Node r, const Index ii,
                                   const HSB_NETWORK_DATA) {
  auto accum = 0.0;

  recursive_tree_delta_stat(udata, r, ii, typename Udata::DC_TAG());
  recursive_eval_delta_score(r, ii);
  recursive_accum_delta_score(udata, r, ii, accum);

  return udata.delta_score;
}

template <typename Udata, typename Node>
const auto& eval_tree_delta_score(Udata& udata, Node r, const Index ii) {
  return _eval_tree_delta_score(udata, r, ii, typename Udata::DATA_TAG());
}

template <typename Udata>
const auto& _eval_tree_delta_score(Udata& udata, const Index ii,
                                   const HSB_NETWORK_DATA) {
  return _eval_tree_delta_score(udata, udata.tree.root_node_obj(), ii,
                                typename Udata::DATA_TAG());
}

template <typename Udata>
const auto& eval_tree_delta_score(Udata& udata, const Index ii) {
  return _eval_tree_delta_score(udata, ii, typename Udata::DATA_TAG());
}

////////////////////////////////////////////////////////////////
// Evaluation of scores
template <typename Node>
Scalar recursive_eval_tree_score(Node r) {
  if (r->is_leaf()) {
    eval_score(r->data);
    return r->data.score;
  } else {
    Scalar sc_left = recursive_eval_tree_score(r->get_left());
    Scalar sc_right = recursive_eval_tree_score(r->get_right());
    eval_score(r->data);
    return sc_left + sc_right + r->data.score;
  }
}

template <typename Udata>
Scalar _eval_tree_score(Udata& udata, const HSB_NETWORK_DATA) {
  return recursive_eval_tree_score(udata.tree.root_node_obj());
}

template <typename Udata>
Scalar eval_tree_score(Udata& udata) {
  return _eval_tree_score(udata, typename Udata::DATA_TAG());
}

// Evaluation of stat scores
template <typename Node>
Scalar recursive_eval_tree_stat_score(Node r) {
  if (r->is_leaf()) {
    eval_stat_score(r->data);
    return r->data.score;
  } else {
    Scalar sc_left = recursive_eval_tree_stat_score(r->get_left());
    Scalar sc_right = recursive_eval_tree_stat_score(r->get_right());
    eval_stat_score(r->data);
    return sc_left + sc_right + r->data.score;
  }
}

template <typename Udata>
Scalar _eval_tree_stat_score(Udata& udata, const HSB_NETWORK_DATA) {
  return recursive_eval_tree_stat_score(udata.tree.root_node_obj());
}

template <typename Udata>
Scalar eval_tree_stat_score(Udata& udata) {
  return _eval_tree_stat_score(udata, typename Udata::DATA_TAG());
}

////////////////////////////////////////////////////////////////
// dump tree data
template <typename Node>
void dump_tree_data(Node r) {
  if (r->is_leaf()) {
    dump(r->data);
  } else {
    dump_tree_data(r->get_left());
    dump_tree_data(r->get_right());
    dump(r->data);
  }
}

template <typename Udata>
void dump_tree_data_udata(Udata& udata) {
  dump_tree_data(udata.tree.root_node_obj());
}

////////////////////////////////////////////////////////////////
// update variational parameters
template <typename Node, typename Scalar>
void recursive_update_tree_var(Node r, const Scalar rate) {
  if (r->is_leaf()) {
    update(r->data, rate);
  } else {
    recursive_update_tree_var(r->get_left(), rate);
    recursive_update_tree_var(r->get_right(), rate);
    update(r->data, rate);
  }
}

template <typename Udata, typename Scalar>
void _update_tree_var(Udata& udata, const Scalar rate, const HSB_NETWORK_DATA) {
  recursive_update_tree_var(udata.tree.root_node_obj(), rate);
}

template <typename Udata, typename Scalar>
void update_tree_var(Udata& udata, const Scalar rate) {
  _update_tree_var(udata, rate, typename Udata::DATA_TAG());
}

/////////////////////////////////////////////
// Note: this doesn't work in delta-update //
/////////////////////////////////////////////

// collect statistics in tree given delta statistics
template <typename Udata, typename Node, typename Op>
inline void recursive_tree_stat(Udata& udata, Node r, const Index ii, Op op,
                                const NON_DEGREE_CORRECTED) {
  const Scalar half = 0.5;

  SpMat& C = udata.C;
  SpMat& Z = udata.Z;
  Vec& ClustSize = udata.ClustSize;
  Vec& Volume = udata.Volume;

  if (r->is_leaf()) {
    const Index kk = r->leaf_idx();

    r->data.delta_stat_dik = C.coeff(kk, ii);
    r->data.delta_stat_nik = Z.coeff(kk, ii);
    r->data.delta_stat_nk = ClustSize(kk);

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

    recursive_tree_stat(udata, left, ii, op);
    recursive_tree_stat(udata, right, ii, op);
    merge_left_right_delta(left->data, right->data, r->data);

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

/////////////////////////////////////////////
// Note: this doesn't work in delta-update //
/////////////////////////////////////////////

// collect statistics with degree correction
template <typename Udata, typename Node, typename Op>
inline void recursive_tree_stat(Udata& udata, Node r, const Index ii, Op op,
                                const DEGREE_CORRECTED) {
  SpMat& C = udata.C;
  SpMat& Z = udata.Z;
  Vec& ClustSize = udata.ClustSize;
  Vec& Volume = udata.Volume;
  const Scalar DegSum = udata.DegSum;
  const Vec& Deg = udata.Deg;
  const Scalar half = 0.5;
  const Scalar d_i = Deg(ii);

  if (r->is_leaf()) {
    const Index kk = r->leaf_idx();

    r->data.delta_stat_dik = C.coeff(kk, ii);
    r->data.delta_stat_nik = Z.coeff(kk, ii);
    r->data.delta_stat_nk = Volume(kk);

    const Scalar d_ik = r->data.delta_stat_dik;
    const Scalar z_ik = r->data.delta_stat_nik;
    const Scalar vk = r->data.delta_stat_nk;

    // E = sum_ij A_ij z_ik z_jk / 2 + sum_i z_ik e_i
    //     sum_i [ d_ik * z_ik / 2 + e_i * z_ik]
    const Scalar dE = d_ik * z_ik * half;

    // T = sum_ij I[i!=j] d_i z_ik d_j z_jk / 2 / 2E
    //   = 1/2 * sum_i d_i z_ik sum_{j!=i} d_j z_jk / 2E
    //   = 1/2 * sum_i d_i z_ik (vol(k) - d_i*z_ik) / 2E
    //   = sum_i [ z_ik * d_i * (vol(k) - d_i * z_ik) /2/2E ]
    const Scalar dT = (z_ik * d_i * vk - d_i * z_ik) * half / DegSum;

    // increase stat edge and stat total
    r->data.stat_edge = op(r->data.stat_edge, dE);
    r->data.stat_total = op(r->data.stat_total, dT);

  } else {
    Node left = r->get_left();
    Node right = r->get_right();

    recursive_tree_stat(udata, left, ii, op, DEGREE_CORRECTED());
    recursive_tree_stat(udata, right, ii, op, DEGREE_CORRECTED());
    merge_left_right_delta(left->data, right->data, r->data);

    const Scalar niL = left->data.delta_stat_nik;
    const Scalar niR = right->data.delta_stat_nik;

    const Scalar diL = left->data.delta_stat_dik;
    const Scalar diR = right->data.delta_stat_dik;

    const Scalar vL = left->data.delta_stat_nk;
    const Scalar vR = right->data.delta_stat_nk;

    // E = sum_ij A_ij sum_{l in L} sum_{r in R} z_il z_jr
    //     sum_i sum_{l in L} z_il sum_{r in R} sum_j A_ij z_jr
    //     sum_i ( n_iL * d_iR )
    const Scalar dE = niL * diR;  //  + niR * diL;

    // T = sum_ij sum_{l in L} sum_{r in R} d_i z_il d_j z_jr /2E
    //   = sum_i d_i n_iL * sum_{r in R} [ vol_r - d_i z_ir ] /2E
    //   = sum_i [ d_i n_iL * (vol_R - d_i n_iR) /2E ]
    const Scalar dT = d_i * niL * (vR - d_i * niR) / DegSum;

    // increase stat edge and stat total
    r->data.stat_edge = op(r->data.stat_edge, dE);
    r->data.stat_total = op(r->data.stat_total, dT);
    // dump(r->data);
  }
}

////////////////////////////////////////////////////////////////
// Dispatch depending on degree correction tag

template <typename Udata, typename Node>
void __increase_tree_stat(Udata& udata, Node r, const Index ii,
                          const HSB_NETWORK_DATA) {
  auto inc_op = [](auto prev, auto delt) { return prev + delt; };
  recursive_tree_stat(udata, r, ii, inc_op, typename Udata::DC_TAG());
}

template <typename Udata>
void _increase_tree_stat(Udata& udata, const Index ii, const HSB_NETWORK_DATA) {
  __increase_tree_stat(udata, udata.tree.root_node_obj(), ii,
                       HSB_NETWORK_DATA());
}

template <typename Udata>
void _increase_tree_stat(Udata& udata, const HSB_NETWORK_DATA) {
  for (Index ii = 0; ii < udata.n; ++ii)
    __increase_tree_stat(udata, udata.tree.root_node_obj(), ii,
                         HSB_NETWORK_DATA());
}

template <typename Udata, typename Node>
void increase_tree_stat(Udata& udata, Node r, const Index ii) {
  __increase_tree_stat(udata, r, ii, typename Udata::DATA_TAG());
}

template <typename Udata>
void increase_tree_stat(Udata& udata) {
  _increase_tree_stat(udata, typename Udata::DATA_TAG());
}

template <typename Udata>
void increase_tree_stat(Udata& udata, const Index ii) {
  _increase_tree_stat(udata, ii, typename Udata::DATA_TAG());
}

////////////////////////////////////////////////////////////////

template <typename Udata, typename Node>
void __decrease_tree_stat(Udata& udata, Node r, const Index ii,
                          const HSB_NETWORK_DATA) {
  auto dec_op = [](auto prev, auto delt) { return prev + delt; };
  recursive_tree_stat(udata, r, ii, dec_op, typename Udata::DC_TAG());
}

template <typename Udata>
void _decrease_tree_stat(Udata& udata, const Index ii, const HSB_NETWORK_DATA) {
  __decrease_tree_stat(udata, udata.tree.root_node_obj(), ii,
                       HSB_NETWORK_DATA());
}

template <typename Udata>
void _decrease_tree_stat(Udata& udata, const HSB_NETWORK_DATA) {
  for (Index ii = 0; ii < udata.n; ++ii)
    __decrease_tree_stat(udata, udata.tree.root_node_obj(), ii,
                         HSB_NETWORK_DATA());
}

template <typename Udata, typename Node>
void decrease_tree_stat(Udata& udata, Node r, const Index ii) {
  __decrease_tree_stat(udata, r, ii, typename Udata::DATA_TAG());
}

template <typename Udata>
void decrease_tree_stat(Udata& udata) {
  _decrease_tree_stat(udata, typename Udata::DATA_TAG());
}

template <typename Udata>
void decrease_tree_stat(Udata& udata, const Index ii) {
  _decrease_tree_stat(udata, ii, typename Udata::DATA_TAG());
}

#endif
