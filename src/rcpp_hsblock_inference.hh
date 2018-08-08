#include <algorithm>
#include <vector>
#include "hsb_func.hh"
#include "hsb_func_dummy.hh"
#include "hsb_update_data.hh"

#ifndef RCPP_HSBLOCK_INFERENCE_HH_
#define RCPP_HSBLOCK_INFERENCE_HH_

template <typename RandDisc, typename Stat, typename RNG, typename... ModelData,
          typename... NullData>
inline void hsblock_latent_inference(const Index numVertex, Stat& zstat,
                                     RandDisc& randK, RNG& rng,
                                     const options_t& opt, bool is_final,
                                     std::tuple<ModelData...>&& model_data_tup,
                                     std::tuple<NullData...>&& null_data_tup) {
  const Index K = randK.K;
  const Scalar ZERO = 0.0;
  const Scalar TOL = 1e-4;

  Index vertex_ii = 0;   // vertex index
  Index cluster_kk = 0;  // cluster index

  Vec logScore(K);

  auto init_data = [](auto&& data) { data.calibrate(); };

  auto discount_matrix = [&vertex_ii](auto&& data) {
    remove_vertex(data, vertex_ii);
  };

  // update matrix by stochastic sampling of z(k,i)
  auto update_matrix_stoch = [&vertex_ii, &cluster_kk](auto&& data) {
    assign_vertex(data, vertex_ii, cluster_kk);
  };

  auto prune_matrix = [&ZERO, &TOL](auto&& data) { data.C.prune(ZERO, TOL); };

  // evaluate delta likelihood functions

  auto accum_eval = [&vertex_ii, &logScore](auto&& data) {
    logScore += eval_tree_delta_score(data, vertex_ii);
  };

  auto accum_eval_null = [&vertex_ii, &logScore](auto&& data) {
    logScore -= eval_tree_delta_score(data, vertex_ii);
  };

  /////////////////////
  // start inference //
  /////////////////////

  func_apply(init_data, std::move(model_data_tup));
  func_apply(init_data, std::move(null_data_tup));

  const Index max_iter = is_final ? (opt.final_inner_iter() + opt.burnin_iter())
                                  : (opt.inner_iter() + opt.burnin_iter());

  for (Index inner_iter = 0; inner_iter < max_iter; ++inner_iter) {
    for (Index ii = 0; ii < numVertex; ++ii) {
      logScore.setZero();
      vertex_ii = ii;

      func_apply(discount_matrix, std::move(model_data_tup));
      func_apply(discount_matrix, std::move(null_data_tup));

      func_apply(accum_eval, std::move(model_data_tup));
      func_apply(accum_eval_null, std::move(null_data_tup));

      cluster_kk = randK(logScore, rng);
      func_apply(update_matrix_stoch, std::move(model_data_tup));
      func_apply(update_matrix_stoch, std::move(null_data_tup));
    }

    if (inner_iter > opt.burnin_iter() &&
        ((inner_iter % opt.record_interval()) == 0)) {
      func_apply(prune_matrix, std::move(model_data_tup));
      auto& data = std::get<0>(model_data_tup);  // Z is shared
      data.Z.prune(ZERO, TOL);
      zstat.add(data.Z);
    }
  }
}

template <typename... ModelData>
inline Scalar hsblock_empirical_score(
    std::tuple<ModelData...>&& model_data_tup) {
  Scalar score = 0.0;

  auto clear_update = [](auto&& data) {
    clear_tree_data(data);
    data.calibrate();
  };

  auto increase_update = [](auto&& data) { increase_tree_stat(data); };

  auto eval_tree_add = [&score](auto&& data) {
    score += eval_tree_stat_score(data);
  };

  func_apply(clear_update, std::move(model_data_tup));
  func_apply(increase_update, std::move(model_data_tup));
  func_apply(eval_tree_add, std::move(model_data_tup));
  func_apply(clear_update, std::move(model_data_tup));

  return score;
}

template <typename... ModelData>
inline Scalar hsblock_var_score(std::tuple<ModelData...>&& model_data_tup) {
  Scalar score = 0.0;

  auto clear_update = [](auto&& data) {
    clear_tree_data(data);
    data.calibrate();
  };

  auto increase_update = [](auto&& data) { increase_tree_stat(data); };

  auto eval_tree_add = [&score](auto&& data) {
    score += eval_tree_score(data);
  };

  func_apply(clear_update, std::move(model_data_tup));
  func_apply(increase_update, std::move(model_data_tup));
  func_apply(eval_tree_add, std::move(model_data_tup));
  func_apply(clear_update, std::move(model_data_tup));

  return score;
}

template <typename... ModelData, typename... NullData, typename Scalar>
inline Scalar hsblock_param_inference(const Scalar rate,
                                      std::tuple<ModelData...>&& model_data_tup,
                                      std::tuple<NullData...>&& null_data_tup) {
  Scalar score = 0.0;

  auto clear_increase_update = [&rate](auto&& data) {
    clear_tree_data(data);
    data.calibrate();
    increase_tree_stat(data);
    update_tree_var(data, rate);
  };

  auto eval_tree_add = [&score](auto&& data) {
    score += eval_tree_score(data);
  };

  auto eval_tree_subtract = [&score](auto&& data) {
    score -= eval_tree_score(data);
  };

  func_apply(clear_increase_update, std::move(model_data_tup));
  func_apply(clear_increase_update, std::move(null_data_tup));
  func_apply(eval_tree_add, std::move(model_data_tup));
  func_apply(eval_tree_subtract, std::move(null_data_tup));

  return score;
}

  ///////////////
  // deprecaed //
  ///////////////

  /////////////////////////////////////////////////////////////////////
  // This could save some time if the graphs are totally assortative //
  /////////////////////////////////////////////////////////////////////

  // Index kMin, kMax;

  // auto narrow_clusters = [&kMin, &kMax, &vertex_ii](auto&& func) {
  //   for (SpMat::InnerIterator kt(func.C, vertex_ii); kt; ++kt) {
  //     if (kMin > kt.index()) kMin = kt.index();
  //     if (kMax < kt.index()) kMax = kt.index();
  //   }
  // };

  // auto update_stat = [](auto&& func) { func.update_z_stat(); };

  // auto accum_eval_econ = [&logScore, &vertex_ii, &kMin, &kMax](auto&& func) {
  //   auto lca = func.tree.get_lca_node_obj(kMin, kMax);
  //   logScore += func.eval_tree_delta_score(lca, vertex_ii);
  // };

  ///////////////
  // deprecaed //
  ///////////////

  // update matrix normalizing z
  // auto update_matrix_norm = [&](auto&& func) {
  //   const SpMat& At = func.At;
  //   SpMat& C = func.C;
  //   SpMat& Z = func.Z;
  //   Vec& ClustSize = func.ClustSize;
  //   Z.col(vertex_ii) +=
  //       randK.take_prob(logScore).sparseView().pruned(ZERO, TOL);
  //   Z.col(vertex_ii) = Z.col(vertex_ii) / Z.col(vertex_ii).sum();
  //   // this: C += Z.col(vertex_ii) * A.row(vertex_ii);
  //   for (SpMat::InnerIterator jt(At, vertex_ii); jt; ++jt) {
  //     const Index jj = jt.index();
  //     const Scalar Aij = jt.value();
  //     for (SpMat::InnerIterator kt(Z, vertex_ii); kt; ++kt) {
  //       const Index kk = kt.index();
  //       const Scalar Zki = kt.value();
  //       C.coeffRef(kk, jj) += Zki * Aij;
  //     }
  //   }
  //   ClustSize += Z.col(vertex_ii);
  // };

#endif
