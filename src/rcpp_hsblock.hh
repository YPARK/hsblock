// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "btree.hh"
#include "cnetwork.hh"
#include "options.hh"
#include "rcpp_util.hh"
#include "sparse_data.hh"
#include "sparse_data_func.hh"
#include "tuple_util.hh"

#include "hsb_data.hh"
#include "mathutil.hh"
#include "eigen_sampler.hh"
#include "eigen_util.hh"

#ifndef RCPP_HSBLOCK_HH_
#define RCPP_HSBLOCK_HH_

using Mat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
using Vec = Eigen::Matrix<float, Eigen::Dynamic, 1>;
using SpMat = Eigen::SparseMatrix<float, Eigen::ColMajor>;
using Scalar = Mat::Scalar;
using Index = Mat::Index;

////////////////////////////////////////////////////////////////
// A collection of functors needed for inference in btree
template <typename Tree>
struct hsb_update_func_t {
  using Node = typename Tree::node_ptr_t;
  using Data = typename Tree::node_data_t;
  using Scalar = typename Data::Unit;

  explicit hsb_update_func_t(Mat& cc, Mat& zz, Vec& nn)
      : C(cc),
        Z(zz),
        N(nn),
        delta_score(zz.rows()),
        n(cc.cols()),
        K(cc.rows()) {
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
  void clear_tree_data(Node r) {
    if (r->is_leaf()) {
      clear(r->data);
    } else {
      clear_tree_data(r->get_left());
      clear_tree_data(r->get_right());
      clear(r->data);
    }
  }

  void increase_tree_stat(Node r, Index ii) {
    auto inc_op = [](auto prev, auto delt) { return prev + delt; };
    recursive_tree_delta_stat(r, ii);
    recursive_tree_stat(r, ii, inc_op);
  }

  void decrease_tree_stat(Node r, Index ii) {
    auto dec_op = [](auto prev, auto delt) { return prev - delt; };
    recursive_tree_delta_stat(r, ii);
    recursive_tree_stat(r, ii, dec_op);
  }

  const Vec& eval_tree_delta_score(Node r, Index ii) {
    Scalar accum = 0.0;
    recursive_tree_delta_stat(r, ii);
    recursive_eval_delta_score(r, ii);
    recursive_accum_delta_score(r, ii, accum);
    return delta_score;
  }

  Mat& C;  // K x n cluster degree matrix
  Mat& Z;  // K x n latent membership matrix
  Vec& N;  // K x 1 size vector

 private:
  Vec delta_score;  // K x 1 delta_score vector

  const Index n;
  const Index K;

 private:
  // accumulate delta scores
  inline void recursive_accum_delta_score(Node r, Index ii, Scalar accum) {
    if (r->is_leaf()) {
      const Index kk = r->leaf_idx();
      const Scalar accum_tot = accum + r->data.delta_score;
      delta_score(kk) = accum_tot;
    } else {
      Scalar accum_left = accum + r->data.delta_score_left;
      Scalar accum_right = accum + r->data.delta_score_right;

      // standardize
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
      auto left = r->get_left();
      auto right = r->get_right();

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
      auto left = r->get_left();
      auto right = r->get_right();
      recursive_tree_delta_stat(left, ii);
      recursive_tree_delta_stat(right, ii);
      merge_left_right_delta(left->data, right->data, r->data);
    }
  }

  // collect statistics in tree given delta statistics
  template <typename Op>
  inline void recursive_tree_stat(Node r, Index ii, Op op) {
    if (r->is_leaf()) {
      const Index kk = r->leaf_idx();

      const Scalar d_ik = r->data.delta_stat_dik;
      const Scalar z_ik = r->data.delta_stat_nik;
      const Scalar n_k = r->data.delta_stat_nk;

      const Scalar half = 0.5;

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
      auto left = r->get_left();
      auto right = r->get_right();

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
      const Scalar dE = niL * diR;

      // T = sum_{i!=j} sum_{l in L} sum_{r in R} z_il z_jr
      //   = sum_i sum_{l in L} z_il sum_{r in R} sum_{j!=i} z_jr
      //   = sum_i n_iL sum_{r in R} [ n_r - z_ir ]
      //   = sum_i [ n_iL (n_R - n_iR) ]
      const Scalar dT = niL * (nR - niR);

      // increase stat edge and stat total
      r->data.stat_edge = op(r->data.stat_edge, dE);
      r->data.stat_total = op(r->data.stat_total, dT);
      // dump(r->data);
    }
  }
};

#endif
