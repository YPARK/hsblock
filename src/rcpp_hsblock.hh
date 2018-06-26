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

#ifndef RCPP_HSBLOCK_HH_
#define RCPP_HSBLOCK_HH_

using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
using Scalar = Mat::Scalar;
using Index = Mat::Index;

////////////////////////////////////////////////////////////////
// A collection of functors needed for inference in btree
template <typename Tree>
struct hsb_update_func_t {
  using Node = typename Tree::node_ptr_t;
  using Data = typename Tree::node_data_t;
  using Scalar = typename Data::Unit;

  explicit hsb_update_func_t(Mat& cc, Mat& zz)
      : C(cc), Z(zz), N(zz.cols()), n(cc.rows()), K(cc.cols()) {
#ifdef DEBUG
    ASSERT(C.rows() == n && C.cols() == K,    // match dim
           "C [" << n << " x " << K << "]");  //
    ASSERT(Z.rows() == n && Z.cols() == K,    // match dim
           "Z [" << n << " x " << K << "]");  //
#endif
  }

  void init() {

    // 1. clear tree data to zero

    N = Z.transpose() * Mat::Ones(n, 1);

    // 2. collect statistics
    // for each network vertex

    // 3. 

  }


  // Initialize tree model with respect to Z
  void clear_tree_data(Node r) {
    if (r->is_leaf()) {
      init(r->data);
    } else {
      init_tree_data(r->get_left());
      init_tree_data(r->get_right());
      init(r->data);
    }
  }

  // TODO: _collect_tree_stat

  // collect statistics in tree
  void init_tree_stat(Node r, Index ii) {
    if (r->is_leaf()) {
      const Index kk = r->leaf_idx();
      const Scalar z_ik = Z(ii, kk);
      const Scalar d_ik = C(ii, kk);
      const Scalar n_k = N(kk);

      // statistics to push up
      r->data.delta_stat_edge = d_ik;
      r->data.delta_stat_total = z_ik;

      const Scalar half = 0.5;

      // E = sum_ij A_ij z_ik z_jk / 2 + sum_i z_ik e_i
      //     sum_i [ d_ik * z_ik / 2 + e_i * z_ik]
      const Scalar dE = d_ik * half;

      // T = sum_ij z_ik z_jk I[i!=j] / 2
      //   = sum_i [ z_ik * (n_k-z_ik) /2 ]
      const Scalar dT = z_ik * (n_k - z_ik) * half;

      // increase stat edge and stat total
      r->data.stat_edge += dE;
      r->data.stat_total += dT;

      return;
    }
  }



  // collect partial statistics with respect to vertex ii
  void collect_delta_stat(Node r, Index ii) {
    if (r->is_leaf()) {
      const Index k = r->leaf_idx();
      const Scalar d_ik = C(ii, k);
      const Scalar nk = Z(ii, k);

      r->data.delta_stat_edge = d_ik;

      // r->data.delta

      return;
    }
    collect_delta_stat(r->get_left(), ii);
    collect_delta_stat(r->get_right(), ii);

    merge_left_right_delta(r->get_left()->data,   // left
                           r->get_right()->data,  // right
                           r->data);

    dump(r->data);

    return;
  }

  // TODO: remove ii from current cluster

  // TODO: assign ii to a new cluster

  Mat& C;  // n x K cluster degree matrix
  Mat& Z;  // n x K latent membership matrix
  Vec N;   // K x 1 size vector

  const Index n;
  const Index K;

  // const SpMat A;  // n x n adjacency matrix for this snapshot
  // Mat ClustDeg;
  // C = A * Z = [n x n] * [n x K]
};

// score_func(sufficient_stat_data, new_edge, new_tot)

// degree_func(G, ClustDeg)

// tot_func(G, Volume, Size)

// TODO: use functor, template

// template<typename TreeNodePtr, typename Graph>
// struct hsblock_score_t {

//   Scalar calc_degree_sum(TreeNodePtr r, Index curr_i) {
//     Scalar ret = 0.0;
//     if(r -> is_leaf()) { // bottom
//       int k = r->leaf_idx();
//       ret = clust_deg_mat(curr_i, k);
//       memo_deg_sum(dim_t(r->hash()), ret);
//       return ret;
//     }

//     Scalar deg2left = calc_degree_sum(r->get_left(), curr_i);
//     Scalar deg2right = calc_degree_sum(r->get_right(), curr_i);

//     memo_deg2left[r->hash()] = deg2left;
//     memo_deg2right[r->hash()] = deg2right;

//     ret = deg2left + deg2right;
//     return ret;
//   }

//   Scalar calc_tot_sum(TreeNodePtr r) {
//   }

//   template<typename F>
//   Scalar calc_partial_score(TreeNodePtr r, F score_func) {
//     Scalar score = 0.0;
//     if(r->is_leaf()){
//       // dim_t rr(r->hash())
//     // Scalar score = score_func(r->data,
//     // memo_deg_sum(rr),
//     // memo_tot_sum(rr));
//       // memo_left_score
//       return score;
//     }

//     // calc_partial_score(r->get_left()

//     return score;
//   }

//   Index curr_vertex_index;

//   // Mat ;

//   // node_data_t& D = r->data;

//   sp_stat_mat_t clust_deg_mat;    // vertex x cluster

//   // memoorary statistics
//   sp_stat_vec_t memo_deg_sum;     // tree node x 1
//   sp_stat_vec_t memo_tot_sum;     // tree node x 1
//   sp_stat_vec_t memo_left_score;  // tree node x 1
//   sp_stat_vec_t memo_right_score; // tree node x 1

//   const Graph& G;
// };

#endif
