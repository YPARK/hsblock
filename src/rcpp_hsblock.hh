// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

#include <Eigen/Eigenvalues>
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

#ifndef RCPP_HSBLOCK_HH_
#define RCPP_HSBLOCK_HH_

using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
using Scalar = Mat::Scalar;
using Index = Mat::Index;

template<typename Tree>
struct hsblock_inference_t {

  using Node = typename Tree::Node;

  void merge_left_right(Node left,
			Node right,
			Node parent) {
    // Data& deg2left = left.data.deg;
    // Data& deg2right = right.data.deg;
    // parent.data.deg = deg2left + deg2right;

  }



};


struct tag_bottom {};
struct tag_inter {};

template<typename Tree,
	 typename LeafFun,
	 typename LRFun,
	 typename Opt>
void eval_local_score(typename Tree::node_ptr_t r,
		      LeafFun& Fleaf,
		      LRFun& Flr,
		      Opt& opt)
{			  

  if(r->is_leaf()) {
    LeafFun(r->data, opt);
    return;
  }

  // in-order traversal
  eval_local_score(r->get_left(), Fleaf, Flr, opt);
  eval_local_score(r->get_right(), Fleaf, Flr, opt);

  // opt.Fun(r->data)

  // accumuate
  //
  // opt.Fun(left, right, r)
  //
  // r.data = (r->get_left()).data + (r->get_right()).data;



  return;
}




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
