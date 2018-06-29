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
#include "options.hh"
#include "rcpp_options.hh"
#include "rcpp_util.hh"
#include "tuple_util.hh"

#include "eigen_sampler.hh"
#include "eigen_util.hh"
#include "hsb_data.hh"
#include "mathutil.hh"

#ifndef RCPP_HSBLOCK_HH_
#define RCPP_HSBLOCK_HH_

using Mat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
using Vec = Eigen::Matrix<float, Eigen::Dynamic, 1>;
using SpMat = Eigen::SparseMatrix<float, Eigen::ColMajor>;
using Scalar = Mat::Scalar;
using Index = Mat::Index;

template <typename TT>
struct hsb_update_func_t;

template <typename RandDisc, typename Stat, typename RNG,
          typename... UpdateFunc>
inline void hsblock_latent_inference(
    const Index numVertex, Stat& zstat, RandDisc& randK, RNG& rng,
    const options_t& opt, std::tuple<UpdateFunc...>&& update_func_tup);

#include "rcpp_hsblock_impl.hh"

#endif
