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

#include "hsb_func.hh"
#include "hsb_func_dummy.hh"
#include "rcpp_hsblock_inference.hh"

// Basic stochastic block model estimation




#endif
