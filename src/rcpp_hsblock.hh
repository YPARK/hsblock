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
#include "mathutil.hh"

#ifndef RCPP_HSBLOCK_HH_
#define RCPP_HSBLOCK_HH_

using Mat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
using Vec = Eigen::Matrix<float, Eigen::Dynamic, 1>;
using SpMat = Eigen::SparseMatrix<float, Eigen::ColMajor>;
using Scalar = Mat::Scalar;
using Index = Mat::Index;

#include "hsb_data.hh"
#include "hsb_func.hh"
#include "hsb_func_dummy.hh"
#include "rcpp_hsblock_inference.hh"

template <typename Tree, typename UpdateData>
Rcpp::List var_em(const SpMat adj, const SpMat latent_init,
                  const options_t& opt);

bool valid_bern_data(const SpMat adj);
bool valid_pois_data(const SpMat adj);

#include "rcpp_hsblock_impl.hh"

#endif
