// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <random>
#include <string>
#include <vector>
#include <numeric>

#include "convergence.hh"
#include "options.hh"
#include "rcpp_util.hh"
#include "tuple_util.hh"

#ifndef RCPP_HSBLOCK_HH_
#define RCPP_HSBLOCK_HH_

using Mat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
using Vec = Eigen::Matrix<float, Eigen::Dynamic, 1>;
using SpMat = Eigen::SparseMatrix<float, Eigen::ColMajor>;
using Scalar = Mat::Scalar;
using Index = Mat::Index;

#endif
