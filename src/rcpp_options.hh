#include "options.hh"

#ifndef RCPP_OPTIONS_HH_
#define RCPP_OPTIONS_HH_

void set_options_from_list(Rcpp::List& _list, options_t& opt) {
  if (_list.containsElementNamed("tree.depth")) opt.TREE_DEPTH = Rcpp::as<int>(_list["tree.depth"]);
  if (_list.containsElementNamed("rseed"))
    opt.RSEED = Rcpp::as<int>(_list["rseed"]);
  if (_list.containsElementNamed("vbiter"))
    opt.VBITER = Rcpp::as<int>(_list["vbiter"]);
  if (_list.containsElementNamed("verbose"))
    opt.VERBOSE = Rcpp::as<bool>(_list["verbose"]);
  if (_list.containsElementNamed("k")) opt.K = Rcpp::as<int>(_list["k"]);
}

#endif
