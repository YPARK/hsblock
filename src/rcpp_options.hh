#include "options.hh"

#ifndef RCPP_OPTIONS_HH_
#define RCPP_OPTIONS_HH_

void set_options_from_list(Rcpp::List& _list, options_t& opt) {
  if (_list.containsElementNamed("tree.depth"))
    opt.TREE_DEPTH = Rcpp::as<int>(_list["tree.depth"]);

  if (_list.containsElementNamed("rseed"))
    opt.RSEED = Rcpp::as<int>(_list["rseed"]);

  if (_list.containsElementNamed("vbiter"))
    opt.VBITER = Rcpp::as<int>(_list["vbiter"]);

  if (_list.containsElementNamed("inner.iter"))
    opt.INNER_ITER = Rcpp::as<int>(_list["inner.iter"]);

  if (_list.containsElementNamed("final.inner.iter"))
    opt.FINAL_INNER_ITER = Rcpp::as<int>(_list["final.inner.iter"]);

  if (_list.containsElementNamed("inner_iter"))
    opt.INNER_ITER = Rcpp::as<int>(_list["inner_iter"]);

  if (_list.containsElementNamed("burnin.iter"))
    opt.BURNIN_ITER = Rcpp::as<int>(_list["burnin.iter"]);

  if (_list.containsElementNamed("record.interval"))
    opt.RECORD_INTERV = Rcpp::as<int>(_list["record.interval"]);

  if (_list.containsElementNamed("verbose"))
    opt.VERBOSE = Rcpp::as<bool>(_list["verbose"]);

  if (_list.containsElementNamed("rate"))
    opt.RATE0 = Rcpp::as<float>(_list["rate"]);

  if (_list.containsElementNamed("decay"))
    opt.DECAY = Rcpp::as<float>(_list["decay"]);

  if (_list.containsElementNamed("delay"))
    opt.DELAY = Rcpp::as<float>(_list["delay"]);
}

#endif
