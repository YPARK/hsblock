#ifndef SPARSE_DATA_FUNC_HH_
#define SPARSE_DATA_FUNC_HH_

#include <functional>
#include "rcpp_util.hh"
#include "sparse_data.hh"

void row_for_each(
    const sp_stat_mat_t &G, const dim_t &r,
    std::function<void(const dim_t &, const dim_t &, double)> func) {
  if (G.has_row(r)) {
    const sp_stat_vec_t &row = G.pull_const_row(r);
    for (sp_stat_vec_t::iterator_t it = row.begin_nonzero();
         it != row.end_nonzero(); ++it)
      func(r, dim_t(it.dim()), it.value());
  }
}

////////////////////////////////////////////////////////////////
// visit column neighbors in G
// accumulate vector of d(i->k)
// G and Z must be sparse, otherwise slow
class func_degree_t {
 public:
  explicit func_degree_t(const sp_stat_mat_t &_mu, sp_stat_mat_t &_D)
      : Z(_mu), D(_D) {}

  ~func_degree_t() {}

  // column name and edge weight
  void operator()(const dim_t &r, const dim_t &c, double edge) {
    if (r != c && Z.has_row(c)) {
      if (!D.has_row(r)) D.add_row(r);
      sp_stat_vec_t &ret = D.pull_row(r);
      for (sp_stat_mat_t::col_iterator_t ki = Z.begin_col(c);
           ki != Z.end_col(c); ++ki) {
        int k = ki.dim();
        double z = ki.value();
        ret.increase(dim_t(k), z * edge);
      }
    }
  }

 private:
  const sp_stat_mat_t &Z;
  sp_stat_mat_t &D;
};

#endif
