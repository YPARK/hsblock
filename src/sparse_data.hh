// sparse data matrix
// (c) Yongjin Park, 2012
#ifndef SPARSE_DATA_HH
#define SPARSE_DATA_HH

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <cmath>
#include <iostream>
#include <numeric>

using namespace boost;
using namespace std;

// dimension index
// strong type checker
struct dim_t {
  explicit dim_t(int _val) : val(_val) {
#ifdef DEBUG
    assert(_val >= 0);
#endif
  }
  int val;
};

bool operator==(const dim_t &lhs, const dim_t &rhs);
bool operator!=(const dim_t &lhs, const dim_t &rhs);

namespace sp_stat {
static const double zero_cutoff = 1e-5;
}

////////////////////////////////////////////////////////////////
// sparse stat vector type
// zero-based indexing
class sp_stat_vec_t {
 public:
  // elementwise arithmetic operations
  // can be done while traversing via
  // iterators
  // e.g., iterate over k-th order stat
  class iterator_t {
   public:
    explicit iterator_t(int _order,
                        const unordered_map<int, double>::const_iterator &_jt)
        : order(_order), jt(_jt) {
#ifdef DEBUG
      assert(_order >= 0);
#endif
    }

    iterator_t &operator=(const iterator_t &sit) {
      jt = sit.jt;
      order = sit.order;
      return *this;
    }

    iterator_t &operator++() {
      ++jt;
      return *this;
    }

    iterator_t operator++(int) {
      iterator_t ret = *this;
      ++(*this);
      return ret;
    }

    iterator_t &operator+(int step) {
#ifdef DEBUG
      assert(step >= 0);
#endif
      std::advance(jt, step);
      return *this;
    }

    bool operator==(const iterator_t &rhs) const { return rhs.jt == jt; }

    bool operator!=(const iterator_t &rhs) const { return rhs.jt != jt; }

    double value() const {
      double val = jt->second;  // raw value
      return pow(val, order);
    }

    double operator*() const { return value(); }

    double dim() const { return jt->first; }

    friend class sp_stat_vec_t;

   private:
    double order;
    unordered_map<int, double>::const_iterator jt;
  };

  ////////////////////////////////////////////////////////////////
  sp_stat_vec_t();
  sp_stat_vec_t(const sp_stat_vec_t &);
  sp_stat_vec_t &operator=(const sp_stat_vec_t &);

  ~sp_stat_vec_t();

  void clear();

  // set and get value
  void operator()(const dim_t &dim, double value);
  void increase(const dim_t &dim, double value);
  double operator()(const dim_t &dim) const;
  bool is_empty(const dim_t &dim) const;

  size_t cardinality() const;
  int get_max_dim() const;
  int get_min_dim() const;
  size_t size() const;

  iterator_t begin_nonzero(int order) const;
  iterator_t begin_nonzero() const;
  iterator_t end_nonzero() const;

  friend class iterator_t;

  ////////////////////////////////////////////////////////////////
  // simple elementwise arithmetic operations
  // += and scalar multiplication
  sp_stat_vec_t &operator+=(const sp_stat_vec_t &);
  sp_stat_vec_t &operator-=(const sp_stat_vec_t &);
  sp_stat_vec_t &operator+=(const double);
  sp_stat_vec_t &operator-=(const double);
  sp_stat_vec_t &operator*=(const double);

  void increase(const sp_stat_vec_t &, double factor);
  void decrease(const sp_stat_vec_t &, double factor);

  // this * (1-eta) + rhs * eta
  sp_stat_vec_t &convex(const sp_stat_vec_t &, double eta);

  void dump() const;

 private:
  void set(int d, double value);
  double get(int d) const;

  void apply_add(const sp_stat_vec_t &, double);
  void apply_sub(const sp_stat_vec_t &, double);
  void apply_add(const sp_stat_vec_t &);
  void apply_sub(const sp_stat_vec_t &);
  void apply_scalar_mult(const double);

  int max_dim;  // maximum dimension index
  int min_dim;  // minimum dimension index
  unordered_map<int, double> data;
};

// inner product of two vectors
// trasverse over fist vector's index
double inner_product(const sp_stat_vec_t &xVec, const sp_stat_vec_t &yVec);

// powered inner-product
double inner_product(const sp_stat_vec_t &xVec, double p,
                     const sp_stat_vec_t &yVec, double q);

// norm
double norm(const sp_stat_vec_t &xVec, int p);

// sum of x^p
double sum(const sp_stat_vec_t &xVec, int p);

// easy computation in log-scale
double log_sum(const sp_stat_vec_t &);

double log_inner_product(const sp_stat_vec_t &xv, const sp_stat_vec_t &yv);

// approximately equal
bool approx_equal(const sp_stat_vec_t &x, const sp_stat_vec_t &y);

////////////////////////////////////////////////////////////////
// sparse stat matrix
// simple inherit smap of all methods of ds.hh
// but with more usability
// * row-major sparse matrix
class sp_stat_mat_t {
 public:
  typedef unordered_map<int, sp_stat_vec_t *> data_map_t;

  sp_stat_mat_t() {}
  sp_stat_mat_t(const sp_stat_mat_t &);
  ~sp_stat_mat_t();

  void operator()(const dim_t &r, const dim_t &c, double val);
  double operator()(const dim_t &r, const dim_t &c) const;

  void add_row(const dim_t &r);
  void add_row(const dim_t &r, const sp_stat_vec_t &row);
  bool has_row(const dim_t &r) const;
  void del_row(const dim_t &r);

  sp_stat_vec_t &pull_row(const dim_t &r);
  const sp_stat_vec_t &pull_const_row(const dim_t &r) const;

  template <typename T>
  class iterator_t {
   public:
    explicit iterator_t(typename T::const_iterator _it) : it(_it) {}

    iterator_t &operator=(const iterator_t &idx) { it = idx.it; }

    iterator_t &operator++() {
      ++it;
      return *this;
    }

    bool operator==(const iterator_t &rhs) const { return it == rhs.it; }

    bool operator!=(const iterator_t &rhs) const { return it != rhs.it; }

    int operator*() { return it->first; }

   private:
    typename T::const_iterator it;
  };

  typedef iterator_t<data_map_t> row_iterator_t;
  typedef sp_stat_vec_t::iterator_t col_iterator_t;

  row_iterator_t begin_row() const { return row_iterator_t(data.cbegin()); }
  row_iterator_t end_row() const { return row_iterator_t(data.cend()); }

  col_iterator_t begin_col(const dim_t &r) const {
    return pull_const_row(r).begin_nonzero();
  }
  col_iterator_t end_col(const dim_t &r) const {
    return pull_const_row(r).end_nonzero();
  }

  void clear();
  size_t size() const;

  void dump() const;

 private:
  data_map_t data;
};

#endif  // EOF
