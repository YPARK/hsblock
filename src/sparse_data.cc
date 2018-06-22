#include "sparse_data.hh"

bool operator==(const dim_t& lhs, const dim_t& rhs) {
  return lhs.val == rhs.val;
}

bool operator!=(const dim_t& lhs, const dim_t& rhs) {
  return lhs.val != rhs.val;
}

////////////////////////////////////////////////////////////////
// sparse stat vector type
sp_stat_vec_t::sp_stat_vec_t() {
  max_dim = -1;
  min_dim = -1;
}

sp_stat_vec_t::sp_stat_vec_t(const sp_stat_vec_t& rhs)
    : data(rhs.data), max_dim(rhs.max_dim), min_dim(rhs.min_dim) {}

sp_stat_vec_t& sp_stat_vec_t::operator=(const sp_stat_vec_t& rhs) {
  if (this == &rhs) return *this;
  clear();
  data = rhs.data;
  max_dim = rhs.max_dim;
  min_dim = rhs.min_dim;
  return *this;
}

sp_stat_vec_t::~sp_stat_vec_t() {
  // nothing
}

void sp_stat_vec_t::clear() {
  max_dim = -1;
  min_dim = -1;
  data.clear();
}

void sp_stat_vec_t::operator()(const dim_t& dim, double value) {
  int d = dim.val;
  set(d, value);
}

void sp_stat_vec_t::increase(const dim_t& dim, double value) {
  int d = dim.val;
  double newVal = value + get(d);
  set(d, newVal);
}

void sp_stat_vec_t::set(int d, double value) {
  if (abs(value) < sp_stat::zero_cutoff)  // ignore effectively zero
  {
    if (data.count(d) > 0)  // this will remove element
    {
      data.erase(d);  // erase this spot

      if (data.size() > 0) {
        // if this happened to be min or max dim
        // then, assign new min and max dim
        // ~ amortized O(1) cost
        if (d == max_dim || d == min_dim) {
          iterator_t it = begin_nonzero();
          min_dim = it.dim();
          max_dim = it.dim();
          for (; it != end_nonzero(); ++it) {
            int _d = it.dim();
            if (_d > max_dim) {
              max_dim = _d;
            }
            if (_d < min_dim) {
              min_dim = _d;
            }
          }
        }
      } else {
        clear();
      }
    }
    // otherwise, just ignore
  } else {
    max_dim = (d > max_dim) ? d : max_dim;
    if (min_dim < 0) {
      min_dim = d;
    } else if (min_dim > d) {
      min_dim = d;
    }
    data[d] = value;
  }
}

double sp_stat_vec_t::operator()(const dim_t& dim) const {
  return get(dim.val);
}

double sp_stat_vec_t::get(int d) const {
  return (data.count(d) > 0) ? data.at(d) : 0.;
}

sp_stat_vec_t::iterator_t sp_stat_vec_t::begin_nonzero(int order) const {
  return iterator_t(order, data.cbegin());
}

sp_stat_vec_t::iterator_t sp_stat_vec_t::begin_nonzero() const {
  return iterator_t(1, data.cbegin());
}

sp_stat_vec_t::iterator_t sp_stat_vec_t::end_nonzero() const {
  return iterator_t(1, data.cend());
}

size_t sp_stat_vec_t::cardinality() const { return (max_dim + 1); }

int sp_stat_vec_t::get_max_dim() const {
#ifdef DEBUG
  assert(max_dim >= 0);
#endif
  return max_dim;
}

int sp_stat_vec_t::get_min_dim() const {
#ifdef DEBUG
  assert(min_dim >= 0);
#endif
  return min_dim;
}

size_t sp_stat_vec_t::size() const { return data.size(); }

////////////////////////////////////////////////////////////////
// simple elementwise arithmetic operation
// +=
sp_stat_vec_t& sp_stat_vec_t::operator+=(const sp_stat_vec_t& rhs_vec) {
  apply_add(rhs_vec);
  return *this;
}

// -=
sp_stat_vec_t& sp_stat_vec_t::operator-=(const sp_stat_vec_t& rhs_vec) {
  apply_sub(rhs_vec);
  return *this;
}

sp_stat_vec_t& sp_stat_vec_t::operator+=(const double val) {
  for (iterator_t it = begin_nonzero(); it != end_nonzero();) {
    int d = it.dim();
    double v = it.value();
    ++it;
    set(d, v + val);
  }
  return *this;
}

sp_stat_vec_t& sp_stat_vec_t::operator-=(const double val) {
  for (iterator_t it = begin_nonzero(); it != end_nonzero();) {
    int d = it.dim();
    double v = it.value();
    ++it;
    set(d, v - val);
  }
  return *this;
}

void sp_stat_vec_t::apply_add(const sp_stat_vec_t& rhs_vec) {
  apply_add(rhs_vec, 1.);
}

void sp_stat_vec_t::apply_add(const sp_stat_vec_t& rhs_vec, double factor) {
  for (iterator_t rhs = rhs_vec.begin_nonzero(); rhs != rhs_vec.end_nonzero();
       ++rhs) {
    int d = rhs.dim();
    double curr = get(d);
    set(d, curr + factor * rhs.value());
  }
}

void sp_stat_vec_t::apply_sub(const sp_stat_vec_t& rhs_vec) {
  apply_sub(rhs_vec, 1.);
}

void sp_stat_vec_t::apply_sub(const sp_stat_vec_t& rhs_vec, double factor) {
  apply_add(rhs_vec, -factor);
}

void sp_stat_vec_t::increase(const sp_stat_vec_t& rhs, double factor) {
  apply_add(rhs, factor);
}

void sp_stat_vec_t::decrease(const sp_stat_vec_t& rhs, double factor) {
  apply_sub(rhs, factor);
}

// scalar multiplication
sp_stat_vec_t& sp_stat_vec_t::operator*=(const double c) {
  apply_scalar_mult(c);
  return *this;
}

void sp_stat_vec_t::apply_scalar_mult(const double c) {
  for (iterator_t it = begin_nonzero(); it != end_nonzero();) {
    iterator_t itx = it;
    ++it;
    int d = itx.dim();
    double val = itx.value();
    set(d, val * c);
  }
}

// convex combination
sp_stat_vec_t& sp_stat_vec_t::convex(const sp_stat_vec_t& rhs_vec, double eta) {
#ifdef DEBUG
  assert(eta >= 0. && eta <= 1.);
#endif

  if (eta < sp_stat::zero_cutoff) return *this;

  apply_scalar_mult(1. - eta);

  for (iterator_t rhs = rhs_vec.begin_nonzero(); rhs != rhs_vec.end_nonzero();
       ++rhs) {
    int d = rhs.dim();
    double val = rhs.value();
    double curr = get(d);
    set(d, curr + val * eta);
  }
  return *this;
}

bool sp_stat_vec_t::is_empty(const dim_t& dim) const {
  return (data.count(dim.val) == 0);
}

void sp_stat_vec_t::dump() const {
  cerr << "|V| " << cardinality();
  cerr << " #{V} " << size();
  for (iterator_t it = begin_nonzero(); it != end_nonzero(); ++it) {
    cerr << "\t[" << it.dim() << "]";
    cerr << " " << it.value();
  }
  cerr << endl;
}

// inner product
double inner_product(const sp_stat_vec_t& xVec, const sp_stat_vec_t& yVec) {
  double ret = 0.;
  for (sp_stat_vec_t::iterator_t x = xVec.begin_nonzero();
       x != xVec.end_nonzero(); ++x)
    ret += x.value() * yVec(dim_t(x.dim()));
  return ret;
}

double inner_product(const sp_stat_vec_t& xVec, double p,
                     const sp_stat_vec_t& yVec, double q) {
  double ret = 0.;
  double x, y;
  for (sp_stat_vec_t::iterator_t xi = xVec.begin_nonzero(p);
       xi != xVec.end_nonzero(); ++xi) {
    x = xi.value();
    y = yVec(dim_t(xi.dim()));
    ret += x * pow(y, q);
  }
  return ret;
}

// norm
double myAbs(double accum, double x) { return accum + abs(x); }

double norm(const sp_stat_vec_t& xVec, int p) {
#ifdef DEBUG
  assert(p >= 0.);
#endif
  if (p == 0) return sum(xVec, 0);
  if (p == 1)
    return accumulate(xVec.begin_nonzero(), xVec.end_nonzero(), 0., myAbs);
  return pow(sum(xVec, p), 1. / ((double)p));
}

// sum of x^p
double sum(const sp_stat_vec_t& xVec, int p) {
  return accumulate(xVec.begin_nonzero(p), xVec.end_nonzero(), 0.);
}

// log-scale computation
double _log_sum(double log_a, double log_b) {
  double v;

  if (log_a < log_b) {
    v = log_b + log(1 + exp(log_a - log_b));
  } else {
    v = log_a + log(1 + exp(log_b - log_a));
  }
  return (v);
}

double log_sum(const sp_stat_vec_t& xv) {
  double ret = 0.;
  int i = 0;
  for (sp_stat_vec_t::iterator_t xi = xv.begin_nonzero();
       xi != xv.end_nonzero(); ++xi) {
    if (i == 0) {
      ret = xi.value();
    } else {
      ret = _log_sum(ret, xi.value());
    }
    ++i;
  }
  return ret;
}

double log_inner_product(const sp_stat_vec_t& xv, const sp_stat_vec_t& yv) {
  double ret = 0.;
  int i = 0;
  for (sp_stat_vec_t::iterator_t xi = xv.begin_nonzero();
       xi != xv.end_nonzero(); ++xi) {
    double log_xy = xi.value() + yv(dim_t(xi.dim()));
    if (i == 0) {
      ret = log_xy;
    } else {
      ret = _log_sum(ret, log_xy);
    }
    ++i;
  }
  return ret;
}

// approximately equal
bool approx_equal(const sp_stat_vec_t& x, const sp_stat_vec_t& y) {
  // different cardinality or size
  if ((x.cardinality() != y.cardinality()) || (x.size() != y.size()))
    return false;

  for (sp_stat_vec_t::iterator_t xi = x.begin_nonzero(); xi != x.end_nonzero();
       ++xi) {
    dim_t ii(xi.dim());
    if (y(ii) < sp_stat::zero_cutoff) return false;
    if (abs(xi.value() - y(ii)) > sp_stat::zero_cutoff) return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////
// sparse data matrix

// copy constructor
sp_stat_mat_t::sp_stat_mat_t(const sp_stat_mat_t& smat) {
  for (sp_stat_mat_t::row_iterator_t ri = smat.begin_row();
       ri != smat.end_row(); ++ri)
    add_row(dim_t(*ri), smat.pull_const_row(dim_t(*ri)));
}

sp_stat_mat_t::~sp_stat_mat_t() {
  for (row_iterator_t it = begin_row(); it != end_row();) {
    int r = *it;
    ++it;
    del_row(dim_t(r));
  }
}

void sp_stat_mat_t::add_row(const dim_t& r) {
  if (!has_row(r)) data[r.val] = new sp_stat_vec_t;
}

void sp_stat_mat_t::add_row(const dim_t& r, const sp_stat_vec_t& row) {
#ifdef DEBUG
  assert(!has_row(r));  // prevent overwriting
#endif
  data[r.val] = new sp_stat_vec_t(row);
}

bool sp_stat_mat_t::has_row(const dim_t& r) const {
  return (data.count(r.val) > 0);
}

void sp_stat_mat_t::del_row(const dim_t& r) {
  if (has_row(r)) {
    delete data[r.val];
    data.erase(r.val);
  }
}

sp_stat_vec_t& sp_stat_mat_t::pull_row(const dim_t& r) {
#ifdef DEBUG
  assert(has_row(r));
#endif
  return *(data.at(r.val));
}

const sp_stat_vec_t& sp_stat_mat_t::pull_const_row(const dim_t& r) const {
#ifdef DEBUG
  assert(has_row(r));
#endif
  return *(data.at(r.val));
}

void sp_stat_mat_t::operator()(const dim_t& r, const dim_t& c, double val) {
  if (!has_row(r)) add_row(r);
  sp_stat_vec_t& row = pull_row(r);
  row(c, val);
}

double sp_stat_mat_t::operator()(const dim_t& r, const dim_t& c) const {
#ifdef DEBUG
  assert(has_row(r));
#endif
  const sp_stat_vec_t& row = *(data.at(r.val));
  return row(c);
}

// pair<sp_stat_mat_t::col_iterator_t, sp_stat_mat_t::col_iterator_t>
// sp_stat_mat_t::column_iterator(const dim_t& r)
// {
//   sp_stat_vec_t& row = pull_row(r);
//   pair<col_iterator_t,col_iterator_t> ret(row.begin_nonzero(),
//   row.end_nonzero()); return ret;
// }

void sp_stat_mat_t::clear() {
  for (row_iterator_t it = begin_row(); it != end_row(); ++it)
    pull_row(dim_t(*it)).clear();
}

size_t sp_stat_mat_t::size() const { return data.size(); }

void sp_stat_mat_t::dump() const {
  for (row_iterator_t it = begin_row(); it != end_row(); ++it) {
    cerr << *it << " : ";
    pull_const_row(dim_t(*it)).dump();
  }
}
