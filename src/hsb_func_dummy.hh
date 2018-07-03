#include "dummy.hh"

#ifndef HSB_FUNC_DUMMY_HH_
#define HSB_FUNC_DUMMY_HH_

struct hsb_func_dummy_t {
  using DATA_TAG = HSB_EMPTY_DATA;

  inline void calibrate() {}

  template <typename T>
  inline const dummy_mat_t& eval_tree_delta_score(T) {
    return dummy_mat;
  }

  template <typename T>
  inline void remove_vertex(T) {}

  template <typename T>
  inline void assign_vertex(T, T) {}

  dummy_mat_t dummy_mat;
  dummy_mat_t Z;
  dummy_mat_t ClustSize;
};

template <typename Index>
inline const dummy_mat_t& _eval_tree_delta_score(hsb_func_dummy_t& dummy,
                                                 const Index ii,
                                                 const HSB_EMPTY_DATA) {
  return dummy.dummy_mat;
}

template <typename T>
inline void _remove_vertex(hsb_func_dummy_t& dummy, T, const HSB_EMPTY_DATA) {}

template <typename T>
inline void _assign_vertex(hsb_func_dummy_t& dummy, T, T,
                           const HSB_EMPTY_DATA) {}

Scalar _eval_tree_score(hsb_func_dummy_t& dummy, const HSB_EMPTY_DATA) {
  return 0.0;
}

Scalar _eval_tree_stat_score(hsb_func_dummy_t& dummy, const HSB_EMPTY_DATA) {
  return 0.0;
}

////////////////////////////////////////////////////////////////

void _update_tree_var(hsb_func_dummy_t& dummy, const Scalar rate,
                      const HSB_EMPTY_DATA) {}

void _clear_tree_data(hsb_func_dummy_t& dummy, const HSB_EMPTY_DATA) {}

////////////////////////////////////////////////////////////////

template <typename Node>
void __increase_tree_stat(hsb_func_dummy_t& dummy, Node r, const Index ii,
                          const HSB_EMPTY_DATA) {}

void _increase_tree_stat(hsb_func_dummy_t& dummy, const HSB_EMPTY_DATA) {}

void _increase_tree_stat(hsb_func_dummy_t& dummy, const Index ii,
                         const HSB_EMPTY_DATA) {}

template <typename Node>
void __decrease_tree_stat(hsb_func_dummy_t& dummy, Node r, const Index ii,
                          const HSB_EMPTY_DATA) {}

void _decrease_tree_stat(hsb_func_dummy_t& dummy, const HSB_EMPTY_DATA) {}

void _decrease_tree_stat(hsb_func_dummy_t& dummy, const Index ii,
                         const HSB_EMPTY_DATA) {}

#endif
