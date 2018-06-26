#ifndef RCPP_HSBLOCK_IMPL_HH_
#define RCPP_HSBLOCK_IMPL_HH_

template <typename T>
void merge_left_right_delta(const hsb_data_t<T>& left,
                            const hsb_data_t<T>& right, hsb_data_t<T>& ret) {
  ret.delta_stat_edge = left.delta_stat_edge + right.delta_stat_edge;
  ret.delta_stat_total = left.delta_stat_total + right.delta_stat_total;
}

template <typename T>
void merge_left_right_stat(const hsb_data_t<T>& left,
                           const hsb_data_t<T>& right, hsb_data_t<T>& ret) {
  ret.stat_edge = left.stat_edge + right.stat_edge;
  ret.stat_total = left.stat_total + right.stat_total;
}


#endif
