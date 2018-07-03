#ifndef HSB_UDPATE_DATA_HH_
#define HSB_UDPATE_DATA_HH_

// Tag for dummy or non-dummy
struct HSB_EMPTY_DATA {};
struct HSB_NETWORK_DATA {};

// Tag for degree correction
struct NON_DEGREE_CORRECTED {};
struct DEGREE_CORRECTED {};

// A collection of temporary data for the inference in btree
template <typename TT, typename _DC_TAG>
struct hsb_update_data_t {
  using DC_TAG = _DC_TAG;
  using DATA_TAG = HSB_NETWORK_DATA;
  using Tree = TT;
  using Node = typename TT::node_ptr_t;
  using Data = typename TT::node_data_t;
  using Scalar = typename Data::Unit;

  explicit hsb_update_data_t(TT& tt, SpMat& adj, SpMat& zz)
      : tree(tt),
        n(zz.cols()),
        K(zz.rows()),
        A(adj),
        At(A.transpose()),
        Z(zz),
        C(K, n),
        ClustSize(K),
        Volume(K),
        Deg(n),
        delta_score(K) {
    calibrate();
  }

  // Initialize tree model with respect to Z
  void calibrate() {
    At = A.transpose();
    Deg = A * Mat::Ones(n, 1);
    DegSum = Deg.sum();
    C = Z * At;                       // sum_j A(i,j) Z(k,j)
    ClustSize = Z * Mat::Ones(n, 1);  // sum_i Z(k,i)
    Volume = Z * Deg.asDiagonal();    // sum_i d(i) Z(k,i)
  }

  TT& tree;  // binary tree object
  const Index n;
  const Index K;
  const SpMat& A;   // n x n adjaency matrix
  SpMat At;         // n x n adjaency transpose matrix
  SpMat& Z;         // K x n latent membership matrix
  SpMat C;          // K x n cluster degree matrix
  Vec ClustSize;    // K x 1 size vector
  Vec Volume;       // K x 1 volume vector
  Vec Deg;          // n x 1 degree vector
  Scalar DegSum;    // total sum of degree sequence
  Vec delta_score;  // K x 1 delta_score vector
};

// factory function

// make hsb

// make dhsb

#endif
