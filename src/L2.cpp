#include "L2.h"

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat getCumSumCpp(const arma::mat& X) {

  int nr = X.n_rows;
  int nc = X.n_cols;

  arma::mat cumsumMat(nr + 1, nc, arma::fill::zeros);
  cumsumMat.rows(1, nr) = arma::cumsum(X, 0);

  return cumsumMat;
}

Cost_L2::Cost_L2(const arma::mat& inputMat) {
  csX = getCumSumCpp(inputMat);
  csXsq = getCumSumCpp(arma::pow(inputMat, 2));
  nr = inputMat.n_rows;
}

double Cost_L2::eval(int start, int end) {

  if (start >= end-1) {
    return 0.0;
  }

  int len = end - start;
  return arma::sum(csXsq.row(end) - csXsq.row(start)) - std::pow(arma::norm(csX.row(end) - csX.row(start), 2), 2) / len;
}

RCPP_EXPOSED_CLASS(Cost_L2)
  RCPP_MODULE(Cost_L2_module) {
    Rcpp::class_<Cost_L2>("Cost_L2")
    .constructor<arma::mat>()
    .method("eval", &Cost_L2::eval,
    "Evaluate L2 cost on interval (start, end]")
    ;
  }
