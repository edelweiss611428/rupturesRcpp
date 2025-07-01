#include "Cost_L2.h"

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat getCumsumCpp(const arma::mat& X) {
  int nr = X.n_rows;
  int nc = X.n_cols;

  arma::mat cumsumMat(nr + 1, nc, arma::fill::zeros);
  cumsumMat.rows(1, nr) = arma::cumsum(X, 0);

  return cumsumMat;
}

Cost_L2::Cost_L2(const arma::mat& inputMat) {
  csX = getCumsumCpp(inputMat);
  csXsq = getCumsumCpp(arma::pow(inputMat, 2));
  nr = inputMat.n_rows;
}

double Cost_L2::effEvalCpp(int start, int end,
                           bool addSmallDiag,
                           double epsilon) const {
  if (start == end - 1) {
    return 0.0;
  }

  int len = end - start;
  double sumErrXsq = arma::sum(csXsq.row(end) - csXsq.row(start));
  double sqErrsumX = std::pow(arma::norm(csX.row(end) - csX.row(start), 2), 2);

  return sumErrXsq - sqErrsumX / len;
}


//Return for testing purposes
RCPP_MODULE(Cost_L2_module) {
  class_<Cost_L2>("Cost_L2")

  .constructor<arma::mat>()

  .method("effEvalCpp", &Cost_L2::effEvalCpp,
  "Evaluate L2 cost on interval (start, end]");
}

