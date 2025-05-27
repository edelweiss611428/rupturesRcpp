#include "Cost.h"

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat getCumsumCpp(const arma::mat& X) {
  int nr = X.n_rows;
  int nc = X.n_cols;

  arma::mat cumsumMat(nr + 1, nc, arma::fill::zeros);
  cumsumMat.rows(1, nr) = arma::cumsum(X, 0);

  return cumsumMat;
}

Cost::Cost(const arma::mat& inputMat) {
  X = inputMat;
  csX = getCumsumCpp(inputMat);
  csXsq = getCumsumCpp(arma::pow(inputMat, 2));
  nr = X.n_rows;
}

double Cost::effEvalCpp(int start, int end) const {
  if (start == end - 1) {
    return 0.0;
  }

  int len = end - start;
  double sumErrXsq = arma::sum(csXsq.row(end) - csXsq.row(start));
  double sqErrsumX = std::pow(arma::norm(csX.row(end) - csX.row(start), 2), 2);

  return sumErrXsq - sqErrsumX / len;
}

