#include "L2.h"

// [[Rcpp::depends(RcppArmadillo)]]

// ========================================================
//                   Utility: getCumSumCpp
// ========================================================

arma::mat getCumSumCpp(const arma::mat& X) {

  int nr = X.n_rows;
  int nc = X.n_cols;

  arma::mat cumsumMat(nr + 1, nc, arma::fill::zeros);
  cumsumMat.rows(1, nr) = arma::cumsum(X, 0);

  return cumsumMat;
}

// ========================================================
//                  Utility: Cost_L2 class
// ========================================================

Cost_L2::Cost_L2(const arma::mat& inputMat, bool warnOnce){
  csX = getCumSumCpp(inputMat);
  csXsq = getCumSumCpp(arma::pow(inputMat, 2));
  nr = inputMat.n_rows;
  warnOnce_ = warnOnce;
  keepWarning = not warnOnce;
}

double Cost_L2::eval(int start, int end) const {

  if (start >= end-1) {
    return 0.0;
  }

  int len = end - start;
  return arma::sum(csXsq.row(end) - csXsq.row(start)) - std::pow(arma::norm(csX.row(end) - csX.row(start), 2), 2) / len;
}

void Cost_L2::resetWarning(bool reset){

  warnOnce_ = reset;
  keepWarning = not reset;

}


RCPP_EXPOSED_CLASS(Cost_L2)
  RCPP_MODULE(Cost_L2_module) {
    Rcpp::class_<Cost_L2>("Cost_L2")
    .constructor<arma::mat, bool>()
    .method("eval", &Cost_L2::eval,
    "Evaluate L2 cost on interval (start, end]")
    .method("resetWarning", &Cost_L2::resetWarning,
    "Set the status of warnOnce_");
  }
