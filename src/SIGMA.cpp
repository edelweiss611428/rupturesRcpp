#include "SIGMA.h"


// [[Rcpp::depends(RcppArmadillo)]]

// covariancePrecomputer constructor
covariancePrecomputer::covariancePrecomputer(const arma::mat& X) {

  nr = X.n_rows;
  nc = X.n_cols;

  S_k.set_size(nr + 1, nc);
  S_k.zeros();

  Q_k.resize(nr + 1, arma::mat(nc, nc, arma::fill::zeros));

  for (int i = 1; i <= nr; i++) {
    Q_k[i] = Q_k[i - 1] + X.row(i - 1).t() * X.row(i - 1);
  }
  S_k.rows(1, nr) = arma::cumsum(X, 0);
}

// covarianceComputer method
arma::mat covariancePrecomputer::covarianceComputer(int start, int end) const {

  int len = end - start;
  if (len <= 1) {
    return arma::mat(nc, nc, arma::fill::zeros);
  }

  arma::rowvec segSum = S_k.row(end) - S_k.row(start);
  arma::mat segOuter = Q_k[end] - Q_k[start];
  arma::vec mean = segSum.t() / len;

  return (segOuter / len) - (mean * mean.t());
}

Cost_SIGMA::Cost_SIGMA(const arma::mat& inputMat,
                       const bool& addSmallDiag, const double& epsilon)

  : preComp(inputMat) {
  addSmallDiag_ = addSmallDiag;
  epsilon_ = epsilon;
  nr = inputMat.n_rows;

}

double Cost_SIGMA::eval(int start, int end) const {

  arma::mat covMat = preComp.covarianceComputer(start, end);
  double sign = 0.0, logDet = 0.0;

  if (start >= end-1) { // If(failed), return 0
    return -std::numeric_limits<double>::max();
  }

  if (addSmallDiag_) {
    covMat.diag() += epsilon_;
  }

  bool success = arma::log_det(logDet, sign, covMat);

  if (success && sign > 0 && std::isfinite(logDet)) {
    return logDet * (end - start);
  }

  return -std::numeric_limits<double>::max(); // If(failed) return the most negative double

}

RCPP_EXPOSED_CLASS(Cost_SIGMA)
  RCPP_MODULE(Cost_SIGMA_module) {
    Rcpp::class_<Cost_SIGMA>("Cost_SIGMA")
    .constructor<arma::mat, bool, double>()
    .method("eval", &Cost_SIGMA::eval,
    "Evaluate SIGMA cost on interval (start, end]")
    ;
  }
