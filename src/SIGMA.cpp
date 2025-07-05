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
                       const bool& addSmallDiag, const double& epsilon,
                       const bool&warnOnce)

  : preComp(inputMat) {
  addSmallDiag_ = addSmallDiag;
  epsilon_ = epsilon;
  nr = inputMat.n_rows;
  nc = inputMat.n_cols;
  lbDet = nc*log(epsilon);
  warnOnce_ = warnOnce;
  const bool keepWarning = not warnOnce;
  (void)keepWarning;

}

double Cost_SIGMA::eval(int start, int end) {

  arma::mat covMat = preComp.covarianceComputer(start, end);
  double logDet = 0.0;

  if (addSmallDiag_) {
    covMat.diag() += epsilon_;
  }

  logDet = arma::log_det_sympd(covMat);

  if (std::isfinite(logDet)) {
    if(not addSmallDiag_){
      return logDet * (end - start);

    } else{
      if(logDet < lbDet){
        return lbDet * (end - start);

      } else{
        return logDet * (end - start);

      }
    }
  } else if(addSmallDiag_ and epsilon_ > 0.0){

    if(warnOnce_){
      warning("`covMat` is singular! Consider increasing either `epsilon` or `minSize`");
      warnOnce_ = false;
      warning("Return the lower-bound `p*log(epsilon)*segLen`!");
    }

    if(keepWarning){
      warning("`covMat` is singular! Consider increasing either `epsilon` or `minSize`");
      warning("Return the lower-bound `p*log(epsilon)*segLen`!");
    }

    return lbDet*(end-start);

  } else{
    stop("`covMat` is singular! Consider using `addSmallDiag` option!");
  }

}

RCPP_EXPOSED_CLASS(Cost_SIGMA)
  RCPP_MODULE(Cost_SIGMA_module) {
    Rcpp::class_<Cost_SIGMA>("Cost_SIGMA")
    .constructor<arma::mat, bool, double, bool>()
    .method("eval", &Cost_SIGMA::eval,
    "Evaluate SIGMA cost on interval (start, end]")
    ;
  }
