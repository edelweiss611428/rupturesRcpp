// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "ogkSIGMA.h"

using namespace Rcpp;
using namespace arma;


// Column-wise MAD (scaled by 1.4826)
inline arma::vec MADcpp(const arma::mat& X) {
  arma::rowvec med = arma::median(X);
  arma::rowvec mad = arma::median(abs(X.each_row() - med));
  return 1.4826 * mad.t();
}

// GK estimator for covariance between two vectors
inline double biVarGK(const arma::vec& x, const arma::vec& y) {
  arma::mat tmpMat(x.n_elem, 2);
  tmpMat.col(0) = x + y;
  tmpMat.col(1) = x - y;

  arma::vec MADxy = MADcpp(tmpMat);
  return (pow(MADxy(0), 2) - pow(MADxy(1), 2)) / 4.0;
}

// [[Rcpp::export]]
arma::mat covMatGK(const arma::mat& X) {
  int p = X.n_cols;
  arma::mat covMat(p, p, arma::fill::zeros);

  for (int i = 0; i < p - 1; ++i) {
    for (int j = i + 1; j < p; ++j) {
      covMat(i, j) = biVarGK(X.col(i), X.col(j));
    }
  }

  covMat = arma::symmatu(covMat);  // fill lower triangle
  arma::vec madVec = MADcpp(X);
  covMat.diag() += madVec % madVec;

  return covMat;
}

// Qn scale estimator (univariate)
double scaleQn(const arma::vec& x) {
  int n = x.n_elem;
  int k = R::choose(std::floor(n / 2.0) + 1, 2);  // kth order statistic

  arma::uvec uppInd = arma::trimatu_ind(arma::size(n, n), 1);
  arma::mat temp = x * arma::ones(1, n);
  arma::vec dist = temp(uppInd);
  temp = temp.t();
  dist -= temp(uppInd);
  dist = arma::sort(abs(dist));

  return 2.2219 * dist(k);  // magic number from Rousseeuw & Croux (1993)
}

arma::mat covOGK(const arma::mat& data) {
  int p = data.n_cols;
  arma::mat res(p, p, arma::fill::zeros);
  arma::vec D(p);
  arma::mat Z = data;

  for (int i = 0; i < p; ++i) {
    D(i) = scaleQn(data.col(i));
  }

  Z.each_row() /= D.t();
  arma::mat U = covMatGK(Z);
  U.diag().ones();  // force correlation matrix

  arma::vec Lambda;
  arma::mat E;
  eig_sym(Lambda, E, U);

  arma::mat A = E;
  A.each_col() %= D;
  Z = Z * E;

  arma::vec Gamma(p);
  for (int i = 0; i < p; ++i) {
    double g = scaleQn(Z.col(i));
    Gamma(i) = g * g;
  }

  res = A * diagmat(Gamma) * A.t();
  return res;
}



// ------------------------------------------------------------
//                     Cost_ogkSIGMA Class
// ------------------------------------------------------------

// Constructor
Cost_ogkSIGMA::Cost_ogkSIGMA(const arma::mat& inputMat, bool warnOnce)
  : X(inputMat), warnOnce_(warnOnce), keepWarning(!warnOnce) {
  nr = inputMat.n_rows;
  nc = inputMat.n_cols;
}

// Evaluate cost on the segment [start, end)
double Cost_ogkSIGMA::eval(int start, int end) const {

  if (start >= end - 1) {
    return 0.0;
  }

  int len = end - start;
  arma::mat submat = X.rows(start, end - 1);  // X[start:(end-1), ]
  arma::mat cov = covOGK(submat);             // robust OGK covariance
  double val = std::log(arma::det(cov)) * len;      // log-determinant cost

  return val;

}

// Reset warning flags
void Cost_ogkSIGMA::resetWarning(bool reset) {
  warnOnce_ = reset;
  keepWarning = !reset;
}


RCPP_EXPOSED_CLASS(Cost_ogkSIGMA)
  RCPP_MODULE(Cost_ogkSIGMA_module) {
    Rcpp::class_<Cost_ogkSIGMA>("Cost_ogkSIGMA")
    .constructor<arma::mat, bool>()
    .method("eval", &Cost_ogkSIGMA::eval,
    "Evaluate ogkSIGMA cost on interval (start, end]")
    .method("resetWarning", &Cost_ogkSIGMA::resetWarning,
    "Set the status of warnOnce_");
  }
