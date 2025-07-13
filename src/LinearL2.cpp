#include "LinearL2.h"

// [[Rcpp::depends(RcppArmadillo)]]

// ========================================================
//                   Cost_LinearL2 class
// ========================================================

Cost_LinearL2::Cost_LinearL2(const arma::mat& Y, const arma::mat& X, bool intercept_, bool warnOnce){

  intercept = intercept_;
  warnOnce_ = warnOnce;
  keepWarning = not warnOnce;

  if (Y.n_rows != X.n_rows) {
    Rcpp::stop("Number of observations in response and covariate matrices must match!");
  }

  nr = Y.n_rows;
  nc = Y.n_cols;  // inherited from CostBase
  pY = nc;
  pX = X.n_cols;

  J = pX + (intercept ? 1 : 0);

  if (nr < J){
    Rcpp::stop("The full dataset contains not enough observations to fit a linear regression model!");
  }

  // Design matrix
  arma::mat Z(nr, J);
  if (intercept) {
    Z.col(0).ones();
    Z.cols(1, J - 1) = X;
  } else {
    Z = X;
  }

  // Allocate cumulative sums
  csXtX.set_size(J, J, nr + 1);
  csXtY.set_size(J, pY, nr + 1);
  csYtY.set_size(pY, pY, nr + 1);

  csXtX.slice(0).zeros();
  csXtY.slice(0).zeros();
  csYtY.slice(0).zeros();

  for (int i = 0; i < nr; ++i) {
    arma::rowvec zi = Z.row(i);
    arma::rowvec yi = Y.row(i);

    csXtX.slice(i + 1) = csXtX.slice(i) + zi.t() * zi;
    csXtY.slice(i + 1) = csXtY.slice(i) + zi.t() * yi;
    csYtY.slice(i + 1) = csYtY.slice(i) + yi.t() * yi;
  }
}

double Cost_LinearL2::eval(int start, int end) const {

  if (start >= end - 1 || end - start < J) {
    return 0.0;
  }

  arma::mat XtX = csXtX.slice(end) - csXtX.slice(start);
  arma::mat XtY = csXtY.slice(end) - csXtY.slice(start);
  arma::mat YtY = csYtY.slice(end) - csYtY.slice(start);

  arma::mat B;
  bool success = arma::solve(B, XtX, XtY, arma::solve_opts::no_approx + arma::solve_opts::likely_sympd);

  if (!success) {
    if (warnOnce_) {
      Rcpp::warning("System is singular. Switching to approximate solve.");
      warnOnce_ = false;
    }
    if (keepWarning) {
      Rcpp::warning("Singular system encountered. Using force_approx.");
    }

    arma::solve(B, XtX, XtY, arma::solve_opts::force_approx);
  }

  double trace_YtY = arma::trace(YtY);
  double trace_BtXtY = arma::trace(B.t() * XtY);

  return std::max(0.0, trace_YtY - trace_BtXtY);
}

void Cost_LinearL2::resetWarning(bool reset) {
  warnOnce_ = reset;
  keepWarning = !reset;
}


// ========================================================
//                 Rcpp Module
// ========================================================

RCPP_EXPOSED_CLASS(Cost_LinearL2)
  RCPP_MODULE(Cost_LinearL2_module) {
    Rcpp::class_<Cost_LinearL2>("Cost_LinearL2")
    .constructor<arma::mat, arma::mat, bool, bool>()
    .method("eval", &Cost_LinearL2::eval,
    "Evaluate linear regression cost on interval (start, end]")
    .method("resetWarning", &Cost_LinearL2::resetWarning,
    "Set the status of warnOnce_");
  }
