#include "LinearL1.h"

// [[Rcpp::depends(RcppArmadillo)]]

// ========================================================
//                   Cost_LinearL1 class
// ========================================================

Cost_LinearL1::Cost_LinearL1(const arma::mat& Y,
                             const arma::mat& X,
                             bool intercept_,
                             bool multitask_,
                             bool warnOnce,
                             int max_iter_,
                             double tol_,
                             double eps_) {

  intercept = intercept_;
  multitask = multitask_;
  warnOnce_ = warnOnce;
  keepWarning = !warnOnce;

  max_iter = max_iter_;
  tol = tol_;
  eps = eps_;

  if (Y.n_rows != X.n_rows) {
    Rcpp::stop("Number of observations in response and covariate matrices must match!");
  }

  nr = Y.n_rows;
  pY = Y.n_cols;
  int pX = X.n_cols;

  J = pX + (intercept ? 1 : 0);

  if (nr < J) {
    Rcpp::stop("Not enough observations to fit LAD regression!");
  }

  Yfull = Y;

  Zfull.set_size(nr, J);
  if (intercept) {
    Zfull.col(0).ones();
    Zfull.cols(1, J - 1) = X;
  } else {
    Zfull = X;
  }
}

double Cost_LinearL1::eval(int start, int end) const {

  int n = end - start;
  if (n < J || start >= end - 1) {
    return 0.0;
  }

  arma::mat Z = Zfull.rows(start, end - 1);
  arma::mat Y = Yfull.rows(start, end - 1);

  arma::mat B;
  // Initial OLS
  bool success = arma::solve(B, Z, Y, arma::solve_opts::likely_sympd);
  if (!success) {
    if (warnOnce_) {
      Rcpp::warning("OLS initialisation failed. Using approximate solve.");
      warnOnce_ = false;
    }
    arma::solve(B, Z, Y, arma::solve_opts::force_approx);
  }

  if (multitask) {
    // ---------------------------
    // Multitask IRLS
    // ---------------------------
    arma::mat B_old;
    for (int k = 0; k < max_iter; ++k) {
      B_old = B;
      arma::mat R = Y - Z * B;                           // residuals per row
      arma::vec w = 1.0 / arma::clamp(arma::sqrt(arma::sum(R % R, 1)), eps, arma::datum::inf);

      arma::mat Zw = Z.each_col() % w;
      arma::mat XtWX = Z.t() * Zw;
      arma::mat XtWY = Z.t() * (Y.each_col() % w);

      bool ok = arma::solve(B, XtWX, XtWY, arma::solve_opts::likely_sympd);
      if (!ok) {
        arma::solve(B, XtWX, XtWY, arma::solve_opts::force_approx);
      }

      if (arma::norm(B - B_old, "fro") < tol) break;
    }

    arma::mat R = Y - Z * B;
    double cost = arma::accu(arma::sqrt(arma::sum(R % R, 1)));
    return cost;

  } else {
    // ---------------------------
    // Single-task IRLS (per column)
    // ---------------------------
    arma::mat B_old = B;
    for (int j = 0; j < pY; ++j) {
      arma::vec b_j = B.col(j);
      arma::vec y_j = Y.col(j);
      for (int k = 0; k < max_iter; ++k) {
        arma::vec r = y_j - Z * b_j;
        arma::vec w = 1.0 / arma::clamp(arma::abs(r), eps, arma::datum::inf);
        arma::mat Zw = Z.each_col() % w;
        arma::vec XtwYw = Z.t() * (w % y_j);
        bool ok = arma::solve(b_j, Z.t() * Zw, XtwYw, arma::solve_opts::likely_sympd);
        if (!ok) {
          arma::solve(b_j, Z.t() * Zw, XtwYw, arma::solve_opts::force_approx);
        }
        if (arma::norm(b_j - B_old.col(j), 2) < tol) break;
        B_old.col(j) = b_j;
      }
      B.col(j) = b_j;
    }

    // Sum absolute residuals across tasks
    arma::mat R = Y - Z * B;
    double cost = arma::accu(arma::abs(R));
    return cost;
  }
}

void Cost_LinearL1::resetWarning(bool reset) {
  warnOnce_ = reset;
  keepWarning = !reset;
}
