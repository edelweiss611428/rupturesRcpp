#include "costs.h"

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

arma::mat getCumSumCpp(const arma::mat& X) {
  arma::mat cumsumMat(X.n_rows + 1, X.n_cols, arma::fill::zeros);
  cumsumMat.rows(1, X.n_rows) = arma::cumsum(X, 0);
  return cumsumMat;
}

Rcpp::NumericVector asNumeric(const arma::rowvec& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}

template <class T>
void resetWarningR(T* obj, bool reset) {
  obj->resetWarning(reset);
}

} // namespace

// ========================================================
//                        CostBase
// ========================================================

void CostBase::checkSegment(int start, int end) const {
  if (start >= end) {
    Rcpp::stop("`start < end` must be true!");
  }
  if (start >= nr || start < 0) {
    Rcpp::stop("`0 <= start < nSamples` must be true!");
  }
  if (end > nr || end <= 0) {
    Rcpp::stop("`0 < end <= nSamples` must be true!");
  }
}

int CostBase::warnLevel() const {
  if (warnOnce_) {
    warnOnce_ = false;
    return 1;
  }
  return keepWarning ? 2 : 0;
}

double CostBase::logDetCost(const arma::mat& covMat, bool addSmallDiag, double epsilon,
                             double lbDet, int len) const {
  double logDet = 0.0;

  if (arma::log_det_sympd(logDet, covMat)) {
    return (addSmallDiag ? std::max(logDet, lbDet) : logDet) * len;
  }

  if (addSmallDiag && epsilon > 0.0) {
    if (warnLevel()) {
      Rcpp::warning("`covMat` is singular! Consider increasing either `epsilon` or `minSize`");
      Rcpp::warning("Return the lower-bound `p*log(epsilon)*segLen`!");
    }
    return lbDet * len;
  }

  Rcpp::stop("`covMat` is singular! Consider using `addSmallDiag` option or increasing `minSize`!");
}

// ========================================================
//                     Cost_L1_cwMed
// ========================================================

Cost_L1_cwMed::Cost_L1_cwMed(const arma::mat& inputMat, bool warnOnce)
  : CostBase(warnOnce), X(inputMat) {
  nr = X.n_rows;
  nc = X.n_cols;
}

double Cost_L1_cwMed::eval(int start, int end) const {
  if (start >= end - 1) {
    return 0.0;
  }
  arma::mat segment = X.rows(start, end - 1);
  arma::rowvec med = arma::median(segment, 0);
  arma::mat residuals = arma::abs(segment.each_row() - med);
  return arma::accu(residuals);
}

Rcpp::List Cost_L1_cwMed::get_params(int start, int end) const {
  checkSegment(start, end);
  arma::rowvec med = arma::median(X.rows(start, end - 1), 0);
  return Rcpp::List::create(Rcpp::Named("median") = asNumeric(med));
}

// ========================================================
//                        Cost_L2
// ========================================================

Cost_L2::Cost_L2(const arma::mat& inputMat, bool warnOnce) : CostBase(warnOnce) {
  nr = inputMat.n_rows;
  nc = inputMat.n_cols;
  csX = getCumSumCpp(inputMat);
  csXsq = getCumSumCpp(arma::pow(inputMat, 2));
}

double Cost_L2::eval(int start, int end) const {
  if (start >= end - 1) {
    return 0.0;
  }
  int len = end - start;
  return arma::sum(csXsq.row(end) - csXsq.row(start)) -
    std::pow(arma::norm(csX.row(end) - csX.row(start), 2), 2) / len;
}

Rcpp::List Cost_L2::get_params(int start, int end) const {
  checkSegment(start, end);
  arma::rowvec mu = (csX.row(end) - csX.row(start)) / (end - start);
  return Rcpp::List::create(Rcpp::Named("mean") = asNumeric(mu));
}

// ========================================================
//                       Cost_SIGMA
// ========================================================

Cost_SIGMA::Cost_SIGMA(const arma::mat& inputMat, bool addSmallDiag, double epsilon, bool warnOnce)
  : CostBase(warnOnce), addSmallDiag_(addSmallDiag), epsilon_(epsilon) {
  nr = inputMat.n_rows;
  nc = inputMat.n_cols;
  lbDet = nc * std::log(epsilon);

  csX = getCumSumCpp(inputMat);
  csXXt.zeros(nc, nc, nr + 1);
  for (int i = 1; i <= nr; i++) {
    csXXt.slice(i) = csXXt.slice(i - 1) + inputMat.row(i - 1).t() * inputMat.row(i - 1);
  }
}

arma::mat Cost_SIGMA::segmentCov(int start, int end) const {
  int len = end - start;
  arma::mat covMat(nc, nc, arma::fill::zeros);
  if (len > 1) {
    arma::rowvec segSum = csX.row(end) - csX.row(start);
    arma::mat segOuter = csXXt.slice(end) - csXXt.slice(start);
    arma::vec mu = segSum.t() / len;
    covMat = (segOuter / len) - (mu * mu.t());
  }
  if (addSmallDiag_) {
    covMat.diag() += epsilon_;
  }
  return covMat;
}

double Cost_SIGMA::eval(int start, int end) const {
  return logDetCost(segmentCov(start, end), addSmallDiag_, epsilon_, lbDet, end - start);
}

Rcpp::List Cost_SIGMA::get_params(int start, int end) const {
  checkSegment(start, end);
  arma::rowvec mu = (csX.row(end) - csX.row(start)) / (end - start);
  return Rcpp::List::create(Rcpp::Named("mean") = asNumeric(mu),
                            Rcpp::Named("cov") = segmentCov(start, end));
}

// ========================================================
//                     RegressionCost
// ========================================================

RegressionCost::RegressionCost(bool warnOnce, const char* msgOnce, const char* msgEvery)
  : CostBase(warnOnce), msgOnce_(msgOnce), msgEvery_(msgEvery) {}

void RegressionCost::precompute(const arma::mat& Z, const arma::mat& Y, int offset) {
  csZtZ.zeros(J, J, nr + 1);
  csZtY.zeros(J, Y.n_cols, nr + 1);
  csYtY.zeros(Y.n_cols, Y.n_cols, nr + 1);

  for (arma::uword r = 0; r < Z.n_rows; ++r) {
    int i = offset + static_cast<int>(r);
    arma::rowvec zi = Z.row(r);
    arma::rowvec yi = Y.row(r);
    csZtZ.slice(i + 1) = csZtZ.slice(i) + zi.t() * zi;
    csZtY.slice(i + 1) = csZtY.slice(i) + zi.t() * yi;
    csYtY.slice(i + 1) = csYtY.slice(i) + yi.t() * yi;
  }
}

arma::mat RegressionCost::solveSegment(const arma::mat& ZtZ, const arma::mat& ZtY) const {
  arma::mat B;
  if (!arma::solve(B, ZtZ, ZtY, arma::solve_opts::no_approx + arma::solve_opts::likely_sympd)) {
    int level = warnLevel();
    if (level) {
      Rcpp::warning(level == 1 ? msgOnce_ : msgEvery_);
    }
    arma::solve(B, ZtZ, ZtY, arma::solve_opts::force_approx);
  }
  return B;
}

double RegressionCost::ssr(int start, int end) const {
  arma::mat ZtZ = csZtZ.slice(end) - csZtZ.slice(start);
  arma::mat ZtY = csZtY.slice(end) - csZtY.slice(start);
  arma::mat YtY = csYtY.slice(end) - csYtY.slice(start);
  arma::mat B = solveSegment(ZtZ, ZtY);
  return std::max(0.0, arma::trace(YtY) - arma::trace(B.t() * ZtY));
}

arma::mat RegressionCost::residualSSR(int start, int end, const arma::mat& B) const {
  arma::mat ZtY = csZtY.slice(end) - csZtY.slice(start);
  arma::mat YtY = csYtY.slice(end) - csYtY.slice(start);
  return YtY - B.t() * ZtY;
}

arma::mat RegressionCost::solveCoef(int start, int end, int nEff) const {
  arma::mat B(J, csZtY.n_cols);
  if (nEff < J) {
    B.fill(NA_REAL);
  } else {
    B = solveSegment(csZtZ.slice(end) - csZtZ.slice(start), csZtY.slice(end) - csZtY.slice(start));
  }
  return B;
}

Rcpp::List RegressionCost::coef(int start, int end, int nEff) const {
  return Rcpp::List::create(Rcpp::Named("coef") = solveCoef(start, end, nEff));
}

// ========================================================
//                     Cost_LinearL2
// ========================================================

Cost_LinearL2::Cost_LinearL2(const arma::mat& Y, const arma::mat& X, bool intercept_, bool warnOnce)
  : RegressionCost(warnOnce,
                   "System is singular. Switching to approximate solve.",
                   "Singular system encountered. Using force_approx."),
    intercept(intercept_) {

  if (Y.n_rows != X.n_rows) {
    Rcpp::stop("Number of observations in response and covariate matrices must match!");
  }

  nr = Y.n_rows;
  nc = Y.n_cols;
  J = X.n_cols + (intercept ? 1 : 0);

  if (nr < J) {
    Rcpp::stop("The full dataset contains not enough observations to fit a linear regression model!");
  }

  if (intercept) {
    precompute(arma::join_rows(arma::ones(nr), X), Y, 0);
  } else {
    precompute(X, Y, 0);
  }
}

double Cost_LinearL2::eval(int start, int end) const {
  if (start >= end - 1 || end - start < J) {
    return 0.0;
  }
  return ssr(start, end);
}

Rcpp::List Cost_LinearL2::get_params(int start, int end) const {
  checkSegment(start, end);
  return coef(start, end, end - start);
}

// ========================================================
//                        Cost_VAR
// ========================================================

Cost_VAR::Cost_VAR(const arma::mat& inputMat, int pVAR, bool warnOnce)
  : RegressionCost(warnOnce,
                   "Some systems seem singular! Switch to the approximate arma::solve()!",
                   "The system seems singular! Switch to the approximate arma::solve()!"),
    p(pVAR) {

  nr = inputMat.n_rows;
  nc = inputMat.n_cols;
  J = 1 + p * nc;

  if (p < 1) {
    Rcpp::stop("pVAR must be >= 1!");
  }
  if (nr - p < J) {
    Rcpp::stop("The full dataset contains not enough observations to fit VAR(p)!");
  }

  arma::mat Z(nr - p, J);
  Z.col(0).ones();
  for (int L = 0; L < p; L++) {
    Z.cols(1 + L * nc, (L + 1) * nc) = inputMat.rows(p - L - 1, nr - L - 2);
  }

  precompute(Z, inputMat.rows(p, nr - 1), p);
}

double Cost_VAR::eval(int start, int end) const {
  if (start > nr - p || end - start < J) {
    return 0.0;
  }
  return ssr(start, end);
}

Rcpp::List Cost_VAR::get_params(int start, int end) const {
  checkSegment(start, end);
  return coef(start, end, end - std::max(start, p));
}

// ========================================================
//                    Cost_LinearSIGMA
// ========================================================

Cost_LinearSIGMA::Cost_LinearSIGMA(const arma::mat& Y, const arma::mat& X, bool intercept_,
                                    bool addSmallDiag, double epsilon, bool warnOnce)
  : RegressionCost(warnOnce,
                   "System is singular. Switching to approximate solve.",
                   "Singular system encountered. Using force_approx."),
    intercept(intercept_), addSmallDiag_(addSmallDiag), epsilon_(epsilon) {

  if (Y.n_rows != X.n_rows) {
    Rcpp::stop("Number of observations in response and covariate matrices must match!");
  }

  nr = Y.n_rows;
  nc = Y.n_cols;
  J = X.n_cols + (intercept ? 1 : 0);
  lbDet = nc * std::log(epsilon);

  if (nr < J) {
    Rcpp::stop("The full dataset contains not enough observations to fit a linear regression model!");
  }

  if (intercept) {
    precompute(arma::join_rows(arma::ones(nr), X), Y, 0);
  } else {
    precompute(X, Y, 0);
  }
}

arma::mat Cost_LinearSIGMA::residualCov(int start, int end, const arma::mat& B) const {
  arma::mat covMat = residualSSR(start, end, B) / (end - start);
  if (addSmallDiag_) {
    covMat.diag() += epsilon_;
  }
  return covMat;
}

double Cost_LinearSIGMA::eval(int start, int end) const {
  if (start >= end - 1 || end - start < J) {
    return 0.0;
  }
  arma::mat B = solveCoef(start, end, end - start);
  return logDetCost(residualCov(start, end, B), addSmallDiag_, epsilon_, lbDet, end - start);
}

Rcpp::List Cost_LinearSIGMA::get_params(int start, int end) const {
  checkSegment(start, end);
  int len = end - start;
  arma::mat B = solveCoef(start, end, len);

  arma::mat covMat(nc, nc);
  if (len < J) {
    covMat.fill(NA_REAL);
  } else {
    covMat = residualCov(start, end, B);
  }

  return Rcpp::List::create(Rcpp::Named("coef") = B, Rcpp::Named("cov") = covMat);
}

// ========================================================
//                     Cost_LinearL1
// ========================================================

Cost_LinearL1::Cost_LinearL1(const arma::mat& Y_, const arma::mat& X, bool intercept,
                              double tol, int maxIter, bool warnOnce)
  : CostBase(warnOnce), tol_(tol), maxIter_(maxIter) {

  if (Y_.n_rows != X.n_rows) {
    Rcpp::stop("Number of observations in response and covariate matrices must match!");
  }

  nr = Y_.n_rows;
  nc = Y_.n_cols;
  J = X.n_cols + (intercept ? 1 : 0);

  if (nr < J) {
    Rcpp::stop("The full dataset contains not enough observations to fit a linear regression model!");
  }
  if (tol <= 0.0) {
    Rcpp::stop("`tol` must be a single positive value!");
  }
  if (maxIter < 1) {
    Rcpp::stop("`maxIter` must be at least 1!");
  }

  Y = Y_;
  Z = intercept ? arma::join_rows(arma::ones(nr), X) : X;
}

arma::vec Cost_LinearL1::fitColumnIRLS(int start, int end, int col, double* costOut) const {
  static const double kWeightFloor = 1e-6;

  arma::mat Zseg = Z.rows(start, end - 1);
  arma::vec yseg = Y.submat(start, col, end - 1, col);

  // Tracked across the initial solve and every IRLS iteration so a singular system is
  // reported once per call (like the other regression costs), not once per iteration --
  // up to maxIter warnings from a single eval() would be excessive, but staying fully
  // silent after the first iteration (unlike every other cost's singular-system handling)
  // would hide a real degradation even in "warn every call" mode.
  bool anySingular = false;

  arma::vec B;
  if (!arma::solve(B, Zseg.t() * Zseg, Zseg.t() * yseg,
                    arma::solve_opts::no_approx + arma::solve_opts::likely_sympd)) {
    anySingular = true;
    arma::solve(B, Zseg.t() * Zseg, Zseg.t() * yseg, arma::solve_opts::force_approx);
  }

  arma::vec resid = yseg - Zseg * B;
  double prevCost = arma::accu(arma::abs(resid));
  bool converged = false;

  for (int iter = 0; iter < maxIter_; ++iter) {
    arma::vec w = 1.0 / arma::clamp(arma::abs(resid), kWeightFloor, arma::datum::inf);
    arma::vec sqrtw = arma::sqrt(w);
    arma::mat Zw = arma::diagmat(sqrtw) * Zseg;
    arma::vec yw = sqrtw % yseg;

    arma::vec Bnew;
    if (!arma::solve(Bnew, Zw.t() * Zw, Zw.t() * yw,
                      arma::solve_opts::no_approx + arma::solve_opts::likely_sympd)) {
      anySingular = true;
      arma::solve(Bnew, Zw.t() * Zw, Zw.t() * yw, arma::solve_opts::force_approx);
    }

    B = Bnew;
    resid = yseg - Zseg * B;
    double cost = arma::accu(arma::abs(resid));

    if (std::abs(prevCost - cost) <= tol_ * (1.0 + prevCost)) {
      prevCost = cost;
      converged = true;
      break;
    }
    prevCost = cost;
  }

  if (anySingular) {
    int level = warnLevel();
    if (level) {
      Rcpp::warning(level == 1 ? "System is singular. Switching to approximate solve."
                                : "Singular system encountered. Using force_approx.");
    }
  }

  if (!converged) {
    int level = warnLevel();
    if (level) {
      Rcpp::warning(level == 1
        ? "IRLS did not converge within `maxIter`! Consider increasing `maxIter` or `tol`."
        : "IRLS did not converge within `maxIter` for this segment!");
    }
  }

  if (costOut) {
    *costOut = prevCost;
  }
  return B;
}

double Cost_LinearL1::eval(int start, int end) const {
  if (start >= end - 1 || end - start < J) {
    return 0.0;
  }
  double total = 0.0;
  for (int col = 0; col < nc; ++col) {
    double colCost = 0.0;
    fitColumnIRLS(start, end, col, &colCost);
    total += colCost;
  }
  return total;
}

Rcpp::List Cost_LinearL1::get_params(int start, int end) const {
  checkSegment(start, end);
  arma::mat B(J, nc);
  if (end - start < J) {
    B.fill(NA_REAL);
  } else {
    for (int col = 0; col < nc; ++col) {
      B.col(col) = fitColumnIRLS(start, end, col, nullptr);
    }
  }
  return Rcpp::List::create(Rcpp::Named("coef") = B);
}

// ========================================================
//                      Rcpp modules
// ========================================================

RCPP_EXPOSED_CLASS(Cost_L1_cwMed)
RCPP_EXPOSED_CLASS(Cost_L2)
RCPP_EXPOSED_CLASS(Cost_SIGMA)
RCPP_EXPOSED_CLASS(Cost_LinearL2)
RCPP_EXPOSED_CLASS(Cost_VAR)
RCPP_EXPOSED_CLASS(Cost_LinearSIGMA)
RCPP_EXPOSED_CLASS(Cost_LinearL1)

RCPP_MODULE(Cost_L1_cwMed_module) {
  Rcpp::class_<Cost_L1_cwMed>("Cost_L1_cwMed")
  .constructor<arma::mat, bool>()
  .method("eval", &Cost_L1_cwMed::eval, "Evaluate L1 cost on interval (start, end]")
  .method("get_params", &Cost_L1_cwMed::get_params, "Coordinate-wise median on (start, end]")
  .method("resetWarning", &resetWarningR<Cost_L1_cwMed>, "Set the status of warnOnce_");
}

RCPP_MODULE(Cost_L2_module) {
  Rcpp::class_<Cost_L2>("Cost_L2")
  .constructor<arma::mat, bool>()
  .method("eval", &Cost_L2::eval, "Evaluate L2 cost on interval (start, end]")
  .method("get_params", &Cost_L2::get_params, "Mean on (start, end]")
  .method("resetWarning", &resetWarningR<Cost_L2>, "Set the status of warnOnce_");
}

RCPP_MODULE(Cost_SIGMA_module) {
  Rcpp::class_<Cost_SIGMA>("Cost_SIGMA")
  .constructor<arma::mat, bool, double, bool>()
  .method("eval", &Cost_SIGMA::eval, "Evaluate SIGMA cost on interval (start, end]")
  .method("get_params", &Cost_SIGMA::get_params, "Mean and covariance on (start, end]")
  .method("resetWarning", &resetWarningR<Cost_SIGMA>, "Set the status of warnOnce_");
}

RCPP_MODULE(Cost_LinearL2_module) {
  Rcpp::class_<Cost_LinearL2>("Cost_LinearL2")
  .constructor<arma::mat, arma::mat, bool, bool>()
  .method("eval", &Cost_LinearL2::eval, "Evaluate linear regression cost on interval (start, end]")
  .method("get_params", &Cost_LinearL2::get_params, "Regression coefficients on (start, end]")
  .method("resetWarning", &resetWarningR<Cost_LinearL2>, "Set the status of warnOnce_");
}

RCPP_MODULE(Cost_VAR_module) {
  Rcpp::class_<Cost_VAR>("Cost_VAR")
  .constructor<arma::mat, int, bool>()
  .method("eval", &Cost_VAR::eval, "Evaluate VAR cost on interval (start, end]")
  .method("get_params", &Cost_VAR::get_params, "VAR coefficients on (start, end]")
  .method("resetWarning", &resetWarningR<Cost_VAR>, "Set the status of warnOnce_");
}

RCPP_MODULE(Cost_LinearSIGMA_module) {
  Rcpp::class_<Cost_LinearSIGMA>("Cost_LinearSIGMA")
  .constructor<arma::mat, arma::mat, bool, bool, double, bool>()
  .method("eval", &Cost_LinearSIGMA::eval, "Evaluate LinearSIGMA cost on interval (start, end]")
  .method("get_params", &Cost_LinearSIGMA::get_params, "Regression coefficients and residual covariance on (start, end]")
  .method("resetWarning", &resetWarningR<Cost_LinearSIGMA>, "Set the status of warnOnce_");
}

RCPP_MODULE(Cost_LinearL1_module) {
  Rcpp::class_<Cost_LinearL1>("Cost_LinearL1")
  .constructor<arma::mat, arma::mat, bool, double, int, bool>()
  .method("eval", &Cost_LinearL1::eval, "Evaluate LinearL1 (IRLS) cost on interval (start, end]")
  .method("get_params", &Cost_LinearL1::get_params, "Regression coefficients on (start, end]")
  .method("resetWarning", &resetWarningR<Cost_LinearL1>, "Set the status of warnOnce_");
}
