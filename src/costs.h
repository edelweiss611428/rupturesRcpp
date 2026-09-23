#ifndef RUPTURES_COSTS_H
#define RUPTURES_COSTS_H

#include <RcppArmadillo.h>

// All costs act on the segment (start, end], i.e. rows start, ..., end - 1 (0-indexed).

// ========================================================
//                        CostBase
// ========================================================

class CostBase {
public:
  int nr = 0;
  int nc = 0;

  explicit CostBase(bool warnOnce = true) { resetWarning(warnOnce); }
  virtual ~CostBase() = default;

  virtual double eval(int start, int end) const = 0;
  virtual Rcpp::List get_params(int start, int end) const = 0;

  // true: warn once (inside fit/predict); false: warn on every call
  void resetWarning(bool reset) {
    warnOnce_ = reset;
    keepWarning = !reset;
  }

  void checkSegment(int start, int end) const;

protected:
  mutable bool warnOnce_;
  bool keepWarning;

  // 0 = silent, 1 = first warning since resetWarning(true), 2 = warn on every call
  int warnLevel() const;

  // len * log|covMat|, falling back to len * lbDet (with a warning) when covMat is singular
  // and addSmallDiag/epsilon allow it; stops otherwise. `covMat` must already include epsilon
  // regularization if addSmallDiag is true. Shared by Cost_SIGMA and Cost_LinearSIGMA.
  double logDetCost(const arma::mat& covMat, bool addSmallDiag, double epsilon,
                     double lbDet, int len) const;
};

// ========================================================
//                     Cost_L1_cwMed
// ========================================================

class Cost_L1_cwMed : public CostBase {
public:
  Cost_L1_cwMed(const arma::mat& inputMat, bool warnOnce = true);
  double eval(int start, int end) const override;
  Rcpp::List get_params(int start, int end) const override;  // median

private:
  arma::mat X;
};

// ========================================================
//                        Cost_L2
// ========================================================

class Cost_L2 : public CostBase {
public:
  Cost_L2(const arma::mat& inputMat, bool warnOnce = true);
  double eval(int start, int end) const override;
  Rcpp::List get_params(int start, int end) const override;  // mean

private:
  arma::mat csX;
  arma::mat csXsq;
};

// ========================================================
//                       Cost_SIGMA
// ========================================================

class Cost_SIGMA : public CostBase {
public:
  Cost_SIGMA(const arma::mat& inputMat, bool addSmallDiag = true,
             double epsilon = 1e-6, bool warnOnce = true);
  double eval(int start, int end) const override;
  Rcpp::List get_params(int start, int end) const override;  // mean, cov

private:
  arma::mat csX;     // cumsum of rows
  arma::cube csXXt;  // cumsum of outer products
  bool addSmallDiag_;
  double epsilon_;
  double lbDet;      // p * log(epsilon)

  arma::mat segmentCov(int start, int end) const;  // MLE covariance (+ epsilon * I if addSmallDiag)
};

// ========================================================
//      RegressionCost: shared base for LinearL2 and VAR
// ========================================================

// SSR of the multivariate least-squares fit Y ~ Z on a segment, computed from
// cumulative sums of Z'Z, Z'Y and Y'Y.
class RegressionCost : public CostBase {
protected:
  int J = 0;  // number of regression coefficients per response
  arma::cube csZtZ;
  arma::cube csZtY;
  arma::cube csYtY;

  RegressionCost(bool warnOnce, const char* msgOnce, const char* msgEvery);

  // Z.row(r) and Y.row(r) belong to observation offset + r
  void precompute(const arma::mat& Z, const arma::mat& Y, int offset);

  double ssr(int start, int end) const;
  arma::mat solveCoef(int start, int end, int nEff) const;  // B_hat (J x pY), NA-filled if nEff < J
  Rcpp::List coef(int start, int end, int nEff) const;  // {coef: solveCoef(...)}
  arma::mat residualSSR(int start, int end, const arma::mat& B) const;  // R'R (q x q) given a precomputed B_hat

private:
  const char* msgOnce_;
  const char* msgEvery_;

  arma::mat solveSegment(const arma::mat& ZtZ, const arma::mat& ZtY) const;
};

class Cost_LinearL2 : public RegressionCost {
public:
  Cost_LinearL2(const arma::mat& Y, const arma::mat& X,
                bool intercept = true, bool warnOnce = true);
  double eval(int start, int end) const override;
  Rcpp::List get_params(int start, int end) const override;  // coef: J x pY, intercept first

private:
  bool intercept;
};

class Cost_VAR : public RegressionCost {
public:
  Cost_VAR(const arma::mat& inputMat, int pVAR, bool warnOnce = true);
  double eval(int start, int end) const override;
  Rcpp::List get_params(int start, int end) const override;  // coef: (1 + p*nc) x nc, intercept then lags 1..p

private:
  int p;
};

// ========================================================
//                    Cost_LinearSIGMA
// ========================================================

// Piecewise linear regression with a segment-varying noise covariance: like Cost_LinearL2,
// but the cost is (b-a)*log|Sigma_hat| on the OLS residual covariance (like Cost_SIGMA),
// instead of the residual sum of squares.
class Cost_LinearSIGMA : public RegressionCost {
public:
  Cost_LinearSIGMA(const arma::mat& Y, const arma::mat& X, bool intercept = true,
                    bool addSmallDiag = true, double epsilon = 1e-6, bool warnOnce = true);
  double eval(int start, int end) const override;
  Rcpp::List get_params(int start, int end) const override;  // coef, cov

private:
  bool intercept;
  bool addSmallDiag_;
  double epsilon_;
  double lbDet;  // nc * log(epsilon)

  // MLE residual covariance (+ epsilon * I if addSmallDiag) given a precomputed B_hat
  arma::mat residualCov(int start, int end, const arma::mat& B) const;
};

// ========================================================
//                     Cost_LinearL1
// ========================================================

// Piecewise linear regression under L1 (least absolute deviations) loss, fit per response
// column via Iteratively Reweighted Least Squares (IRLS): at each iteration, observations are
// weighted by 1/max(|residual|, delta) and a weighted least-squares problem is solved, which
// converges to the L1 minimiser. Unlike the other regression costs, this has no O(1)-per-segment
// closed form derivable from cumulative sums (the IRLS weights are segment- and fit-specific) --
// each eval()/get_params() call re-fits IRLS on the segment's raw rows, costing
// O(len * J^2 * iterations), not O(1).
class Cost_LinearL1 : public CostBase {
public:
  Cost_LinearL1(const arma::mat& Y, const arma::mat& X, bool intercept = true,
                double tol = 1e-6, int maxIter = 50, bool warnOnce = true);
  double eval(int start, int end) const override;
  Rcpp::List get_params(int start, int end) const override;  // coef: J x pY, intercept first

private:
  arma::mat Z;  // design matrix (intercept column prepended if requested)
  arma::mat Y;  // response
  int J;
  double tol_;
  int maxIter_;

  // IRLS fit of one response column on rows [start, end); returns the coefficient vector and,
  // if costOut is non-null, the resulting sum of absolute residuals.
  arma::vec fitColumnIRLS(int start, int end, int col, double* costOut) const;
};

// ========================================================
//                      Cost_RFunc
// ========================================================

// User-defined cost, backed by an R closure supplied at runtime. `costFun` is called as
// costFun(segment, start, end): `segment` is the raw rows start, ..., end - 1 of `inputMat` (as an
// R matrix), and `start`/`end` are the same 0-indexed, half-open (start, end] bounds passed to
// eval() -- this lets a closure align `segment` against any externally-captured, position-indexed
// data (e.g. `externalSeries[(start+1):end]`) without the package needing to know that data exists.
// costFun must return a single numeric value; it is called as-is, with no special-casing of trivial
// (length <= 1) segments. If `paramFun` is supplied, get_params() calls it the same way and wraps
// its return value under "params"; otherwise get_params() returns an empty list. Because every call
// to eval()/get_params() crosses back into R, this cost is substantially slower per call than the
// built-in, closed-form costs -- expected, since the cost logic itself lives in R.
class Cost_RFunc : public CostBase {
public:
  Cost_RFunc(const arma::mat& inputMat, Rcpp::Function costFun,
             Rcpp::Nullable<Rcpp::Function> paramFun = R_NilValue, bool warnOnce = true);
  double eval(int start, int end) const override;
  Rcpp::List get_params(int start, int end) const override;

private:
  arma::mat X;
  mutable Rcpp::Function costFun_;
  Rcpp::Nullable<Rcpp::Function> paramFun_;
};

#endif // RUPTURES_COSTS_H
