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
  Rcpp::List coef(int start, int end, int nEff) const;  // NA if nEff < J
  arma::mat residualSSR(int start, int end) const;  // full R'R (q x q), R = Y - Z*coef

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

  arma::mat residualCov(int start, int end) const;  // MLE residual covariance (+ epsilon * I if addSmallDiag)
};

#endif // RUPTURES_COSTS_H
