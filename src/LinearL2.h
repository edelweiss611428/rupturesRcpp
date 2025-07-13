#ifndef COST_LINEARL2_H
#define COST_LINEARL2_H

#include <RcppArmadillo.h>
#include "baseClass.h"

// ========================================================
//                  Cost_LinearL2 class
// ========================================================

class Cost_LinearL2 : public CostBase {

private:
  int pY;         // number of response variables
  int pX;         // number of covariates (excluding intercept)
  int J;          // number of regression parameters
  bool intercept; // whether intercept is included
  bool keepWarning;

  arma::cube csXtX; // cumulative sum of XtX
  arma::cube csXtY; // cumulative sum of XtY
  arma::cube csYtY; // cumulative sum of YtY

public:
  mutable bool warnOnce_;  // mutable to allow updates inside const method

  // Constructor
  Cost_LinearL2(const arma::mat& Y,
                const arma::mat& X,
                bool intercept = true,
                bool warnOnce = true);

  // Evaluate cost on interval (start, end]
  double eval(int start, int end) const override;

  // Reset warning behavior
  void resetWarning(bool reset) override;

  // Destructor (use default)
  virtual ~Cost_LinearL2() = default;
};

#endif // COST_LINEARL2_H
