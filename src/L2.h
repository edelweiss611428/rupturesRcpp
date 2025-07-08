// L2.h

#ifndef COST_L2_H
#define COST_L2_H

#include "baseClass.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

// ========================================================
//                   Utility: getCumSumCpp
// ========================================================

arma::mat getCumSumCpp(const arma::mat& X);

// ========================================================
//                       Cost_L2 class
// ========================================================

class Cost_L2 : public CostBase {

private:
  arma::mat csX;    // cumsum(X)
  arma::mat csXsq;  // cumsum(X^2)

public:

  mutable bool warnOnce_;
  bool keepWarning;

  Cost_L2(const arma::mat& inputMat, bool warnOnce = true);
  double eval(int start, int end) const override;
  void resetWarning(bool reset) override;

};


#endif // COST_L2_H
