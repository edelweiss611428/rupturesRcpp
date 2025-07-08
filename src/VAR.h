// VAR.h

#ifndef COST_VAR_H
#define COST_VAR_H

#include "baseClass.h"
#include <RcppArmadillo.h>
using namespace Rcpp;


// ========================================================
//                      Cost_VAR class
// ========================================================

class Cost_VAR : public CostBase {

private:

  int p;
  int J;
  arma::mat Z; //Stacked design matrix; each i-th row is the `p` previous obs of the (p+i)th obs
  std::vector<arma::mat> csZtZ;
  std::vector<arma::mat> csZtY;
  std::vector<arma::mat> csYtY;

public:

  mutable bool warnOnce_;
  bool keepWarning;

  Cost_VAR(const arma::mat& inputMat, int pVAR = 1, bool warnOnce = true);

  double eval(int start, int end) const override;
  void resetWarning(bool reset) override;
};


#endif // COST_VAR_H
