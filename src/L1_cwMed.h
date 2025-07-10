#ifndef COST_L1_CW_MEDIAN_H
#define COST_L1_CW_MEDIAN_H

#include "baseClass.h"
#include <RcppArmadillo.h>

class Cost_L1_cwMed : public CostBase {
private:
  arma::mat X;
  bool warnOnce_;
  bool keepWarning;

public:

  Cost_L1_cwMed(const arma::mat& inputMat, bool warnOnce = true);
  double eval(int start, int end) const override;
  void resetWarning(bool reset) override;

};

#endif
