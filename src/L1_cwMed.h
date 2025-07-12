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

  //Constructor
  Cost_L1_cwMed(const arma::mat& inputMat, bool warnOnce = true);

  // Evaluate cost on interval (start, end]
  double eval(int start, int end) const override;

  // Reset warning behavior
  void resetWarning(bool reset) override;

  // Destructor (use default)
  virtual ~Cost_L1_cwMed() = default;

};

#endif
