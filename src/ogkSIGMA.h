#ifndef COST_OGKSIGMA_H
#define COST_OGKSIGMA_H

#include "baseClass.h"
#include <RcppArmadillo.h>

// ========================================================
//                   Cost_ogkSIGMA class
// ========================================================

class Cost_ogkSIGMA : public CostBase {
private:
  arma::mat X;             // Input data matrix
  bool warnOnce_;          // Single warning mode
  bool keepWarning;        // Control repeated warnings

public:
  // Constructor
  Cost_ogkSIGMA(const arma::mat& inputMat, bool warnOnce = true);

  // Cost evaluation between [start, end)
  double eval(int start, int end) const override;

  // Reset warning flag
  void resetWarning(bool reset) override;
};

#endif  // COST_OGKSIGMA_H
