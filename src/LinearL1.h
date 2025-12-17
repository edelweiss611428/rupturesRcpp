#ifndef COST_LINEAR_L1_H
#define COST_LINEAR_L1_H

#include <RcppArmadillo.h>
#include "baseClass.h"

class Cost_LinearL1 : public CostBase {

private:
  arma::mat Yfull;
  arma::mat Zfull;

  int nr, pY, J;
  bool intercept;
  bool keepWarning;
  bool multitask;

  int max_iter;
  double tol, eps;

public:

  mutable bool warnOnce_;  // mutable to allow updates inside const method

  Cost_LinearL1(const arma::mat& Y,
                const arma::mat& X,
                bool intercept_ = true,
                bool multitask_ = false,
                bool warnOnce = true,
                int max_iter_ = 50,
                double tol_ = 1e-6,
                double eps_ = 1e-8);

  // Evaluate cost on interval (start, end]
  double eval(int start, int end) const override;

  // Reset warning behavior
  void resetWarning(bool reset) override;

  // Destructor (use default)
  virtual ~Cost_LinearL1() = default;

};

#endif
