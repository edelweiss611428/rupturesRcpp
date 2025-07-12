// SIGMA.h

#ifndef COST_SIGMA_H
#define COST_SIGMA_H

#include "baseClass.h"
#include <RcppArmadillo.h>
using namespace Rcpp;


// ========================================================
//           Utility: covariancePrecomputer class
// ========================================================

struct covariancePrecomputer {
  arma::mat S_k;                  // cumsum of data by row
  std::vector<arma::mat> Q_k;     // cumsum of outer products

  int nr;  // Number of observations
  int nc;  // Number of features

  covariancePrecomputer(const arma::mat& X);

  arma::mat covarianceComputer(int start, int end) const;
};

// ========================================================
//                    Cost_SIGMA class
// ========================================================

class Cost_SIGMA : public CostBase {

private:
  covariancePrecomputer preComp;
  bool addSmallDiag_;
  double epsilon_;


public:

  mutable bool warnOnce_;
  bool keepWarning;
  double lbDet; //lower bound for the determinant


  // Constructor
  Cost_SIGMA(const arma::mat& inputMat,
             bool addSmallDiag = true, double epsilon = 1e-6, bool warnOnce = true);

  // Evaluate cost on interval (start, end]
  double eval(int start, int end) const override;

  // Reset warning behavior
  void resetWarning(bool reset) override;

  // Destructor (use default)
  virtual ~Cost_SIGMA() = default;
};


#endif // COST_SIGMA_H
