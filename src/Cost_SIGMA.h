// Cost_SIGMA.h

#ifndef COST_SIGMA_H
#define COST_SIGMA_H

#include "CostBase.h"
#include <RcppArmadillo.h>
#include <limits>

using namespace Rcpp;

struct covariancePrecomputer {
  arma::mat S_k;                  // cumsum of data by row
  std::vector<arma::mat> Q_k;     // cumsum of outer products

  int nr;  // number of data points
  int nc;  // number of features

  covariancePrecomputer(const arma::mat& X);

  arma::mat covarianceComputer(int start, int end) const;
};

class Cost_SIGMA : public CostBase {
private:
  covariancePrecomputer preComp;

public:
  Cost_SIGMA(const arma::mat& inputMat);

  double effEvalCpp(int start, int end,
                    bool addSmallDiag = true,
                    double epsilon = 1e-6) const;
};

#endif // COST_SIGMA_H
