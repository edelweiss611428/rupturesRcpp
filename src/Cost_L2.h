// Cost_L2.h

#ifndef COST_L2_H
#define COST_L2_H

#include "CostBase.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat getCumsumCpp(const arma::mat& X);

class Cost_L2 : public CostBase {

private:
  arma::mat csX;    // cumsum(X)
  arma::mat csXsq;  // cumsum(X^2)

public:

  Cost_L2(const arma::mat& inputMat);

  double effEvalCpp(int start, int end,
                    bool addSmallDiag = true,
                    double epsilon = 1e-6) const;
};


#endif // COST_L2_H
