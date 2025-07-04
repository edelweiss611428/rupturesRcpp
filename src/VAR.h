// VAR.h

#ifndef COST_VAR_H
#define COST_VAR_H

#include "baseClass.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

class Cost_VAR : public CostBase {

private:

  int p;
  int J;
  arma::mat X;
  arma::mat Z; //Stacked design matrix; each i-th row is the `p` previous obs of the (p+i)th obs
  std::vector<arma::mat> csZtZ;
  std::vector<arma::mat> csZtY;
  std::vector<arma::mat> csYtY;

public:

  Cost_VAR(const arma::mat& inputMat, const int& pVAR = 1);

  double eval(int start, int end) const;
};


#endif // COST_VAR_H
