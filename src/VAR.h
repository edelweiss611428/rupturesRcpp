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
  arma::mat Z_full; //Stacked design matrix; each i-th row is `p` previous observations of the (p+i)th observation

public:

  Cost_VAR(const arma::mat& inputMat, const int& pVAR = 1);

  double eval(int start, int end) const;
};


#endif // COST_VAR_H
