// Cost_VAR.h

#ifndef Cost_VAR_H
#define Cost_VAR_H

#include "CostBase.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

class Cost_VAR : public CostBase {

private:
  int p;
  int J;
  arma::mat X;
  arma::mat Z_full; //stacked design matrix

public:

  Cost_VAR(const arma::mat& inputMat, const int& pVar);

  double effEvalCpp(int start, int end,
                    bool addSmallDiag = true, //unused
                    double epsilon = 1e-6) const;
};


#endif // Cost_VAR_H
