// Cost.h

#ifndef COST_H
#define COST_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat getCumsumCpp(const arma::mat& X);

class Cost {
private:
  arma::mat X;
  arma::mat csX;    // cumsum(X)
  arma::mat csXsq;  // cumsum(X^2)

public:
  int nr;

  Cost(const arma::mat& inputMat);

  double effEvalCpp(int start, int end) const;
};


#endif // COST_H
