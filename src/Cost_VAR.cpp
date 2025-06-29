#include "Cost_VAR.h"

// Cost_VAR constructor
Cost_VAR::Cost_VAR(const arma::mat& tsMat, const int& pVAR)  : p(pVAR), X(tsMat) {
  nr = X.n_rows;
  nc = X.n_cols;
  J = 1 + p * nc;

  if (nr - p < J){
    stop("Not enough observations to fit VAR(p)");
  }

  Z_full.set_size(nr-p, J);
  Z_full.col(0).ones();

  for (int L = 0; L < p; L++) {
    Z_full.cols(1 + L * nc, (L + 1) * nc) = X.rows(p - L - 1, nr - L - 2);
  }


};

double Cost_VAR::effEvalCpp(int start, int end,
                            bool addSmallDiag,      // ignored
                            double epsilon) const {  // ignored


  if(start < p){ //to account for the fact that the longest valid section is y[p:(nr-1)] - here idx from 0
    start = p;
  }

  end = end - 1;

  int len = end - start;
  if(len < J){
    return 0; //in case numerical failure will occur, return 0
  }


  arma::mat Z_sub = Z_full.rows(start - p, end-p);
  arma::mat Y_sub = X.rows(start, end);
  arma::mat regCoefs;
  try {
    regCoefs = arma::solve(Z_sub, Y_sub);
  } catch (...) {
    return 0;  //in case failure will occur, return 0
  }
  arma::mat errMat = Y_sub - Z_sub * regCoefs;

  return arma::accu(arma::square(errMat));

}


RCPP_MODULE(CostVAR_module) {
  class_<Cost_VAR>("CostVAR")

  .constructor<arma::mat, int>()

  .method("effEvalCpp", &Cost_VAR::effEvalCpp,
  "Evaluate VAR cost on interval (start, end]");
}
