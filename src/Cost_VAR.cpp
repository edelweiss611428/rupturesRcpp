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

double Cost_VAR::effEvalCpp(int start, int end,  //Evaluate the cost of the segment (start,end] but counts from 1
                            bool addSmallDiag,      // unused
                            double epsilon) const {  // unused

  //p-step-behind
  // if(start < p){ //not used (1)
  //   start = p;
  // }

  //p-step-ahead
  if(start > nr-p-1){
    return 0;
  }

  end = end - 1;

  int len = end - start;
  if(len < J){ // (start, end] must have at least J points (for OLS)
    return 0; //Return 0
  }

  //p-step-ahead
  arma::mat Z_sub = Z_full.rows(start, end-p);
  arma::mat Y_sub = X.rows(start+p, end); //The longest valid section is y[p:(nr-1)]

  //p-step-behind
  // arma::mat Z_sub = Z_full.rows(start - p, end-p); //not used (1)
  // arma::mat Y_sub = X.rows(start, end); //not used (1)

  arma::mat regCoefs;

  try {
    // Standard solve()
    regCoefs = arma::solve(Z_sub, Y_sub);
  } catch (std::exception) {
    Rcpp::Rcout << "Standard solve failed! Return 0" << std::endl;
  }

  arma::mat errMat = Y_sub - Z_sub * regCoefs;

  return arma::accu(arma::square(errMat));

}


RCPP_EXPOSED_CLASS(Cost_VAR)
  RCPP_MODULE(Cost_VAR_module) {
    class_<Cost_VAR>("Cost_VAR")
      .constructor<arma::mat, int>()
      .method("effEvalCpp", &Cost_VAR::effEvalCpp,
    "Evaluate VAR cost on interval (start, end]")
    ;
  }

