#include "L1_cwMed.h"

using namespace Rcpp;

// ========================================================
//       Cost_L1 Class based on Coordinate-wise Median
// ========================================================

Cost_L1_cwMed::Cost_L1_cwMed(const arma::mat& inputMat, bool warnOnce) {
  X = inputMat;
  nr = X.n_rows;
  warnOnce_ = warnOnce;
  keepWarning = !warnOnce;
}

double Cost_L1_cwMed::eval(int start, int end) const {

  if (start >= end - 1) {
    return 0.0;
  }

  arma::mat segment = X.rows(start, end - 1);
  arma::rowvec med = arma::median(segment, 0);  // coordinate-wise median
  arma::mat residuals = arma::abs(segment.each_row() - med);
  return arma::accu(residuals);
}

void Cost_L1_cwMed::resetWarning(bool reset) {
  warnOnce_ = reset;
  keepWarning = !reset;
}



RCPP_EXPOSED_CLASS(Cost_L1_cwMed)
  RCPP_MODULE(Cost_L1_cwMed_module) {
    Rcpp::class_<Cost_L1_cwMed>("Cost_L1_cwMed")
    .constructor<arma::mat, bool>()
    .method("eval", &Cost_L1_cwMed::eval,
    "Evaluate L2 cost on interval (start, end]")
    .method("resetWarning", &Cost_L1_cwMed::resetWarning,
    "Set the status of warnOnce_");
  }


