#include "VAR.h"


// [[Rcpp::depends(RcppArmadillo)]]

Cost_VAR::Cost_VAR(const arma::mat& inputMat, const int& pVAR,
                   const bool&warnOnce){

  nr = inputMat.n_rows;
  nc = inputMat.n_cols;
  warnOnce_ = warnOnce;
  bool keepWarning = not warnOnce;
  (void)keepWarning;

  p = pVAR;
  J = 1 + p * nc;

  if(p < 1){
    stop("pVAR must be >= 1!");
  }

  if (nr - p < J){
    stop("Not enough observations to fit VAR(p)");
  }

  Z.set_size(nr-p, J);
  Z.col(0).ones();

  //Initialise precomputed matrices

  csZtZ.resize(nr+1, arma::mat(J,J,arma::fill::zeros));
  csZtY.resize(nr+1, arma::mat(J,nc,arma::fill::zeros));
  csYtY.resize(nr+1, arma::mat(nc,nc,arma::fill::zeros));

  for (int L = 0; L < p; L++) {
    Z.cols(1 + L * nc, (L + 1) * nc) = inputMat.rows(p - L - 1, nr - L - 2);
  }


  //Precomputation
  arma::rowvec Zi;
  arma::rowvec Yi;

  for (int i = p; i < nr; i++) {
    Zi = Z.row(i-p);
    Yi = inputMat.row(i);
    csZtZ[i+1] = csZtZ[i] + Zi.t() * Zi;
    csZtY[i+1] = csZtY[i] + Zi.t() * Yi;
    csYtY[i+1] = csYtY[i] + Yi.t() * Yi;
  }

}

double Cost_VAR::eval(int start, int end) {


  //p-step-ahead
  if(start > nr-p){
    return 0.0;
  }

  int len = end - start;

  if(len < J){ // (start, end] must have at least J points
    return 0.0; //Return 0
  }

  //p-step-ahead
  arma::mat G = csZtZ[end] - csZtZ[start];
  arma::mat H = csZtY[end] - csZtY[start];
  arma::mat Syy = csYtY[end] - csYtY[start];

  // Solve G * B = H
  arma::mat B;

  bool success = arma::solve(B, G, H);
  if (!success) {

    if(warnOnce_){
      warning("Some systems seem singular! Consider increasing either `epsilon` or `minSize`!");
      warnOnce_ = false;
    }

    if(keepWarning){
      warning("Some systems seem singular! Consider increasing either `epsilon` or `minSize`!");
    }

    return 0.0;

  }

  double trace_Syy = arma::trace(Syy);
  double trace_BtH = arma::trace(B.t() * H);

  return trace_Syy - trace_BtH;

}

RCPP_EXPOSED_CLASS(Cost_VAR)
  RCPP_MODULE(Cost_VAR_module) {
    Rcpp::class_<Cost_VAR>("Cost_VAR")
    .constructor<arma::mat, int, bool>()
    .method("eval", &Cost_VAR::eval,
    "Evaluate VAR cost on interval (start, end]")
    ;
  }
