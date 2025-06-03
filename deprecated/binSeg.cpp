#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include "../src/Cost.h"

inline List miniOptCpp(const Cost& Xnew, const int& start, const int& end) {

  int len = end - start;

  if(len == 1 or len == 0){
    return List::create(
      Named("valid") = false,
      Named("lErr") = std::numeric_limits<double>::infinity(),
      Named("rErr") = std::numeric_limits<double>::infinity()
    );
  } else if(len == 2){
    return List::create(
      Named("valid") = true,
      Named("err") = 0,
      Named("lErr") = 0,
      Named("rErr") = 0,
      Named("start") = start,
      Named("cp") = start+1,
      Named("end") = end
    );
  }

  double minErr = std::numeric_limits<double>::infinity();
  int cp;
  int tempCp;
  double err;
  double lErr;
  double rErr;
  double minlErr;
  double minrErr;

  for(int i = 0; i < (len-1); i++){
    tempCp = start+i+1;
    lErr = Xnew.effEvalCpp(start,tempCp);
    rErr = Xnew.effEvalCpp(tempCp,end);
    err = lErr + rErr;

    if(err < minErr){
      minErr = err;
      minlErr = lErr;
      minrErr = rErr;
      cp = tempCp;
    }
  }

  return List::create(
    Named("valid") = true,
    Named("err") = minErr,
    Named("lErr") = minlErr,
    Named("rErr") = minrErr,
    Named("start") = start,
    Named("cp") = cp,
    Named("end") = end
  );

}

// Slow version: considers all k partitions in each iteration,
// without checking whether a partition has already been split before.

// // [[Rcpp::export]]
List binSegCpp_v1(const arma::mat& tsMat, const int& maxNRegimes) {

  Cost Xnew(tsMat);
  int nr = Xnew.nr;
  List cpd0 = miniOptCpp(Xnew, 0, nr);
  if(maxNRegimes > nr){
    stop("The maximum number of regimes must be less than or equal to the number of observations.");;
  }

  if(nr == 1){
    stop("There is no changepoint as there is only one observation!");
  }

  if (maxNRegimes == 2){
    return cpd0;
  }


  IntegerVector changePoints(maxNRegimes-1);
  IntegerMatrix regimes(2, maxNRegimes);
  IntegerVector savedCps(maxNRegimes);

  NumericVector currentErrs(maxNRegimes);
  NumericVector tempErrs(maxNRegimes);
  NumericVector leftErrs(maxNRegimes);
  NumericVector rightErrs(maxNRegimes);
  NumericVector gains(maxNRegimes);


  changePoints[0] = Rcpp::as<int>(cpd0["cp"]);

  currentErrs[0] = Rcpp::as<double>(cpd0["lErr"]);
  currentErrs[1] = Rcpp::as<double>(cpd0["rErr"]);

  //1d indexing for matrices
  regimes[0] = 0;
  regimes[1] = Rcpp::as<int>(cpd0["cp"]);
  regimes[2] = Rcpp::as<int>(cpd0["cp"]);
  regimes[3] = nr;

  int nRegimes = 2;

  while(nRegimes < maxNRegimes){

    int bestIdx;

    double maxGain = -std::numeric_limits<double>::infinity();

    for(int i = 0; i < nRegimes; i++){

      List cpdi = miniOptCpp(Xnew, regimes[2*i], regimes[2*i+1]);


      if(not cpdi["valid"]){
        gains[i] = -std::numeric_limits<double>::infinity();
        continue;
      }

      savedCps[i] = Rcpp::as<int>(cpdi["cp"]);
      leftErrs[i] = Rcpp::as<double>(cpdi["lErr"]);
      rightErrs[i] = Rcpp::as<double>(cpdi["rErr"]);
      gains[i] = currentErrs[i] - Rcpp::as<double>(cpdi["err"]);

      if(gains[i] > maxGain){
        maxGain = gains[i];
        bestIdx = i;
      }
    }

    changePoints[nRegimes-1] = savedCps[bestIdx];
    IntegerVector aIdx = IntegerVector::create(bestIdx, nRegimes);

    currentErrs[aIdx] = NumericVector::create(leftErrs[bestIdx], rightErrs[bestIdx]);
    int tEnd = regimes[2*bestIdx+1];
    regimes[2*bestIdx+1] = savedCps[bestIdx];
    regimes[2*nRegimes] = savedCps[bestIdx];
    regimes[2*nRegimes+1] = tEnd;
    nRegimes++;

  }

  return List::create(
    Named("regimes") = regimes,
    Named("changePoints") = changePoints
  );
}


//Faster version: considers all k partitions in each iteration but
// also checks whether a partition has already been split before.
// // [[Rcpp::export]]
List binSegCpp_v2(const arma::mat& tsMat, const int& maxNRegimes) {

  Cost Xnew(tsMat);
  int nr = Xnew.nr;
  List cpd0 = miniOptCpp(Xnew, 0, nr);
  if(maxNRegimes > nr){
    stop("The maximum number of regimes must be less than or equal to the number of observations.");;
  }

  if(nr == 1){
    stop("There is no changepoint as there is only one observation!");
  }

  if (maxNRegimes == 2){
    return cpd0;
  }


  IntegerVector changePoints(maxNRegimes-1);
  IntegerMatrix regimes(2, maxNRegimes);
  IntegerVector savedCps(maxNRegimes);

  LogicalVector visited(maxNRegimes, false);
  NumericVector cost(maxNRegimes);

  cost[0] = Xnew.effEvalCpp(0,nr);
  cost[1] = Rcpp::as<double>(cpd0["err"]);

  NumericVector currentErrs(maxNRegimes);
  NumericVector tempErrs(maxNRegimes);
  NumericVector leftErrs(maxNRegimes);
  NumericVector rightErrs(maxNRegimes);
  NumericVector gains(maxNRegimes);

  changePoints[0] = Rcpp::as<int>(cpd0["cp"]);

  currentErrs[0] = Rcpp::as<double>(cpd0["lErr"]);
  currentErrs[1] = Rcpp::as<double>(cpd0["rErr"]);

  //1d indexing for matrices
  regimes[0] = 0;
  regimes[1] = Rcpp::as<int>(cpd0["cp"]);
  regimes[2] = Rcpp::as<int>(cpd0["cp"]);
  regimes[3] = nr;

  int nRegimes = 2;

  while(nRegimes < maxNRegimes){

    int bestIdx;

    double maxGain = -std::numeric_limits<double>::infinity();

    for(int i = 0; i < nRegimes; i++){


      if(not visited[i]){

        visited[i] = true;
        List cpdi = miniOptCpp(Xnew, regimes[2*i], regimes[2*i+1]);


        if(not cpdi["valid"]){
          gains[i] = -std::numeric_limits<double>::infinity();
          continue;
        }

        savedCps[i] = Rcpp::as<int>(cpdi["cp"]);
        leftErrs[i] = Rcpp::as<double>(cpdi["lErr"]);
        rightErrs[i] = Rcpp::as<double>(cpdi["rErr"]);
        gains[i] = currentErrs[i] - Rcpp::as<double>(cpdi["err"]);
      }

      if(gains[i] > maxGain){
        maxGain = gains[i];
        bestIdx = i;
      }
    }

    changePoints[nRegimes-1] = savedCps[bestIdx];
    IntegerVector aIdx = IntegerVector::create(bestIdx, nRegimes);
    visited[aIdx] = false;
    currentErrs[aIdx] = NumericVector::create(leftErrs[bestIdx], rightErrs[bestIdx]);
    int tEnd = regimes[2*bestIdx+1];
    regimes[2*bestIdx+1] = savedCps[bestIdx];
    regimes[2*nRegimes] = savedCps[bestIdx];
    regimes[2*nRegimes+1] = tEnd;

    cost[nRegimes] = cost[nRegimes-1] - maxGain;
    nRegimes++;

  }

  return List::create(
    Named("regimes") = regimes,
    Named("changePoints") = changePoints,
    Named("cost") = cost
  );


}


