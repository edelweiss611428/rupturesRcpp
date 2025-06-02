#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include "Cost.h"

// New structure for a segment
struct Segment {
  int start;
  int end;
  bool valid;
  int cp;
  double gain;
  double lErr;
  double rErr;
  double err;

  // Max-heap: higher gain has higher priority
  bool operator<(const Segment& other) const {
    return gain < other.gain;
  }
};


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

// [[Rcpp::export]]
List slowBinSegCpp(const arma::mat& tsMat, const int& maxNRegimes) {

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



// [[Rcpp::export]]
List fastBinSegCpp(const arma::mat& tsMat, const int& maxNRegimes) {

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



inline Segment miniOptHeapCpp(const Cost& Xnew, const int& start, const int& end,
                              double totalErr = -1) {

  int len = end - start;

  if(totalErr < 0){
    totalErr = Xnew.effEvalCpp(start, end);
  }

  if(len == 1 or len == 0){

    return Segment{start, end, true, start,
                   - std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};
  } else if(len == 2){

    return Segment{start, end, true, start+1,
                   totalErr,
                   0,
                   0,
                   0};
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

  return Segment{start, end, true, cp,
                 totalErr - minErr,
                 minlErr,
                 minrErr,
                 minErr};
}

// [[Rcpp::export]]
List fasterBinSegCpp(const arma::mat& tsMat, const int& maxNRegimes) {

  Cost Xnew(tsMat);
  int nr = Xnew.nr;

  NumericVector cost(maxNRegimes);
  Segment seg0 = miniOptHeapCpp(Xnew, 0, nr, Xnew.effEvalCpp(0,nr));

  IntegerVector changePoints(maxNRegimes-1);

  std::priority_queue<Segment> heap;
  heap.push(Segment{0, nr, true, seg0.cp, seg0.gain, seg0.lErr, seg0.rErr, seg0.err});

  int nRegimes = 2;
  int idx = 0;

  do {

    Segment bestSeg = heap.top();
    heap.pop();

    changePoints[idx] = bestSeg.cp;
    idx++;
    Segment leftSeg = miniOptHeapCpp(Xnew, bestSeg.start, bestSeg.cp, bestSeg.lErr);
    Segment rightSeg = miniOptHeapCpp(Xnew, bestSeg.cp, bestSeg.end, bestSeg.rErr);

    heap.push(leftSeg);
    heap.push(rightSeg);

    nRegimes++;

  } while(nRegimes < maxNRegimes);

  Segment bestSeg = heap.top();
  changePoints[idx] = bestSeg.cp;


  return List::create(
    Named("cp") = changePoints
  );


}
