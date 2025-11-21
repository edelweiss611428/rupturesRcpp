#include <RcppArmadillo.h>
#include <queue>
#include <limits>
#include "VAR.h"
#include "L2.h"
#include "SIGMA.h"
#include "L1_cwMed.h"
#include "LinearL2.h"
#include "baseClass.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// ========================================================
//                  Utility: binSegPredCpp
// ========================================================

IntegerVector binSegPredCpp(const IntegerVector& bkps,
                            const NumericVector& cost,
                            double penalty = 0){

  // # nocov start
  if(penalty <0){
    stop("Penalty should be non-negative!"); //tested in R
  }
  // # nocov end

  NumericVector penCost = clone(cost);

  for(int i = 0; i < cost.size(); i++){
    penCost[i] = cost[i] + penalty*i;
  }

  int minIdx = which_min(penCost);

  if(minIdx == 0){
    return IntegerVector();
  }

  return bkps[Range(0, minIdx-1)];

}


// ========================================================
//                  Utility: Segment class
// ========================================================

struct Segment {
  int start;
  int end;
  bool valid;
  int cp;
  double gain;
  double lErr; //error of left segment
  double rErr; //error of right segment
  double err; //total error

  // Max-heap: higher gain has higher priority
  bool operator<(const Segment& other) const {
    return gain < other.gain;
  }
};



// ========================================================
//                Utility: miniOptHeapCpp
// ========================================================

inline Segment miniOptHeapCpp(const CostBase& costModule, int start, int end, int minLen,
                              int minSize = 1,  int jump = 1, double totalErr = -1) {

  int len = end - start;

  if(totalErr < 0){
    totalErr = costModule.eval(start, end);
  }

  if(len < minLen){

    return Segment{start, end, true, start,
                   -std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};

  } else if(len == minLen){

    int cp = start + len/2;
    double lErr = costModule.eval(start, cp);
    double rErr = costModule.eval(cp, end);
    double err = lErr + rErr;
    return Segment{start, end, true, cp,
                   totalErr - err, //gain
                   lErr,
                   rErr,
                   err};
  }

  double minErr = std::numeric_limits<double>::infinity();
  int cp;
  int tempCp = start;
  double err;
  double lErr;
  double rErr;
  double minlErr;
  double minrErr;

  auto allBkps = arma::regspace<arma::ivec>(start, jump, end); //all breakpoinbts
  allBkps = allBkps(arma::find((allBkps - start >= minSize) % (end - allBkps >= minSize)));
  //(start, End]
  for(int i = 0; i < allBkps.n_elem; i++){

    tempCp = allBkps(i);
    lErr = costModule.eval(start,tempCp);
    rErr = costModule.eval(tempCp,end);
    err = lErr + rErr;

    if(err < minErr){
      minErr = err;
      minlErr = lErr;
      minrErr = rErr;
      cp = tempCp;
    }
  }

  return Segment{start, end, true, cp,
                 totalErr - minErr, //gain
                 minlErr,
                 minrErr,
                 minErr};
}


// ========================================================
//                    binSegCppTmpl Class
// ========================================================

// Fast implementation based on heap, which involves an O(1) search instead of O(k) one

template<typename CostType>
class binSegCppTmpl {

  static_assert(std::is_base_of<CostBase, CostType>::value,
                "CostType must inherit from CostBase!");

public:
  CostType costModule;
  int minSize;
  int jump;
  int nSamples;
  IntegerVector bkpsVec;
  NumericVector costVec;
  int minLen;

  // Declare generic constructors (empty here)
  // The actual definitions will be specialized outside.

  // For VAR: constructor with (mat, pVAR, minSize, jump)
  binSegCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_);

  // For L1, L2: constructor with (mat, minSize, jump)
  binSegCppTmpl(const arma::mat& tsMat, int minSize_, int jump_);

  // For SIGMA: constructor with (mat, addSmallDiag, epsilon, minSize, jump)
  binSegCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_);


  // For LinearL2: constructor with (tsMat, covariates, intercept, minSize, jump)
  binSegCppTmpl(const arma::mat& tsMat, const arma::mat& covariates, bool intercept_, int minSize_, int jump_);

  //.fit() method
  void fit(){

    costModule.resetWarning(true); //Only output warning once - unnecessary if fit() only run once

    int nr = costModule.nr;

    const int& maxNRegimes = std::floor(nr / minSize);
    NumericVector cost(maxNRegimes);
    double initCost = costModule.eval(0,nr);  //(0, nr]
    Segment seg0 = miniOptHeapCpp(costModule, 0, nr, minLen, minSize, jump, initCost);
    IntegerVector changePoints(maxNRegimes-1);

    cost[0] = initCost;

    std::priority_queue<Segment> heap;
    heap.push(Segment{0, nr, true, seg0.cp, seg0.gain, seg0.lErr, seg0.rErr, seg0.err});

    int nRegimes = 2;
    int idx = 0;
    bool failed = false;

    do {

      Segment bestSeg = heap.top();
      heap.pop();

      changePoints[idx] = bestSeg.cp;
      idx++;
      cost[idx] = cost[idx-1] - bestSeg.gain;

      Segment leftSeg = miniOptHeapCpp(costModule, bestSeg.start, bestSeg.cp, minLen, minSize, jump, bestSeg.lErr);
      Segment rightSeg = miniOptHeapCpp(costModule, bestSeg.cp, bestSeg.end, minLen,  minSize, jump, bestSeg.rErr);
      heap.push(leftSeg);
      heap.push(rightSeg);

      if(heap.top().gain < 0){ //Negative best gain = no valid interval to split (by design)
        failed = true;
        break;
      }

      nRegimes++;

    } while(nRegimes <= maxNRegimes);

    if(not failed){ //If negative best gain not detected, add the last detected change points & corresponding cost
      changePoints[idx] = heap.top().cp;
      cost[idx] = cost[idx-1] - heap.top().gain;
    }

    bkpsVec = changePoints[Range(0,nRegimes-2)];
    costVec = cost[Range(0,nRegimes-1)];

    costModule.resetWarning(false);

  }

  //.predict() method

  IntegerVector predict(double penalty){

    if(penalty < 0){
      Rcpp::stop("`penalty must be non-negative!`");
    }

    return binSegPredCpp(bkpsVec, costVec, penalty);

  }

  //.eval() method
  double eval(int start, int end) {

    costModule.resetWarning(false);

    if(start >= end){
      Rcpp::stop("`start < end` must be true!");
    }

    if(start >= nSamples or start < 0){
      Rcpp::stop("`0 <= start < nSamples` must be true!");
    }

    if(end > nSamples or end <= 0){
      Rcpp::stop("`0 < end <= nSamples` must be true!");
    }

    return costModule.eval(start, end);

  }

};



// ========================================================
//            L1 class based on piecewise median
// ========================================================

static void L1_cwMedian() {
  // intentionally empty
}


template<>
binSegCppTmpl<Cost_L1_cwMed>::binSegCppTmpl(const arma::mat& tsMat, int minSize_, int jump_)
  : costModule(tsMat, true), minSize(minSize_), jump(jump_) {
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump; //to make sure the mid point is always of the form start + k*jump

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

}

RCPP_EXPOSED_CLASS(binSegCpp_L1_cwMed)
  RCPP_MODULE(binSegCpp_L1_cwMed_module) {
    Rcpp::class_<binSegCppTmpl<Cost_L1_cwMed>>("binSegCpp_L1_cwMed")
    .constructor<arma::mat, int, int>()       // mat, minSize, jump
    .method("fit", &binSegCppTmpl<Cost_L1_cwMed>::fit)
    .method("predict", &binSegCppTmpl<Cost_L1_cwMed>::predict)
    .method("eval", &binSegCppTmpl<Cost_L1_cwMed>::eval);
  }


// ========================================================
//                        L2 class
// ========================================================

static void L2() {
  // intentionally empty
}


template<>
binSegCppTmpl<Cost_L2>::binSegCppTmpl(const arma::mat& tsMat, int minSize_, int jump_)
  : costModule(tsMat, true), minSize(minSize_), jump(jump_) {
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump; //to make sure the mid point is always of the form start + k*jump

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

}

RCPP_EXPOSED_CLASS(binSegCpp_L2)
  RCPP_MODULE(binSegCpp_L2_module) {
    Rcpp::class_<binSegCppTmpl<Cost_L2>>("binSegCpp_L2")
    .constructor<arma::mat, int, int>()       // mat, minSize, jump
    .method("fit", &binSegCppTmpl<Cost_L2>::fit)
    .method("predict", &binSegCppTmpl<Cost_L2>::predict)
    .method("eval", &binSegCppTmpl<Cost_L2>::eval)
    .field("bkpsVec", &binSegCppTmpl<Cost_L2>::bkpsVec)
    .field("costVec", &binSegCppTmpl<Cost_L2>::costVec);
  }



// ========================================================
//                        VAR class
// ========================================================

static void VAR() {
  // intentionally empty
}


template<>
binSegCppTmpl<Cost_VAR>::binSegCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_)
  : costModule(tsMat, pVAR, true), minSize(minSize_), jump(jump_){
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump; //to make sure the mid point is always of the form start + k*jump

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

}


RCPP_EXPOSED_CLASS(binSegCpp_VAR)
  RCPP_MODULE(binSegCpp_VAR_module) {
    Rcpp::class_<binSegCppTmpl<Cost_VAR>>("binSegCpp_VAR")
    .constructor<arma::mat, int, int, int>()  // mat, pVAR, minSize, jump
    .method("fit", &binSegCppTmpl<Cost_VAR>::fit)
    .method("predict", &binSegCppTmpl<Cost_VAR>::predict)
    .method("eval", &binSegCppTmpl<Cost_VAR>::eval);
  }



// ========================================================
//                       SIGMA class
// ========================================================

static void SIGMA() {
  // intentionally empty
}


template<>
binSegCppTmpl<Cost_SIGMA>::binSegCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_)
  : costModule(tsMat, addSmallDiag, epsilon, true), minSize(minSize_), jump(jump_){
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump; //to make sure the mid point is always of the form start + k*jump

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

}

RCPP_EXPOSED_CLASS(binSegCpp_SIGMA)
  RCPP_MODULE(binSegCpp_SIGMA_module) {
    Rcpp::class_<binSegCppTmpl<Cost_SIGMA>>("binSegCpp_SIGMA")
    .constructor<arma::mat, bool, double, int, int>()  // mat, addSmallDiag, epsilon, minSize, jump
    .method("fit", &binSegCppTmpl<Cost_SIGMA>::fit)
    .method("predict", &binSegCppTmpl<Cost_SIGMA>::predict)
    .method("eval", &binSegCppTmpl<Cost_SIGMA>::eval);
  }



// ========================================================
//                     LinearL2 class
// ========================================================

static void LinearL2() {
  // intentionally empty
}


template<>
binSegCppTmpl<Cost_LinearL2>::binSegCppTmpl(const arma::mat& tsMat,  const arma::mat& covariates,
                                            bool intercept_, int minSize_, int jump_)
  : costModule(tsMat, covariates, intercept_, true), minSize(minSize_), jump(jump_){
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump; //to make sure the mid point is always of the form start + k*jump

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

}

// For LinearL2: constructor with (tsMat, covariates, intercept, minSize, jump, h)

RCPP_EXPOSED_CLASS(binSegCpp_LinearL2)
  RCPP_MODULE(binSegCpp_LinearL2_module) {
    Rcpp::class_<binSegCppTmpl<Cost_LinearL2>>("binSegCpp_LinearL2")
    .constructor<arma::mat, arma::mat, bool, int, int>()  // mat, covariates, intercept, minSize, jump
    .method("fit", &binSegCppTmpl<Cost_LinearL2>::fit)
    .method("predict", &binSegCppTmpl<Cost_LinearL2>::predict)
    .method("eval", &binSegCppTmpl<Cost_LinearL2>::eval);
  }
