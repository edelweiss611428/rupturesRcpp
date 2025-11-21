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
//                 Utility: findLocalMaxima
// ========================================================


inline arma::uvec findLocalMaxima(const arma::vec& gains,
                           const arma::uvec& candidates,
                           int order) {
  int nCandidates = gains.n_elem;
  std::vector<unsigned int> peaks;

  for (int i = 0; i < nCandidates; ++i) {
    bool isMax = true;
    double centeredVal = gains[i];

    // Compare with neighbors in [i - order, i + order]
    for (int j = std::max(0, i - order); j <= std::min(nCandidates - 1, i + order); j++) {
      if (j == i) continue;
      if (gains[j] > centeredVal) {
        isMax = false;
        break;
      }
    }

    if (isMax) {
      peaks.push_back(candidates[i]);
    }
  }

  // Convert to arma::uvec
  return arma::uvec(peaks);
}

// ========================================================
//                    windowCppTmpl class
// ========================================================


template<typename CostType>
class windowCppTmpl {

  static_assert(std::is_base_of<CostBase, CostType>::value,
                "CostType must inherit from CostBase!");

public:
  CostType costModule;
  int minSize;
  int jump;
  int nSamples;
  int h;
  int minLen;

  arma::uvec candidates;
  arma::vec gains;
  int nMaxima;

  arma::uvec sortedPeaks;
  arma::vec cumGains;

  // Declare generic constructors (empty here)
  // The actual definitions will be specialized outside.

  // For VAR: constructor with (tsMat, pVAR, minSize, jump, radius)
  windowCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_, int h_);

  // For L1, L2: constructor with (tsMat, minSize, jump, radius)
  windowCppTmpl(const arma::mat& tsMat, int minSize_, int jump_, int h_);

  // For SIGMA: constructor with (tsMat, addSmallDiag, epsilon, minSize, jump, radius)
  windowCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_, int h_);

  // For LinearL2: constructor with (tsMat, covariates, intercept, minSize, jump, radius)
  windowCppTmpl(const arma::mat& tsMat, const arma::mat& covariates, bool intercept_, int minSize_, int jump_, int h_);

  //.fit() method
  void fit(){

    costModule.resetWarning(true); //Only output warning once - unnecessary if fit() only run once

    int nCandidates = (nSamples - 2 * h)/jump + 1; //Integer division
    candidates.resize(nCandidates);
    gains.resize(nCandidates);

    double err;
    double lErr; //left error
    double rErr;

    for (int i = 0; i < nCandidates; ++i) {
      int center = h + i * jump;
      candidates[i] = center;

      err = costModule.eval(center - h, center + h);
      lErr  = costModule.eval(center - h, center);
      rErr = costModule.eval(center, center + h);

      gains[i] = err - lErr - rErr;
    }

    int k = std::max(std::max(h*2, 2 * minSize) / (2 * jump), 1);
    arma::uvec localMaxima = findLocalMaxima(gains, candidates, k);
    nMaxima = localMaxima.n_elem;

    // No local maximum
    if (nMaxima == 0) {
      Rcpp::warning("There is no valid breakpoint!");
    }

    // Extract local maxima

    arma::uvec validPeaks(nMaxima);
    arma::vec validGains(nMaxima);

    for (int i = 0; i < nMaxima; ++i) {
      arma::uvec match = arma::find(candidates == localMaxima[i]);
      validPeaks[i] = match[0];
    }

    for (int i = 0; i < nMaxima; ++i) {
      validGains[i] = gains[validPeaks[i]];
    }

    // Sorted peaks/gains
    arma::uvec sortedIdx = arma::sort_index(validGains, "descend");
    sortedPeaks = localMaxima.elem(sortedIdx);
    arma::vec sortedGains = validGains.elem(sortedIdx);

    // cumulative gains
    cumGains = arma::cumsum(sortedGains);

    costModule.resetWarning(false);

  }

  //.predict() method

  Rcpp::IntegerVector predict(double penalty) {

    if (penalty < 0) {
      Rcpp::stop("Penalty must be non-negative!");
    }

    // No local maximum
    if (nMaxima == 0) {
      return Rcpp::IntegerVector::create(nSamples);
    }

    // penalties vector
    arma::vec penalties = arma::regspace(1, nMaxima) * penalty;

    // Penalised cumulative gains

    arma::vec penCumGains = cumGains - penalties;
    arma::uword bestK = penCumGains.index_max();

    std::vector<int> selectedBkps;
    for (arma::uword i = 0; i <= bestK; ++i) {
      selectedBkps.push_back(sortedPeaks[i]);
    }

    selectedBkps.push_back(nSamples);

    // Sort ascending
    std::sort(selectedBkps.begin(), selectedBkps.end());

    return Rcpp::wrap(selectedBkps);

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
windowCppTmpl<Cost_L1_cwMed>::windowCppTmpl(const arma::mat& tsMat, int minSize_, int jump_, int h_)
  : costModule(tsMat, true), minSize(minSize_), jump(jump_), h(h_) {
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

  if(nSamples <= 2*h){
    Rcpp::stop("Number of observations must be larger than `2*radius`!");
  }

  if(h < 1){
    Rcpp::stop("Radius must be at least 1!");
  }

  //Removed

  // if(2*h <= minSize){
  //   Rcpp::warning("Diameter should be at least `minSize`");
  // }

}

RCPP_EXPOSED_CLASS(windowCpp_L1_cwMed)
  RCPP_MODULE(windowCpp_L1_cwMed_module) {
    Rcpp::class_<windowCppTmpl<Cost_L1_cwMed>>("windowCpp_L1_cwMed")
    .constructor<arma::mat, int, int, int>()       // mat, minSize, jump, h
    .method("fit", &windowCppTmpl<Cost_L1_cwMed>::fit)
    .method("predict", &windowCppTmpl<Cost_L1_cwMed>::predict)
    .method("eval", &windowCppTmpl<Cost_L1_cwMed>::eval);
  }


// ========================================================
//                        L2 class
// ========================================================

static void L2() {
  // intentionally empty
}


template<>
windowCppTmpl<Cost_L2>::windowCppTmpl(const arma::mat& tsMat, int minSize_, int jump_, int h_)
  : costModule(tsMat, true), minSize(minSize_), jump(jump_), h(h_) {
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

  if(nSamples <= 2*h){
    Rcpp::stop("Number of observations must be larger than `2*radius`!");
  }

  if(h < 1){
    Rcpp::stop("Radius must be at least 1!");
  }

}

RCPP_EXPOSED_CLASS(windowCpp_L2)
  RCPP_MODULE(windowCpp_L2_module) {
    Rcpp::class_<windowCppTmpl<Cost_L2>>("windowCpp_L2")
    .constructor<arma::mat, int, int, int>()       // mat, minSize, jump, h
    .method("fit", &windowCppTmpl<Cost_L2>::fit)
    .method("predict", &windowCppTmpl<Cost_L2>::predict)
    .method("eval", &windowCppTmpl<Cost_L2>::eval);
  }



// ========================================================
//                        VAR class
// ========================================================

static void VAR() {
  // intentionally empty
}


template<>
windowCppTmpl<Cost_VAR>::windowCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_, int h_)
  : costModule(tsMat, pVAR, true), minSize(minSize_), jump(jump_), h(h_){
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

  if(nSamples <= 2*h){
    Rcpp::stop("Number of observations must be larger than `2*radius`!");
  }

  if(h < 1){
    Rcpp::stop("Radius must be at least 1!");
  }


}


RCPP_EXPOSED_CLASS(windowCpp_VAR)
  RCPP_MODULE(windowCpp_VAR_module) {
    Rcpp::class_<windowCppTmpl<Cost_VAR>>("windowCpp_VAR")
    .constructor<arma::mat, int, int, int, int>()  // mat, pVAR, minSize, jump, h
    .method("fit", &windowCppTmpl<Cost_VAR>::fit)
    .method("predict", &windowCppTmpl<Cost_VAR>::predict)
    .method("eval", &windowCppTmpl<Cost_VAR>::eval);
  }



// ========================================================
//                       SIGMA class
// ========================================================

static void SIGMA() {
  // intentionally empty
}


template<>
windowCppTmpl<Cost_SIGMA>::windowCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_,
                                         int h_)
  : costModule(tsMat, addSmallDiag, epsilon, true), minSize(minSize_), jump(jump_), h(h_){
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

  if(nSamples <= 2*h){
    Rcpp::stop("Number of observations must be larger than `2*radius`!");
  }

  if(h < 1){
    Rcpp::stop("Radius must be at least 1!");
  }


}

RCPP_EXPOSED_CLASS(windowCpp_SIGMA)
  RCPP_MODULE(windowCpp_SIGMA_module) {
    Rcpp::class_<windowCppTmpl<Cost_SIGMA>>("windowCpp_SIGMA")
    .constructor<arma::mat, bool, double, int, int, int>()  // mat, addSmallDiag, epsilon, minSize, jump, h
    .method("fit", &windowCppTmpl<Cost_SIGMA>::fit)
    .method("predict", &windowCppTmpl<Cost_SIGMA>::predict)
    .method("eval", &windowCppTmpl<Cost_SIGMA>::eval);
  }



// ========================================================
//                     LinearL2 class
// ========================================================

static void LinearL2() {
  // intentionally empty
}


template<>
windowCppTmpl<Cost_LinearL2>::windowCppTmpl(const arma::mat& tsMat,  const arma::mat& covariates,
                                            bool intercept_, int minSize_, int jump_, int h_)
  : costModule(tsMat, covariates, intercept_, true), minSize(minSize_), jump(jump_), h(h_){
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

  if(nSamples <= 2*h){
    Rcpp::stop("Number of observations must be larger than `2*radius`!");
  }

  if(h < 1){
    Rcpp::stop("Radius must be at least 1!");
  }


}

// For LinearL2: constructor with (tsMat, covariates, intercept, minSize, jump, h)

RCPP_EXPOSED_CLASS(windowCpp_LinearL2)
  RCPP_MODULE(windowCpp_LinearL2_module) {
    Rcpp::class_<windowCppTmpl<Cost_LinearL2>>("windowCpp_LinearL2")
    .constructor<arma::mat, arma::mat, bool, int, int, int>()  // mat, covariates, intercept, minSize, jump, h
    .method("fit", &windowCppTmpl<Cost_LinearL2>::fit)
    .method("predict", &windowCppTmpl<Cost_LinearL2>::predict)
    .method("eval", &windowCppTmpl<Cost_LinearL2>::eval);
  }
