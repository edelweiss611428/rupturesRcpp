#include <RcppArmadillo.h>
#include <queue>
#include <limits>
#include "VAR.h"
#include "L2.h"
#include "SIGMA.h"
#include "L1_cwMed.h"
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
  arma::uvec candidates;
  arma::vec gains;
  arma::uvec localMaxima;
  int h;

  // Declare generic constructors (empty here)
  // The actual definitions will be specialized outside.

  // For VAR: constructor with (mat, pVAR, minSize, jump)
  windowCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_, int h_);

  // For L1, L2: constructor with (mat, minSize, jump)
  windowCppTmpl(const arma::mat& tsMat, int minSize_, int jump_, int h_);

  // For SIGMA: constructor with (mat, addSmallDiag, epsilon, minSize, jump)
  windowCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_, int h_);

  //.fit() method
  void fit(){

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
    localMaxima = findLocalMaxima(gains, candidates, k);

  }

  //.predict() method

  Rcpp::IntegerVector predict(double penalty) {

    if (penalty < 0) {
      Rcpp::stop("Penalty must be non-negative!");
    }

    int nMaxima = localMaxima.n_elem;

    // Special case: no local maxima → return only the final index
    if (nMaxima == 0) {
      return Rcpp::IntegerVector::create(nSamples);
    }

    // 1. Map localMaxima positions to their indices in `candidates`
    arma::uvec validPeaks(nMaxima);
    for (int i = 0; i < nMaxima; ++i) {
      arma::uvec match = arma::find(candidates == localMaxima[i]);
      validPeaks[i] = match[0];
    }

    // 2. Extract gains at these positions
    arma::vec validGains(nMaxima);
    for (int i = 0; i < nMaxima; ++i) {
      validGains[i] = gains[validPeaks[i]];
    }

    // 3. Sort by descending gain
    arma::uvec sortedIdx = arma::sort_index(validGains, "descend");
    arma::uvec sortedPeaks = localMaxima.elem(sortedIdx);
    arma::vec sortedGains = validGains.elem(sortedIdx);

    // 4. Compute cumulative sums of gains
    arma::vec cumGains = arma::cumsum(sortedGains);

    // 5. Create penalty vector: 1, 2, ..., nMaxima
    arma::vec penalties = arma::regspace(1, nMaxima) * penalty;

    // Compute penalized gains vector

    arma::vec penGains = cumGains - penalties;
    arma::uword bestK;
    penGains.max(bestK);

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
      Rcpp::stop("`0 <= end <= nSamples` must be true!");
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations < 2*minSize!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations <= jump!");
  }

}

RCPP_EXPOSED_CLASS(windowCpp_L1_cwMed)
  RCPP_MODULE(windowCpp_L1_cwMed_module) {
    Rcpp::class_<windowCppTmpl<Cost_L1_cwMed>>("windowCpp_L1_cwMed")
    .constructor<arma::mat, int, int, int>()       // mat, minSize, jump
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations < 2*minSize!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations <= jump!");
  }

}

RCPP_EXPOSED_CLASS(windowCpp_L2)
  RCPP_MODULE(windowCpp_L2_module) {
    Rcpp::class_<windowCppTmpl<Cost_L2>>("windowCpp_L2")
    .constructor<arma::mat, int, int, int>()       // mat, minSize, jump
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations < 2*minSize!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations <= jump!");
  }

}


RCPP_EXPOSED_CLASS(windowCpp_VAR)
  RCPP_MODULE(windowCpp_VAR_module) {
    Rcpp::class_<windowCppTmpl<Cost_VAR>>("windowCpp_VAR")
    .constructor<arma::mat, int, int, int, int>()  // mat, pVAR, minSize, jump
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations < 2*minSize!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations <= jump!");
  }

}

RCPP_EXPOSED_CLASS(windowCpp_SIGMA)
  RCPP_MODULE(windowCpp_SIGMA_module) {
    Rcpp::class_<windowCppTmpl<Cost_SIGMA>>("windowCpp_SIGMA")
    .constructor<arma::mat, bool, double, int, int, int>()  // mat, addSmallDiag, epsilon, minSize, jump
    .method("fit", &windowCppTmpl<Cost_SIGMA>::fit)
    .method("predict", &windowCppTmpl<Cost_SIGMA>::predict)
    .method("eval", &windowCppTmpl<Cost_SIGMA>::eval);
  }
