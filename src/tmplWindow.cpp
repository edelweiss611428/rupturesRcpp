#include <RcppArmadillo.h>
#include <queue>
#include <limits>
#include "costs.h"

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

    int lo = std::max(0, i - order);
    int hi = std::min(nCandidates - 1, i + order);

    // Ties are broken towards the earlier index: a left neighbor that is equal-or-greater
    // disqualifies i, but a right neighbor only disqualifies i if it is strictly greater.
    // Without this asymmetry, every point in a plateau of equal gains (common on
    // piecewise-constant data) would pass as its own "local maximum", since none of them
    // is ever *strictly* beaten -- defeating the `order`/`radius` separation this function
    // exists to enforce. This way, only the plateau's leftmost point survives.
    for (int j = lo; j < i; j++) {
      if (gains[j] >= centeredVal) {
        isMax = false;
        break;
      }
    }

    if (isMax) {
      for (int j = i + 1; j <= hi; j++) {
        if (gains[j] > centeredVal) {
          isMax = false;
          break;
        }
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

  IntegerVector bkpsVec;  // local maxima in descending-gain order (the order $predict() adds them)
  NumericVector costVec;  // total cost after 0, 1, ..., nMaxima of those change-points; length nMaxima + 1

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

  // For LinearSIGMA: constructor with (tsMat, covariates, intercept, addSmallDiag, epsilon, segParams),
  // where segParams = c(minSize, jump, radius). Bundled into one vector because Rcpp Modules'
  // `.constructor<...>()` supports at most 7 template arguments, and this cost type already needs
  // 5 for its own parameters.
  windowCppTmpl(const arma::mat& tsMat, const arma::mat& covariates, bool intercept_,
                bool addSmallDiag, double epsilon, Rcpp::IntegerVector segParams);

  // For LinearL1: constructor with (tsMat, covariates, intercept, tol, maxIter, segParams), where
  // segParams = c(minSize, jump, radius) -- see the LinearSIGMA constructor above for why.
  windowCppTmpl(const arma::mat& tsMat, const arma::mat& covariates, bool intercept_,
                double tol, int maxIter, Rcpp::IntegerVector segParams);

  // For RFunc: constructor with (tsMat, costFun, paramFun, minSize, jump, radius)
  windowCppTmpl(const arma::mat& tsMat, Rcpp::Function costFun,
                Rcpp::Nullable<Rcpp::Function> paramFun, int minSize_, int jump_, int h_);

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

    // Cost history: costVec[k] is the total cost after adding the top k change-points (by gain);
    // bkpsVec[k-1] is the change-point added at that step. Same shape/semantics as binSeg's own
    // bkpsVec/costVec, computed here from the (independently scored, not re-evaluated) gains above.
    double baseCost = costModule.eval(0, nSamples);
    bkpsVec = Rcpp::wrap(sortedPeaks);
    costVec = NumericVector(nMaxima + 1);
    costVec[0] = baseCost;
    for (int kk = 1; kk <= nMaxima; kk++) {
      costVec[kk] = baseCost - cumGains[kk - 1];
    }

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
    costModule.checkSegment(start, end);
    return costModule.eval(start, end);
  }

  //.get_params() method
  Rcpp::List get_params(int start, int end) {
    costModule.resetWarning(false);
    costModule.checkSegment(start, end);
    return costModule.get_params(start, end);
  }

};



// ========================================================
//            L1 class based on piecewise median
// ========================================================



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
    .method("eval", &windowCppTmpl<Cost_L1_cwMed>::eval)
    .method("get_params", &windowCppTmpl<Cost_L1_cwMed>::get_params)
    .field("bkpsVec", &windowCppTmpl<Cost_L1_cwMed>::bkpsVec)
    .field("costVec", &windowCppTmpl<Cost_L1_cwMed>::costVec);
  }


// ========================================================
//                        L2 class
// ========================================================



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
    .method("eval", &windowCppTmpl<Cost_L2>::eval)
    .method("get_params", &windowCppTmpl<Cost_L2>::get_params)
    .field("bkpsVec", &windowCppTmpl<Cost_L2>::bkpsVec)
    .field("costVec", &windowCppTmpl<Cost_L2>::costVec);
  }



// ========================================================
//                        VAR class
// ========================================================


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
    .method("eval", &windowCppTmpl<Cost_VAR>::eval)
    .method("get_params", &windowCppTmpl<Cost_VAR>::get_params)
    .field("bkpsVec", &windowCppTmpl<Cost_VAR>::bkpsVec)
    .field("costVec", &windowCppTmpl<Cost_VAR>::costVec);
  }



// ========================================================
//                       SIGMA class
// ========================================================


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
    .method("eval", &windowCppTmpl<Cost_SIGMA>::eval)
    .method("get_params", &windowCppTmpl<Cost_SIGMA>::get_params)
    .field("bkpsVec", &windowCppTmpl<Cost_SIGMA>::bkpsVec)
    .field("costVec", &windowCppTmpl<Cost_SIGMA>::costVec);
  }



// ========================================================
//                     LinearL2 class
// ========================================================


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
    .method("eval", &windowCppTmpl<Cost_LinearL2>::eval)
    .method("get_params", &windowCppTmpl<Cost_LinearL2>::get_params)
    .field("bkpsVec", &windowCppTmpl<Cost_LinearL2>::bkpsVec)
    .field("costVec", &windowCppTmpl<Cost_LinearL2>::costVec);
  }



// ========================================================
//                   LinearSIGMA class
// ========================================================


template<>
windowCppTmpl<Cost_LinearSIGMA>::windowCppTmpl(const arma::mat& tsMat, const arma::mat& covariates,
                                                bool intercept_, bool addSmallDiag, double epsilon,
                                                Rcpp::IntegerVector segParams)
  : costModule(tsMat, covariates, intercept_, addSmallDiag, epsilon, true),
    minSize(segParams.at(0)), jump(segParams.at(1)), h(segParams.at(2)){
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

// For LinearSIGMA: constructor with (tsMat, covariates, intercept, addSmallDiag, epsilon, minSize, jump, h)

RCPP_EXPOSED_CLASS(windowCpp_LinearSIGMA)
  RCPP_MODULE(windowCpp_LinearSIGMA_module) {
    Rcpp::class_<windowCppTmpl<Cost_LinearSIGMA>>("windowCpp_LinearSIGMA")
    .constructor<arma::mat, arma::mat, bool, bool, double, Rcpp::IntegerVector>()  // mat, covariates, intercept, addSmallDiag, epsilon, c(minSize, jump, h)
    .method("fit", &windowCppTmpl<Cost_LinearSIGMA>::fit)
    .method("predict", &windowCppTmpl<Cost_LinearSIGMA>::predict)
    .method("eval", &windowCppTmpl<Cost_LinearSIGMA>::eval)
    .method("get_params", &windowCppTmpl<Cost_LinearSIGMA>::get_params)
    .field("bkpsVec", &windowCppTmpl<Cost_LinearSIGMA>::bkpsVec)
    .field("costVec", &windowCppTmpl<Cost_LinearSIGMA>::costVec);
  }



// ========================================================
//                    LinearL1 class
// ========================================================


template<>
windowCppTmpl<Cost_LinearL1>::windowCppTmpl(const arma::mat& tsMat, const arma::mat& covariates,
                                             bool intercept_, double tol, int maxIter,
                                             Rcpp::IntegerVector segParams)
  : costModule(tsMat, covariates, intercept_, tol, maxIter, true),
    minSize(segParams.at(0)), jump(segParams.at(1)), h(segParams.at(2)){
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

// For LinearL1: constructor with (tsMat, covariates, intercept, tol, maxIter, c(minSize, jump, h))

RCPP_EXPOSED_CLASS(windowCpp_LinearL1)
  RCPP_MODULE(windowCpp_LinearL1_module) {
    Rcpp::class_<windowCppTmpl<Cost_LinearL1>>("windowCpp_LinearL1")
    .constructor<arma::mat, arma::mat, bool, double, int, Rcpp::IntegerVector>()  // mat, covariates, intercept, tol, maxIter, c(minSize, jump, h)
    .method("fit", &windowCppTmpl<Cost_LinearL1>::fit)
    .method("predict", &windowCppTmpl<Cost_LinearL1>::predict)
    .method("eval", &windowCppTmpl<Cost_LinearL1>::eval)
    .method("get_params", &windowCppTmpl<Cost_LinearL1>::get_params)
    .field("bkpsVec", &windowCppTmpl<Cost_LinearL1>::bkpsVec)
    .field("costVec", &windowCppTmpl<Cost_LinearL1>::costVec);
  }



// ========================================================
//              RFunc class (user-defined cost)
// ========================================================

template<>
windowCppTmpl<Cost_RFunc>::windowCppTmpl(const arma::mat& tsMat, Rcpp::Function costFun,
                                          Rcpp::Nullable<Rcpp::Function> paramFun,
                                          int minSize_, int jump_, int h_)
  : costModule(tsMat, costFun, paramFun, true), minSize(minSize_), jump(jump_), h(h_){
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

RCPP_EXPOSED_CLASS(windowCpp_RFunc)
  RCPP_MODULE(windowCpp_RFunc_module) {
    Rcpp::class_<windowCppTmpl<Cost_RFunc>>("windowCpp_RFunc")
    .constructor<arma::mat, Rcpp::Function, Rcpp::Nullable<Rcpp::Function>, int, int, int>()  // tsMat, costFun, paramFun, minSize, jump, h
    .method("fit", &windowCppTmpl<Cost_RFunc>::fit)
    .method("predict", &windowCppTmpl<Cost_RFunc>::predict)
    .method("eval", &windowCppTmpl<Cost_RFunc>::eval)
    .method("get_params", &windowCppTmpl<Cost_RFunc>::get_params)
    .field("bkpsVec", &windowCppTmpl<Cost_RFunc>::bkpsVec)
    .field("costVec", &windowCppTmpl<Cost_RFunc>::costVec);
  }
