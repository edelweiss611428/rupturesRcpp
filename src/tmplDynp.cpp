#include <RcppArmadillo.h>
#include <limits>
#include <vector>
#include "costs.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// ========================================================
//                    DynpCppTmpl Class
// ========================================================

// Exact dynamic programming: for every k = 0, ..., nBkpsMax, finds the segmentation into exactly
// k change-points (k+1 segments) that globally minimises the summed cost, over the same
// (minSize, jump)-admissible grid used by PELTCppTmpl/binSegCppTmpl/windowCppTmpl. Unlike PELT
// (exact for a penalty, not a change-point count) and binSeg (greedy, not globally optimal for a
// fixed count), this is exact for a *specified* number of change-points -- the price is
// O(nBkpsMax * M^2) work and O(nBkpsMax * M) memory, where M is the number of admissible grid
// points (M ~= nSamples / jump in the worst case).

template<typename CostType>
class DynpCppTmpl {
  static_assert(std::is_base_of<CostBase, CostType>::value,
                "CostType must inherit from CostBase!");
public:
  CostType costModule;
  int minSize;
  int jump;
  int nSamples;
  int nBkpsMax;
  int minLen;

  std::vector<int> grid;  // admissible positions; grid[0] = 0, grid.back() = nSamples
  arma::mat D;            // (nBkpsMax+1) x M; D(k,i) = min cost covering [0, grid[i]) with exactly k bkps
  arma::imat P;           // backpointers; P(k,i) = predecessor grid-index achieving D(k,i)
  bool fitted_ = false;

  // Declare generic constructors (empty here)
  // The actual definitions will be specialized outside.

  // For VAR: constructor with (mat, pVAR, minSize, jump, nBkpsMax)
  DynpCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_, int nBkpsMax_);

  // For L1, L2: constructor with (mat, minSize, jump, nBkpsMax)
  DynpCppTmpl(const arma::mat& tsMat, int minSize_, int jump_, int nBkpsMax_);

  // For SIGMA: constructor with (mat, addSmallDiag, epsilon, minSize, jump, nBkpsMax)
  DynpCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_, int nBkpsMax_);

  // For LinearL2: constructor with (tsMat, covariates, intercept, minSize, jump, nBkpsMax)
  DynpCppTmpl(const arma::mat& tsMat, const arma::mat& covariates, bool intercept_, int minSize_, int jump_, int nBkpsMax_);

  // For LinearSIGMA: constructor with (tsMat, covariates, intercept, addSmallDiag, epsilon, minSize, jump, nBkpsMax)
  DynpCppTmpl(const arma::mat& tsMat, const arma::mat& covariates, bool intercept_,
              bool addSmallDiag, double epsilon, int minSize_, int jump_, int nBkpsMax_);

  // For LinearL1: constructor with (tsMat, covariates, intercept, tol, maxIter, minSize, jump, nBkpsMax)
  DynpCppTmpl(const arma::mat& tsMat, const arma::mat& covariates, bool intercept_,
              double tol, int maxIter, int minSize_, int jump_, int nBkpsMax_);

  // For RFunc: constructor with (tsMat, costFun, paramFun, minSize, jump, nBkpsMax)
  DynpCppTmpl(const arma::mat& tsMat, Rcpp::Function costFun,
              Rcpp::Nullable<Rcpp::Function> paramFun, int minSize_, int jump_, int nBkpsMax_);

  //.fit() method: builds the full D/P tables for k = 0, ..., nBkpsMax
  void fit() {

    costModule.resetWarning(true);  // Only output warning once - unnecessary if fit() only run once

    grid.clear();
    grid.push_back(0);
    for (int k = jump; k < nSamples; k += jump) {
      if (k >= minSize) grid.push_back(k);
    }
    grid.push_back(nSamples);
    const int M = static_cast<int>(grid.size());

    const double inf = std::numeric_limits<double>::infinity();

    // Precompute every admissible segment's cost once -- reused across all nBkpsMax outer
    // iterations below, instead of being recomputed per k (which would be especially wasteful
    // for costs where eval() is expensive, e.g. LinearL1's IRLS refit or Custom's R callback).
    arma::mat C(M, M, arma::fill::value(inf));
    for (int i = 1; i < M; i++) {
      for (int j = 0; j < i; j++) {
        if (grid[i] - grid[j] >= minSize) {
          C(j, i) = costModule.eval(grid[j], grid[i]);
        }
      }
    }

    D.set_size(nBkpsMax + 1, M);
    D.fill(inf);
    P.set_size(nBkpsMax + 1, M);
    P.fill(-1);

    for (int i = 1; i < M; i++) {
      D(0, i) = C(0, i);  // k = 0: the single segment (0, grid[i]]
    }

    for (int k = 1; k <= nBkpsMax; k++) {
      for (int i = 1; i < M; i++) {
        double best = inf;
        int bestJ = -1;
        for (int j = 0; j < i; j++) {
          if (!std::isfinite(D(k - 1, j)) || !std::isfinite(C(j, i))) continue;
          double cand = D(k - 1, j) + C(j, i);
          if (cand < best) {
            best = cand;
            bestJ = j;
          }
        }
        D(k, i) = best;
        P(k, i) = bestJ;
      }
    }

    fitted_ = true;
    costModule.resetWarning(false);
  }

  // Traceback from (k, grid-index i): returns the k breakpoints in ascending order, plus nSamples
  // -- same "sorted end-points, last one is nSamples" convention as PELT/binSeg/Window.
  std::vector<int> traceback(int k, int i) const {
    std::vector<int> bkps;
    while (k > 0) {
      int j = P(k, i);
      bkps.push_back(grid[j]);
      i = j;
      k--;
    }
    std::reverse(bkps.begin(), bkps.end());
    bkps.push_back(nSamples);
    return bkps;
  }

  void checkFitted() const {
    if (!fitted_) {
      Rcpp::stop("`$fit()` must be run before this method can be used!");
    }
  }

  //.predictK() method: the exact optimal segmentation using exactly `nBkps` change-points
  IntegerVector predictK(int nBkps) {
    checkFitted();
    if (nBkps < 0 || nBkps > nBkpsMax) {
      Rcpp::stop("`nBkps` must be between 0 and `nBkpsMax`!");
    }
    int lastIdx = static_cast<int>(grid.size()) - 1;
    if (!std::isfinite(D(nBkps, lastIdx))) {
      Rcpp::stop("No valid segmentation exists for this `nBkps`, given `minSize`/`jump`!");
    }
    return Rcpp::wrap(traceback(nBkps, lastIdx));
  }

  //.predictPen() method: penalised selection over k = 0, ..., nBkpsMax (same convention as
  // PELT/binSeg/Window's predict(pen)), restricted to the range the DP table was built for.
  IntegerVector predictPen(double penalty) {
    checkFitted();
    if (penalty < 0) {
      Rcpp::stop("`penalty` must be non-negative!");
    }
    int lastIdx = static_cast<int>(grid.size()) - 1;
    arma::vec costs = D.col(lastIdx);
    arma::vec penalties = arma::regspace(0, nBkpsMax) * penalty;
    arma::uword kBest = (costs + penalties).index_min();
    return Rcpp::wrap(traceback(static_cast<int>(kBest), lastIdx));
  }

  //.costPath() method: the exact minimal cost for every k = 0, ..., nBkpsMax (elbow-method data)
  NumericVector costPath() {
    checkFitted();
    int lastIdx = static_cast<int>(grid.size()) - 1;
    return Rcpp::wrap(D.col(lastIdx));
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
DynpCppTmpl<Cost_L1_cwMed>::DynpCppTmpl(const arma::mat& tsMat, int minSize_, int jump_, int nBkpsMax_)
  : costModule(tsMat, true), minSize(minSize_), jump(jump_), nBkpsMax(nBkpsMax_) {
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

  if(nBkpsMax < 0){
    Rcpp::stop("`nBkpsMax` must be non-negative!");
  }

}

RCPP_EXPOSED_CLASS(DynpCpp_L1_cwMed)
  RCPP_MODULE(DynpCpp_L1_cwMed_module) {
    Rcpp::class_<DynpCppTmpl<Cost_L1_cwMed>>("DynpCpp_L1_cwMed")
    .constructor<arma::mat, int, int, int>()       // mat, minSize, jump, nBkpsMax
    .method("fit", &DynpCppTmpl<Cost_L1_cwMed>::fit)
    .method("predictK", &DynpCppTmpl<Cost_L1_cwMed>::predictK)
    .method("predictPen", &DynpCppTmpl<Cost_L1_cwMed>::predictPen)
    .method("costPath", &DynpCppTmpl<Cost_L1_cwMed>::costPath)
    .method("eval", &DynpCppTmpl<Cost_L1_cwMed>::eval)
    .method("get_params", &DynpCppTmpl<Cost_L1_cwMed>::get_params);
  }


// ========================================================
//                        L2 class
// ========================================================

template<>
DynpCppTmpl<Cost_L2>::DynpCppTmpl(const arma::mat& tsMat, int minSize_, int jump_, int nBkpsMax_)
  : costModule(tsMat, true), minSize(minSize_), jump(jump_), nBkpsMax(nBkpsMax_) {
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump;

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

  if(nBkpsMax < 0){
    Rcpp::stop("`nBkpsMax` must be non-negative!");
  }

}

RCPP_EXPOSED_CLASS(DynpCpp_L2)
  RCPP_MODULE(DynpCpp_L2_module) {
    Rcpp::class_<DynpCppTmpl<Cost_L2>>("DynpCpp_L2")
    .constructor<arma::mat, int, int, int>()       // mat, minSize, jump, nBkpsMax
    .method("fit", &DynpCppTmpl<Cost_L2>::fit)
    .method("predictK", &DynpCppTmpl<Cost_L2>::predictK)
    .method("predictPen", &DynpCppTmpl<Cost_L2>::predictPen)
    .method("costPath", &DynpCppTmpl<Cost_L2>::costPath)
    .method("eval", &DynpCppTmpl<Cost_L2>::eval)
    .method("get_params", &DynpCppTmpl<Cost_L2>::get_params);
  }


// ========================================================
//                        VAR class
// ========================================================

template<>
DynpCppTmpl<Cost_VAR>::DynpCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_, int nBkpsMax_)
  : costModule(tsMat, pVAR, true), minSize(minSize_), jump(jump_), nBkpsMax(nBkpsMax_){
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump;

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

  if(nBkpsMax < 0){
    Rcpp::stop("`nBkpsMax` must be non-negative!");
  }

}

RCPP_EXPOSED_CLASS(DynpCpp_VAR)
  RCPP_MODULE(DynpCpp_VAR_module) {
    Rcpp::class_<DynpCppTmpl<Cost_VAR>>("DynpCpp_VAR")
    .constructor<arma::mat, int, int, int, int>()  // mat, pVAR, minSize, jump, nBkpsMax
    .method("fit", &DynpCppTmpl<Cost_VAR>::fit)
    .method("predictK", &DynpCppTmpl<Cost_VAR>::predictK)
    .method("predictPen", &DynpCppTmpl<Cost_VAR>::predictPen)
    .method("costPath", &DynpCppTmpl<Cost_VAR>::costPath)
    .method("eval", &DynpCppTmpl<Cost_VAR>::eval)
    .method("get_params", &DynpCppTmpl<Cost_VAR>::get_params);
  }


// ========================================================
//                       SIGMA class
// ========================================================

template<>
DynpCppTmpl<Cost_SIGMA>::DynpCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_, int nBkpsMax_)
  : costModule(tsMat, addSmallDiag, epsilon, true), minSize(minSize_), jump(jump_), nBkpsMax(nBkpsMax_){
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump;

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

  if(nBkpsMax < 0){
    Rcpp::stop("`nBkpsMax` must be non-negative!");
  }

}

RCPP_EXPOSED_CLASS(DynpCpp_SIGMA)
  RCPP_MODULE(DynpCpp_SIGMA_module) {
    Rcpp::class_<DynpCppTmpl<Cost_SIGMA>>("DynpCpp_SIGMA")
    .constructor<arma::mat, bool, double, int, int, int>()  // mat, addSmallDiag, epsilon, minSize, jump, nBkpsMax
    .method("fit", &DynpCppTmpl<Cost_SIGMA>::fit)
    .method("predictK", &DynpCppTmpl<Cost_SIGMA>::predictK)
    .method("predictPen", &DynpCppTmpl<Cost_SIGMA>::predictPen)
    .method("costPath", &DynpCppTmpl<Cost_SIGMA>::costPath)
    .method("eval", &DynpCppTmpl<Cost_SIGMA>::eval)
    .method("get_params", &DynpCppTmpl<Cost_SIGMA>::get_params);
  }


// ========================================================
//                     LinearL2 class
// ========================================================

template<>
DynpCppTmpl<Cost_LinearL2>::DynpCppTmpl(const arma::mat& tsMat, const arma::mat& covariates,
                                         bool intercept_, int minSize_, int jump_, int nBkpsMax_)
  : costModule(tsMat, covariates, intercept_, true), minSize(minSize_), jump(jump_), nBkpsMax(nBkpsMax_){
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump;

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

  if(nBkpsMax < 0){
    Rcpp::stop("`nBkpsMax` must be non-negative!");
  }

}

RCPP_EXPOSED_CLASS(DynpCpp_LinearL2)
  RCPP_MODULE(DynpCpp_LinearL2_module) {
    Rcpp::class_<DynpCppTmpl<Cost_LinearL2>>("DynpCpp_LinearL2")
    .constructor<arma::mat, arma::mat, bool, int, int, int>()  // mat, covariates, intercept, minSize, jump, nBkpsMax
    .method("fit", &DynpCppTmpl<Cost_LinearL2>::fit)
    .method("predictK", &DynpCppTmpl<Cost_LinearL2>::predictK)
    .method("predictPen", &DynpCppTmpl<Cost_LinearL2>::predictPen)
    .method("costPath", &DynpCppTmpl<Cost_LinearL2>::costPath)
    .method("eval", &DynpCppTmpl<Cost_LinearL2>::eval)
    .method("get_params", &DynpCppTmpl<Cost_LinearL2>::get_params);
  }


// ========================================================
//                   LinearSIGMA class
// ========================================================

template<>
DynpCppTmpl<Cost_LinearSIGMA>::DynpCppTmpl(const arma::mat& tsMat, const arma::mat& covariates,
                                            bool intercept_, bool addSmallDiag, double epsilon,
                                            int minSize_, int jump_, int nBkpsMax_)
  : costModule(tsMat, covariates, intercept_, addSmallDiag, epsilon, true),
    minSize(minSize_), jump(jump_), nBkpsMax(nBkpsMax_){
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump;

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

  if(nBkpsMax < 0){
    Rcpp::stop("`nBkpsMax` must be non-negative!");
  }

}

RCPP_EXPOSED_CLASS(DynpCpp_LinearSIGMA)
  RCPP_MODULE(DynpCpp_LinearSIGMA_module) {
    Rcpp::class_<DynpCppTmpl<Cost_LinearSIGMA>>("DynpCpp_LinearSIGMA")
    .constructor<arma::mat, arma::mat, bool, bool, double, int, int, int>()  // mat, covariates, intercept, addSmallDiag, epsilon, minSize, jump, nBkpsMax
    .method("fit", &DynpCppTmpl<Cost_LinearSIGMA>::fit)
    .method("predictK", &DynpCppTmpl<Cost_LinearSIGMA>::predictK)
    .method("predictPen", &DynpCppTmpl<Cost_LinearSIGMA>::predictPen)
    .method("costPath", &DynpCppTmpl<Cost_LinearSIGMA>::costPath)
    .method("eval", &DynpCppTmpl<Cost_LinearSIGMA>::eval)
    .method("get_params", &DynpCppTmpl<Cost_LinearSIGMA>::get_params);
  }


// ========================================================
//                    LinearL1 class
// ========================================================

template<>
DynpCppTmpl<Cost_LinearL1>::DynpCppTmpl(const arma::mat& tsMat, const arma::mat& covariates,
                                         bool intercept_, double tol, int maxIter,
                                         int minSize_, int jump_, int nBkpsMax_)
  : costModule(tsMat, covariates, intercept_, tol, maxIter, true),
    minSize(minSize_), jump(jump_), nBkpsMax(nBkpsMax_){
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump;

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

  if(nBkpsMax < 0){
    Rcpp::stop("`nBkpsMax` must be non-negative!");
  }

}

RCPP_EXPOSED_CLASS(DynpCpp_LinearL1)
  RCPP_MODULE(DynpCpp_LinearL1_module) {
    Rcpp::class_<DynpCppTmpl<Cost_LinearL1>>("DynpCpp_LinearL1")
    .constructor<arma::mat, arma::mat, bool, double, int, int, int, int>()  // mat, covariates, intercept, tol, maxIter, minSize, jump, nBkpsMax
    .method("fit", &DynpCppTmpl<Cost_LinearL1>::fit)
    .method("predictK", &DynpCppTmpl<Cost_LinearL1>::predictK)
    .method("predictPen", &DynpCppTmpl<Cost_LinearL1>::predictPen)
    .method("costPath", &DynpCppTmpl<Cost_LinearL1>::costPath)
    .method("eval", &DynpCppTmpl<Cost_LinearL1>::eval)
    .method("get_params", &DynpCppTmpl<Cost_LinearL1>::get_params);
  }



// ========================================================
//              RFunc class (user-defined cost)
// ========================================================

template<>
DynpCppTmpl<Cost_RFunc>::DynpCppTmpl(const arma::mat& tsMat, Rcpp::Function costFun,
                                      Rcpp::Nullable<Rcpp::Function> paramFun,
                                      int minSize_, int jump_, int nBkpsMax_)
  : costModule(tsMat, costFun, paramFun, true), minSize(minSize_), jump(jump_), nBkpsMax(nBkpsMax_) {
  nSamples = costModule.nr;

  if(minSize < 1){
    Rcpp::stop("`minSize` must be at least 1!");
  }

  if(jump < 1){
    Rcpp::stop("`jump` must be at least 1!");
  }

  int k = static_cast<int>(std::ceil(static_cast<double>(minSize) / jump));
  minLen = 2 * k * jump;

  if(nSamples < minLen){
    Rcpp::stop("Number of observations must be at least `2*jump*ceiling(minSize/jump)`!");
  }

  if(nSamples <= jump){
    Rcpp::stop("Number of observations must be larger than `jump`!");
  }

  if(nBkpsMax < 0){
    Rcpp::stop("`nBkpsMax` must be non-negative!");
  }

}

RCPP_EXPOSED_CLASS(DynpCpp_RFunc)
  RCPP_MODULE(DynpCpp_RFunc_module) {
    Rcpp::class_<DynpCppTmpl<Cost_RFunc>>("DynpCpp_RFunc")
    .constructor<arma::mat, Rcpp::Function, Rcpp::Nullable<Rcpp::Function>, int, int, int>()  // tsMat, costFun, paramFun, minSize, jump, nBkpsMax
    .method("fit", &DynpCppTmpl<Cost_RFunc>::fit)
    .method("predictK", &DynpCppTmpl<Cost_RFunc>::predictK)
    .method("predictPen", &DynpCppTmpl<Cost_RFunc>::predictPen)
    .method("costPath", &DynpCppTmpl<Cost_RFunc>::costPath)
    .method("eval", &DynpCppTmpl<Cost_RFunc>::eval)
    .method("get_params", &DynpCppTmpl<Cost_RFunc>::get_params);
  }
