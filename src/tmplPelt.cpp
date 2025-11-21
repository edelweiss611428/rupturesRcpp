#include <RcppArmadillo.h>
#include "VAR.h"
#include "L2.h"
#include "SIGMA.h"
#include "L1_cwMed.h"
#include "LinearL2.h"
#include "baseClass.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// ========================================================
//                     Utility: readPath
// ========================================================

inline std::vector<int> readPath(const arma::ivec& pathVec)
{
  const int nSamples = static_cast<int>(pathVec.n_elem) - 1;
  std::vector<int> bkps;
  int end = nSamples;

  while (end > 0)
  {
    bkps.push_back(end);
    end = static_cast<int>(pathVec[end]);
  }

  std::reverse(bkps.begin(), bkps.end());
  return bkps;
}


// ========================================================
//                     PELTCppTmpl Class
// ========================================================


template<typename CostType>
class PELTCppTmpl {
  static_assert(std::is_base_of<CostBase, CostType>::value,
                "CostType must inherit from CostBase!");
public:
  CostType costModule;
  int minSize;
  int jump;
  int nSamples;
  int minLen;

  // Declare generic constructors (empty here)
  // The actual definitions will be specialized outside.

  // For VAR: constructor with (mat, pVAR, minSize, jump)
  PELTCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_);

  // For L2: constructor with (mat, minSize, jump)
  PELTCppTmpl(const arma::mat& tsMat, int minSize_, int jump_);

  // For SIGMA: constructor with (mat, addSmallDiag, epsilon, minSize, jump)
  PELTCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_);


  // For LinearL2: constructor with (tsMat, covariates, intercept, minSize, jump)
  PELTCppTmpl(const arma::mat& tsMat, const arma::mat& covariates, bool intercept_, int minSize_, int jump_);


  // predict() method: Perform PELT segmentation
  std::vector<int> predict(double penalty) {

    costModule.resetWarning(true); // Only output warning once

    if (penalty < 0) {
      Rcpp::stop("`penalty` must be non-negative!");
    }

    // Initialization

    arma::vec socVec(nSamples + 1);     // Total cost up to each point

    socVec[0] = -penalty;

    for (int k = 1; k < minSize; k++){
      socVec[k] = std::numeric_limits<double>::infinity();
    }

    arma::ivec pathVec = arma::zeros<arma::ivec>(nSamples + 1);  // Backpointers
    std::vector<int> admissibleBkps;             // Admissible previous breakpoints
    std::vector<double> tmpCostVec;   // Temporary cost storage


    // Build initial admissible set: all multiples of jump >= minSize
    std::vector<int> all_ends;

    for (int k = 0; k < nSamples; k += jump) {
      if (k >= minSize) all_ends.push_back(k);
    }
    all_ends.push_back(nSamples);

    for (int i = 0; i < static_cast<int>(all_ends.size()); ++i) {
      int end = all_ends[i];

      // Add new admissible point
      int new_adm_pt = static_cast<int>(std::floor((end - minSize) / double(jump))) * jump;
      admissibleBkps.push_back(new_adm_pt);

      tmpCostVec.resize(admissibleBkps.size());
      double minSoc = std::numeric_limits<double>::infinity();
      int bestBkp = -1;

      for (int j = 0; j < admissibleBkps.size(); ++j) {

        int lastBkp = admissibleBkps[j];
        // # nocov start
        if (end - lastBkp < minSize) continue;
        // # nocov end
        double segCost = costModule.eval(lastBkp, end);
        double totalCost = socVec[lastBkp] + segCost + penalty;
        tmpCostVec[j] = segCost;

        if (totalCost < minSoc) {
          minSoc = totalCost;
          bestBkp = lastBkp;
        }
      }

      socVec[end] = minSoc;
      pathVec[end] = bestBkp;

      // Prune admissible set

      std::vector<int> pruned_admissibleBkps;
      for (size_t j = 0; j < admissibleBkps.size(); ++j) {
        int lastBkp = admissibleBkps[j];
        if (socVec[lastBkp] + tmpCostVec[j] <= socVec[end]){
          pruned_admissibleBkps.push_back(lastBkp);
        }
      }
      admissibleBkps = std::move(pruned_admissibleBkps);
      tmpCostVec.resize(admissibleBkps.size());
    }

    costModule.resetWarning(false);
    return readPath(pathVec);
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
//         L1 class based on coordinate-wise median
// ========================================================

static void L1_cwMedian() {
  // intentionally empty
}


template<>
PELTCppTmpl<Cost_L1_cwMed>::PELTCppTmpl(const arma::mat& tsMat, int minSize_, int jump_)
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

RCPP_EXPOSED_CLASS(PELTCpp_L1_cwMed)
  RCPP_MODULE(PELTCpp_L1_cwMed_module) {
    Rcpp::class_<PELTCppTmpl<Cost_L1_cwMed>>("PELTCpp_L1_cwMed")
    .constructor<arma::mat, int, int>()       // tsMat, minSize, jump
    .method("predict", &PELTCppTmpl<Cost_L1_cwMed>::predict)
    .method("eval", &PELTCppTmpl<Cost_L1_cwMed>::eval);
  }


// ========================================================
//                        L2 class
// ========================================================

static void L2() {
  // intentionally empty
}


template<>
PELTCppTmpl<Cost_L2>::PELTCppTmpl(const arma::mat& tsMat, int minSize_, int jump_)
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

RCPP_EXPOSED_CLASS(PELTCpp_L2)
  RCPP_MODULE(PELTCpp_L2_module) {
    Rcpp::class_<PELTCppTmpl<Cost_L2>>("PELTCpp_L2")
    .constructor<arma::mat, int, int>()       // tsMat, minSize, jump
    .method("predict", &PELTCppTmpl<Cost_L2>::predict)
    .method("eval", &PELTCppTmpl<Cost_L2>::eval);
  }



// ========================================================
//                        VAR class
// ========================================================

static void VAR() {
  // intentionally empty
}


template<>
PELTCppTmpl<Cost_VAR>::PELTCppTmpl(const arma::mat& tsMat, int pVAR, int minSize_, int jump_)
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


RCPP_EXPOSED_CLASS(PELTCpp_VAR)
  RCPP_MODULE(PELTCpp_VAR_module) {
    Rcpp::class_<PELTCppTmpl<Cost_VAR>>("PELTCpp_VAR")
    .constructor<arma::mat, int, int, int>()  // tsMat, pVAR, minSize, jump
    .method("predict", &PELTCppTmpl<Cost_VAR>::predict)
    .method("eval", &PELTCppTmpl<Cost_VAR>::eval);
  }



// ========================================================
//                       SIGMA class
// ========================================================

static void SIGMA() {
  // intentionally empty
}


template<>
PELTCppTmpl<Cost_SIGMA>::PELTCppTmpl(const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_)
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

RCPP_EXPOSED_CLASS(PELTCpp_SIGMA)
  RCPP_MODULE(PELTCpp_SIGMA_module) {
    Rcpp::class_<PELTCppTmpl<Cost_SIGMA>>("PELTCpp_SIGMA")
    .constructor<arma::mat, bool, double, int, int>()  // tsMat, addSmallDiag, epsilon, minSize, jump
    .method("predict", &PELTCppTmpl<Cost_SIGMA>::predict)
    .method("eval", &PELTCppTmpl<Cost_SIGMA>::eval);
  }



// ========================================================
//                     LinearL2 class
// ========================================================

static void LinearL2() {
  // intentionally empty
}


template<>
PELTCppTmpl<Cost_LinearL2>::PELTCppTmpl(const arma::mat& tsMat,  const arma::mat& covariates,
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

RCPP_EXPOSED_CLASS(PELTCpp_LinearL2)
  RCPP_MODULE(PELTCpp_LinearL2_module) {
    Rcpp::class_<PELTCppTmpl<Cost_LinearL2>>("PELTCpp_LinearL2")
    .constructor<arma::mat, arma::mat, bool, int, int>()  // mat, covariates, intercept, minSize, jump
    .method("predict", &PELTCppTmpl<Cost_LinearL2>::predict)
    .method("eval", &PELTCppTmpl<Cost_LinearL2>::eval);
  }
