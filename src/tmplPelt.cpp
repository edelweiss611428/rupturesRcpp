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

    costModule.resetWarning(true);

    if (penalty < 0) {
      Rcpp::stop("`penalty` must be non-negative!");
    }

    const double inf = std::numeric_limits<double>::infinity();
    const int never = std::numeric_limits<int>::max();

    arma::vec socVec(nSamples + 1, arma::fill::value(inf));
    socVec[0] = -penalty;
    arma::ivec pathVec = arma::zeros<arma::ivec>(nSamples + 1);

    std::vector<int> ends;
    for (int k = 0; k < nSamples; k += jump) {
      if (k >= minSize) ends.push_back(k);
    }
    ends.push_back(nSamples);

    std::vector<int> admissibleBkps, pruneTime, keepB, keepP;
    std::vector<double> tmpCostVec;
    int lastAdded = -1;

    for (int end : ends) {

      int newPt = ((end - minSize) / jump) * jump;
      if (newPt != lastAdded) {
        admissibleBkps.push_back(newPt);
        pruneTime.push_back(never);
        lastAdded = newPt;
      }

      keepB.clear();
      keepP.clear();
      for (size_t j = 0; j < admissibleBkps.size(); ++j) {
        if (pruneTime[j] > end - minSize) {
          keepB.push_back(admissibleBkps[j]);
          keepP.push_back(pruneTime[j]);
        }
      }
      admissibleBkps.swap(keepB);
      pruneTime.swap(keepP);

      tmpCostVec.assign(admissibleBkps.size(), inf);
      double minSoc = inf;
      int bestBkp = -1;

      for (size_t j = 0; j < admissibleBkps.size(); ++j) {
        int t = admissibleBkps[j];
        if (socVec[t] == inf) continue;
        tmpCostVec[j] = costModule.eval(t, end);
        double total = socVec[t] + tmpCostVec[j] + penalty;
        if (total < minSoc) {
          minSoc = total;
          bestBkp = t;
        }
      }

      socVec[end] = minSoc;
      pathVec[end] = bestBkp;

      for (size_t j = 0; j < admissibleBkps.size(); ++j) {
        if (pruneTime[j] == never &&
            socVec[admissibleBkps[j]] + tmpCostVec[j] > socVec[end]) {
          pruneTime[j] = end;
        }
      }
    }

    costModule.resetWarning(false);
    return readPath(pathVec);
  };

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
