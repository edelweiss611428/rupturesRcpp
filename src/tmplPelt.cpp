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


  // The predict() method as before...
  std::vector<int> predict(double penalty) {

    costModule.resetWarning(true); //Only output warning once

    if(penalty < 0){
      Rcpp::stop("`penalty must be non-negative!`");
    }

    // Initialization
    arma::vec socVec = arma::zeros(nSamples + 1);
    arma::ivec pathVec = arma::zeros<arma::ivec>(nSamples + 1);
    arma::ivec admissibleBkps = arma::zeros<arma::ivec>(nSamples + 1);

    arma::vec tmpCostVec = arma::zeros(nSamples + 1);
    socVec[0] = 0.0;
    pathVec[0] = 0;
    admissibleBkps[0] = 0;
    int nAdmissibleBkps = 1;

    auto all_ends = arma::regspace<arma::ivec>(0, jump, nSamples - 1);
    all_ends = all_ends(arma::find(all_ends >= 2 * minSize));
    all_ends.insert_rows(all_ends.n_rows, 1);
    all_ends(all_ends.n_rows - 1) = nSamples;

    for (arma::uword i = 0; i < all_ends.n_elem; ++i)
    {
      const arma::uword end = all_ends[i];
      double minSoc = std::numeric_limits<double>::infinity();
      int bestBkp = 0;

      for (arma::uword kLastBkp = 0; kLastBkp < nAdmissibleBkps; ++kLastBkp)
      {
        const int lastBkp = admissibleBkps[kLastBkp];
        // LCOV_EXCL_START
        if(end - lastBkp < minSize){ //Error message! However, by design, this should not happen.
          Rcpp::Rcout << "end - lastBkp < minSize!"<< std::endl;
          continue;
        }
        // LCOV_EXCL_STOP
        const double currentCost = costModule.eval(lastBkp, end);

        const double currentSoc = socVec[lastBkp] + currentCost + penalty;
        tmpCostVec[kLastBkp] = currentCost;

        if (currentSoc < minSoc)
        {
          minSoc = currentSoc;
          bestBkp = lastBkp;
        }
      }

      socVec[end] = minSoc;
      pathVec[end] = bestBkp;

      // Pruning
      int nAdmissibleBkpsNew = 0;
      for (int kLastBkp = 0; kLastBkp < nAdmissibleBkps; ++kLastBkp)
      {
        const int lastBkp = admissibleBkps[kLastBkp];
        const double currentCost = tmpCostVec[kLastBkp];
        if (socVec[lastBkp] + currentCost <= minSoc) {
          admissibleBkps[nAdmissibleBkpsNew++] = lastBkp;
        }
      }

      // Add "end" to the set of admissible change-points
      admissibleBkps[nAdmissibleBkpsNew++] = static_cast<int>(std::floor(
        static_cast<double>(end - minSize + 1) / jump) * jump);
      nAdmissibleBkps = nAdmissibleBkpsNew;
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations must be at least `2*minSize`!");
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations must be at least `2*minSize`!");
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations must be at least `2*minSize`!");
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations must be at least `2*minSize`!");
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

  if(nSamples < 2*minSize){
    Rcpp::stop("Number of observations must be at least `2*minSize`!");
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
