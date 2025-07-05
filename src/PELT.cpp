#include <RcppArmadillo.h>
#include "L2.h"
#include "SIGMA.h"
#include "VAR.h"
#include "baseClass.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// Function to read the path from the path vector
std::vector<int> readPath(const arma::ivec& pathVec)
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


// Function to implement the PELT algorithm
// [[Rcpp::export]]
std::vector<int> PELTCpp(const arma::mat& tsMat, const double penalty, const int minSize, const int jump,
                        const Rcpp::List& costFuncObj = R_NilValue) {
  std::string costFunc;

  if (Rf_isNull(costFuncObj)) {
    stop("costFuncObj list must be provided!");

  } else{
    if(costFuncObj.containsElementNamed("costFunc")){

      if(!Rf_isString(costFuncObj["costFunc"])){
        stop("costFunc must be a single character!");

      } else{
        costFunc = as<std::string>(costFuncObj["costFunc"]);

      }
    } else{
      stop("costFunc is missing!");

    }
  }

  //Prevent memory leak
  std::unique_ptr<CostBase> Xnewptr;

  if(costFunc == "VAR"){

    int pVAR;

    if(costFuncObj.containsElementNamed("pVAR")){

      if(!Rf_isInteger(costFuncObj["pVAR"])){
        stop("pVAR must be a single positive integer!");

      } else{
        pVAR = as<int>(costFuncObj["pVAR"]);

        if(pVAR <= 0){
          stop("pVAR must be a single positive integer!");
        }

        Xnewptr = std::make_unique<Cost_VAR>(tsMat, pVAR, true);

      }
    } else{
      stop("pVAR is missing!");

    }

  } else if(costFunc == "SIGMA"){

    bool addSmallDiag;
    double epsilon;

    if(costFuncObj.containsElementNamed("addSmallDiag") and costFuncObj.containsElementNamed("epsilon")){

      if(!Rf_isLogical(costFuncObj["addSmallDiag"])){
        stop("addSmallDiag must be a single boolean value!");

      } else{
        addSmallDiag = as<bool>(costFuncObj["addSmallDiag"]);

      }

      if(!Rf_isNumeric(costFuncObj["epsilon"])){
        stop("epsilon must be a single non-negative numeric value!");

      } else{
        epsilon = as<double>(costFuncObj["epsilon"]);

        if(epsilon  <=  0.0){
          stop("epsilon must be a single positive numeric value!");

        }
      }

      Xnewptr = std::make_unique<Cost_SIGMA>(tsMat, addSmallDiag, epsilon, true);

    } else{
      stop("Either addSmallDiag or epsilon (or both) is missing!");

    }

  } else if(costFunc== "L2"){

    Xnewptr = std::make_unique<Cost_L2>(tsMat);

  } else {
    Rcpp::stop("Cost function not supported!");

  }

  CostBase& Xnew = *Xnewptr;

  const arma::uword nSamples = Xnew.nr;

  if(nSamples < 2*minSize){
    stop("Number of observations < 2*minSize!");
  }

  if(nSamples <= jump){
    stop("Number of observations <= jump!");
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

      if(end - lastBkp < minSize){ //Error message! However, by design, this should not happen.
        Rcpp::Rcout << "end - lastBkp < minSize!"<< std::endl;
        continue;
      }

      const double currentCost = Xnew.eval(lastBkp, end);

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
  return readPath(pathVec);
}
