#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include "Cost.h"


// Function to return the maximum of two integers
inline int maxInt(const int a, const int b)
{
  return (a > b) ? a : b;
}

// Function to read the path from the path vector
std::vector<int> readPathL2(const arma::ivec& pathVec)
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
std::vector<int> peltL2(const arma::mat& tsMat, const double penalty, const int minSize, const int jump)
{

  Cost Xnew(tsMat);
  const arma::uword nSamples = Xnew.nr;

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
      const double currentCost = Xnew.effEvalCpp(lastBkp, end);
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
      if (const double currentCost = tmpCostVec[kLastBkp]; socVec[lastBkp] + currentCost <= minSoc)
      {
        admissibleBkps[nAdmissibleBkpsNew++] = lastBkp;
      }
    }

    // Add "end" to the set of admissible change-points
    admissibleBkps[nAdmissibleBkpsNew++] = static_cast<int>(std::floor(
      static_cast<double>(end - minSize + 1) / jump)
                                                              * jump);
    nAdmissibleBkps = nAdmissibleBkpsNew;
  }
  return readPathL2(pathVec);
}
