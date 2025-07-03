#include <RcppArmadillo.h>
#include <queue>
#include <limits>

#include "L2.h"
#include "SIGMA.h"
#include "VAR.h"
#include "baseClass.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// New structure for a segment
struct Segment {
  int start;
  int end;
  bool valid;
  int cp;
  double gain;
  double lErr;
  double rErr;
  double err;

  // Max-heap: higher gain has higher priority
  bool operator<(const Segment& other) const {
    return gain < other.gain;
  }
};

// Fast implementation based on heap, which involves an O(1) search instead
// of O(k) (see ../deprecated); also supports minSize, and jump.

inline Segment miniOptHeapCpp(const CostBase& Xnew, const int& start, const int& end,
                               const int& minSize = 1,  const int& jump = 1,
                               double totalErr = -1) {

  int len = end - start;

  if(totalErr < 0){
    totalErr = Xnew.eval(start, end);
  }

  if(len < 2*minSize){

    return Segment{start, end, true, start,
                   -9999, //If(segment length < 2*minSize)
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};
  } else if(len == 2*minSize){

    int cp = start + jump;
    double lErr = Xnew.eval(start, cp);
    double rErr = Xnew.eval(cp, end);
    double err = lErr + rErr;
    return Segment{start, end, true, start+jump,
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
    lErr = Xnew.eval(start,tempCp);
    rErr = Xnew.eval(tempCp,end);
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

// [[Rcpp::export]]
List binSegCpp(const arma::mat& tsMat, const int& minSize = 1,  const int& jump = 1,
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

        Xnewptr = std::make_unique<Cost_VAR>(tsMat, pVAR);

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

        if(epsilon  < 0){
          stop("epsilon must be a single non-negative numeric value!");

        }
      }

      Xnewptr = std::make_unique<Cost_SIGMA>(tsMat, addSmallDiag, epsilon);

    } else{
      stop("Either addSmallDiag or epsilon (or both) is missing!");

    }

  } else if(costFunc== "L2"){

    Xnewptr = std::make_unique<Cost_L2>(tsMat);

  } else {
    Rcpp::stop("Cost function not supported!");

  }

  CostBase& Xnew = *Xnewptr;

  int nr = Xnew.nr;

  if(nr < 2*minSize){
    stop("Number of observations < 2*minSize!");
  }

  if(nr <= jump){
    stop("Number of observations <= jump!");
  }

  const int& maxNRegimes = std::floor(nr / minSize);
  NumericVector cost(maxNRegimes);
  double initCost = Xnew.eval(0,nr);  //(0, nr]
  Segment seg0 = miniOptHeapCpp(Xnew, 0, nr, minSize, jump, initCost);
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

    Segment leftSeg = miniOptHeapCpp(Xnew, bestSeg.start, bestSeg.cp, minSize, jump, bestSeg.lErr);
    Segment rightSeg = miniOptHeapCpp(Xnew, bestSeg.cp, bestSeg.end, minSize, jump, bestSeg.rErr);
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


  return List::create(
    Named("bkps") = changePoints[Range(0,nRegimes-2)],
    Named("cost") = cost[Range(0,nRegimes-1)]
  );


}


// [[Rcpp::export]]
IntegerVector binSegPredCpp(const IntegerVector& bkps,
                            const NumericVector& cost,
                            const double& penalty = 0){

  if(penalty <0){
    stop("Penalty should be non-negative!");
  }
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

