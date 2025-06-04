#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#include "Cost.h"

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
// of O(k) (see ../deprecated); also supports minSize, jump, and linear penalty

inline Segment miniOptHeapCpp(const Cost& Xnew, const int& start, const int& end,
                               const int& minSize = 1,  const int& jump = 1,
                               double totalErr = -1) {

  int len = end - start;

  if(totalErr < 0){
    totalErr = Xnew.effEvalCpp(start, end);
  }

  if(len < 2*minSize){

    return Segment{start, end, true, start,
                   -9999,
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};
  } else if(len == 2*minSize){

    int cp = start + jump;
    double lErr = Xnew.effEvalCpp(start, cp);
    double rErr = Xnew.effEvalCpp(cp, end);
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

  for(int i = 0; i < allBkps.n_elem; i++){

    tempCp = allBkps(i);
    lErr = Xnew.effEvalCpp(start,tempCp);
    rErr = Xnew.effEvalCpp(tempCp,end);
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
List binSegCpp(const arma::mat& tsMat, const double& penalty = 0,
                const int& minSize = 1,  const int& jump = 1) {

  Cost Xnew(tsMat);
  int nr = Xnew.nr;

  if(nr < 2*minSize){
    stop("Number of observations < 2*minSize");
  }

  const int& maxNRegimes = std::floor(nr / minSize);
  NumericVector cost(maxNRegimes);
  double initCost = Xnew.effEvalCpp(0,nr);
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
    cost[idx] = cost[idx-1] - bestSeg.gain + penalty;
    Segment leftSeg = miniOptHeapCpp(Xnew, bestSeg.start, bestSeg.cp, minSize, jump, bestSeg.lErr);
    Segment rightSeg = miniOptHeapCpp(Xnew, bestSeg.cp, bestSeg.end, minSize, jump, bestSeg.rErr);
    heap.push(leftSeg);
    heap.push(rightSeg);

    if(heap.top().gain < 0){ //negative gain is equiv to no valid interval to split
      //check if nRegimes is correctly computed? should be nRegimes -=1;
      failed = true;
      break;
    }

    nRegimes++;

  } while(nRegimes <= maxNRegimes);
  //Need to check if the last cost has been computed? NO! => future feature

  if(not failed){
    changePoints[idx] = heap.top().cp;
    cost[idx] = cost[idx-1] - heap.top().gain + penalty;
  }


  return List::create(
    Named("nBkps") = nRegimes-1,
    Named("Bkps") = changePoints[Range(0,nRegimes-2)],
    Named("cost") = cost[Range(0,nRegimes-1)]
  );


}
