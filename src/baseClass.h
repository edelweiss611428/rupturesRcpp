#ifndef COST_BASE_H
#define COST_BASE_H

// ========================================================
//                 CostBase abstract class
// ========================================================

class CostBase {
protected:

public:

  int nr;  // Number of observations
  int nc;  // Number of features

  virtual ~CostBase() = default;

  virtual double eval(int start, int end) const = 0;
  // LCOV_EXCL_START
  virtual void resetWarning(bool reset) {}
  // LCOV_EXCL_STOP
};

#endif // COST_BASE_H
