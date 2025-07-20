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

  virtual void resetWarning(bool reset) {}// LCOV_EXCL_LINE
};

#endif // COST_BASE_H
