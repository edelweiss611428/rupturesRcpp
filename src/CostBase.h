#ifndef COST_BASE_H
#define COST_BASE_H

class CostBase {
protected:

public:

  int nr;  // number of obs
  int nc;  // number of features

  virtual ~CostBase() = default;

  virtual double effEvalCpp(int start, int end,
                            bool addSmallDiag = true,
                            double epsilon = 1e-6) const = 0;
};

#endif // COST_BASE_H
