#ifndef COST_L1_CW_MEDIAN_H
#define COST_L1_CW_MEDIAN_H

#include "baseClass.h"
#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <unordered_map>

// Persistent Segment Tree Node
struct Node {
  Node* left = nullptr;
  Node* right = nullptr;
  int count = 0;
  double sum = 0.0;

  Node() = default;
  Node(Node* l, Node* r, int c, double s) : left(l), right(r), count(c), sum(s) {}
};

// Persistent Segment Tree class (univariate)
class PersistentSegmentTree {
private:
  std::vector<double> sorted_vals;
  std::unordered_map<int, Node*> version_roots;
  int n;

  Node* build(int l, int r);
  Node* insert(Node* prev, int l, int r, int val_idx, double val);
  std::pair<int, double> query(Node* u, Node* v, int l, int r, int k);
  void query_sum(Node* u, Node* v, int l, int r, int ql, int qr, double& count_out, double& sum_out);

public:
  void build_tree(const arma::vec& x);
  double query_L1_cost(int l, int r);
};

// PSTManager manages multiple PSTrees (one per column)
class PSTManager {
private:
  std::vector<PersistentSegmentTree> trees;
  int n_cols = 0;

public:
  void build(const arma::mat& X);
  double query_all(int l, int r);
};

// ========================================================
// Cost_L1_cwMed using PSTManager internally
// ========================================================

class Cost_L1_cwMed : public CostBase {
private:
  arma::mat X;
  bool warnOnce_;
  bool keepWarning;

  mutable PSTManager pst_manager;  // PSTManager is a private member

public:
  explicit Cost_L1_cwMed(const arma::mat& inputMat, bool warnOnce = true);

  double eval(int start, int end) const override;

  void resetWarning(bool reset) override;

  virtual ~Cost_L1_cwMed() = default;
};

#endif // COST_L1_CW_MEDIAN_H
