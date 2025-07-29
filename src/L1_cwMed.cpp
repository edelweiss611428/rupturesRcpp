#include "L1_cwMed.h"

using namespace Rcpp;

// =========== Persistent Segment Tree and PSTManager methods ===========

// Node build
Node* PersistentSegmentTree::build(int l, int r) {
  return new Node(nullptr, nullptr, 0, 0.0);
}

Node* PersistentSegmentTree::insert(Node* prev, int l, int r, int val_idx, double val) {
  if (l + 1 == r) {
    return new Node(nullptr, nullptr, prev->count + 1, prev->sum + val);
  }
  int m = (l + r) / 2;
  if (val_idx < m) {
    return new Node(
        insert(prev->left ? prev->left : build(l, m), l, m, val_idx, val),
        prev->right,
        prev->count + 1,
        prev->sum + val
    );
  } else {
    return new Node(
        prev->left,
        insert(prev->right ? prev->right : build(m, r), m, r, val_idx, val),
        prev->count + 1,
        prev->sum + val
    );
  }
}

std::pair<int, double> PersistentSegmentTree::query(Node* u, Node* v, int l, int r, int k) {
  if (l + 1 == r) {
    return {l, 0.0};
  }
  int m = (l + r) / 2;
  int left_count = (v->left ? v->left->count : 0) - (u->left ? u->left->count : 0);
  if (k <= left_count) {
    return query(u->left ? u->left : build(l, m), v->left ? v->left : build(l, m), l, m, k);
  } else {
    return query(u->right ? u->right : build(m, r), v->right ? v->right : build(m, r), m, r, k - left_count);
  }
}

void PersistentSegmentTree::query_sum(Node* u, Node* v, int l, int r, int ql, int qr, double& count_out, double& sum_out) {
  if (r <= ql || qr <= l) return;
  if (ql <= l && r <= qr) {
    count_out += (v ? v->count : 0) - (u ? u->count : 0);
    sum_out += (v ? v->sum : 0.0) - (u ? u->sum : 0.0);
    return;
  }
  int m = (l + r) / 2;
  query_sum(u ? u->left : nullptr, v ? v->left : nullptr, l, m, ql, qr, count_out, sum_out);
  query_sum(u ? u->right : nullptr, v ? v->right : nullptr, m, r, ql, qr, count_out, sum_out);
}

void PersistentSegmentTree::build_tree(const arma::vec& x) {
  n = x.n_elem;
  sorted_vals = arma::conv_to<std::vector<double>>::from(x);
  std::sort(sorted_vals.begin(), sorted_vals.end());
  sorted_vals.erase(std::unique(sorted_vals.begin(), sorted_vals.end()), sorted_vals.end());

  Node* root = build(0, (int)sorted_vals.size());
  version_roots[-1] = root;

  for (int i = 0; i < n; ++i) {
    int idx = std::lower_bound(sorted_vals.begin(), sorted_vals.end(), x[i]) - sorted_vals.begin();
    version_roots[i] = insert(version_roots[i - 1], 0, (int)sorted_vals.size(), idx, x[i]);
  }
}

double PersistentSegmentTree::query_L1_cost(int l, int r) {
  int len = r - l + 1;
  int k = (len + 1) / 2; // median position

  Node* left = (l > 0) ? version_roots[l - 1] : version_roots[-1];
  Node* right = version_roots[r];

  auto [median_idx, _] = query(left, right, 0, (int)sorted_vals.size(), k);
  double median = sorted_vals[median_idx];

  double count_left = 0, sum_left = 0, count_right = 0, sum_right = 0;
  query_sum(left, right, 0, (int)sorted_vals.size(), 0, median_idx + 1, count_left, sum_left);
  query_sum(left, right, 0, (int)sorted_vals.size(), median_idx + 1, (int)sorted_vals.size(), count_right, sum_right);

  return median * count_left - sum_left + sum_right - median * count_right;
}

// PSTManager methods

void PSTManager::build(const arma::mat& X) {
  n_cols = X.n_cols;
  trees.clear();
  trees.resize(n_cols);
  for (int j = 0; j < n_cols; ++j) {
    trees[j].build_tree(X.col(j));
  }
}

double PSTManager::query_all(int l, int r) {
  double cost = 0.0;
  for (auto& tree : trees) {
    cost += tree.query_L1_cost(l, r);
  }
  return cost;
}

// =========== Cost_L1_cwMed methods ===========

Cost_L1_cwMed::Cost_L1_cwMed(const arma::mat& inputMat, bool warnOnce) {
  X = inputMat;
  nr = X.n_rows;
  warnOnce_ = warnOnce;
  keepWarning = !warnOnce;
  pst_manager.build(X);  // Build PST for all columns once
}

double Cost_L1_cwMed::eval(int start, int end) const {
  if (start >= end - 1) {
    return 0.0;
  }
  int l = start;
  int r = end - 1;
  return pst_manager.query_all(l, r);
}

void Cost_L1_cwMed::resetWarning(bool reset) {
  warnOnce_ = reset;
  keepWarning = !reset;
}

// =========== RCPP Module ===========

RCPP_EXPOSED_CLASS(Cost_L1_cwMed)
  RCPP_MODULE(Cost_L1_cwMed_module) {
    Rcpp::class_<Cost_L1_cwMed>("Cost_L1_cwMed")
    .constructor<arma::mat, bool>()
    .method("eval", &Cost_L1_cwMed::eval, "Evaluate L1 cost on interval (start, end]")
    .method("resetWarning", &Cost_L1_cwMed::resetWarning, "Set the status of warnOnce_");
  }
