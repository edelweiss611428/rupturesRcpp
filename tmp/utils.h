#ifndef UTILS_H
#define UTILS_H
#include <armadillo>

arma::mat loadSignal(const std::string& filePath);
arma::mat makeAutoRegressors(const arma::mat &signal, const int nLags, const  bool intercept);
arma::mat getCumSum(const arma::mat& X);

#endif // UTILS_H