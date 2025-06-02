#ifndef PELTL2_H
#define PELTL2_H

#include <vector>
#include <armadillo>

std::vector<int> peltL2(const arma::mat& signal, double penalty, int minSize, int jump);

#endif //PELTL2_H
