#ifndef COSTL2_H
#define COSTL2_H
#include <armadillo>
#include "utils.h"


class CostL2
{
public:
    explicit CostL2(const arma::mat& signal)
        : nSamples_(signal.n_rows), nDims_(signal.n_cols),
          csSignal_(getCumSum(signal)), csSquaredSignal_(getCumSum(arma::square(signal)))

    {
    }

    [[nodiscard]] double eval(arma::uword start, arma::uword end) const;

private:
    const arma::uword nSamples_;
    const arma::uword nDims_;
    const arma::mat csSignal_;
    const arma::mat csSquaredSignal_;
};


#endif //COSTL2_H
