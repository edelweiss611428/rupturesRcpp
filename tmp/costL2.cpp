#include "costL2.h"


double CostL2::eval(const arma::uword start, const arma::uword end) const
{
    if (end > nSamples_ || start >= end)
    {
        throw std::invalid_argument("Invalid start or end indices.");
    }
    const double segmentLength = static_cast<double>(end - start);
    const double res = arma::sum(csSquaredSignal_.row(end) - csSquaredSignal_.row(start)) -
        arma::sum(arma::square(csSignal_.row(end) - csSignal_.row(start))) / segmentLength;

    return res;
}
