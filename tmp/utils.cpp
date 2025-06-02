#include <armadillo>
#include "utils.h"

arma::mat loadSignal(const std::string& filePath) {
    arma::mat signal;
    signal.load(filePath, arma::csv_ascii);
    return signal;
}

/**
 * @brief Creates a matrix of autoregressors from the input signal.
 *
 * This function takes an input matrix `signal` and generates a matrix of autoregressors
 * with the specified number of lags. If `intercept` is true, an intercept term is added.
 *
 * @param signal The input matrix of shape (nSamples, nDims).
 * @param nLags The number of lags to include in the autoregressors.
 * @param intercept A boolean indicating whether to include an intercept term.
 * @return A matrix of shape (nSamples, nDims * nRegs) containing the autoregressors.
 */
arma::mat makeAutoRegressors(const arma::mat &signal, const int nLags, const bool intercept) {
    const arma::uword nSamples = signal.n_rows;
    const arma::uword nDims = signal.n_cols;
    const arma::uword nRegs = intercept ? nLags + 1 : nLags;

    arma::mat out(nSamples, nDims * nRegs, arma::fill::zeros);

    for (arma::uword t = nLags; t < nSamples; ++t) {
        for (int lag = 0; lag < nLags; ++lag) {
            out.submat(t, lag * nDims, t, (lag + 1) * nDims - 1) = signal.row(t - lag - 1);
        }
    }

    if (intercept) {
        out.submat(0, nDims * nLags, nSamples - 1, nDims * nRegs - 1).ones();
    }

    return out;
}

/**
 * @brief Computes the cumulative sum of the rows of a matrix.
 *
 * This function takes a matrix X of shape (nSamples, nDims) and returns a matrix
 * of shape (nSamples + 1, nDims) where each row is the cumulative sum of the rows
 * of X up to that point. The first row of the output matrix is initialized to zero.
 *
 * @param X The input matrix of shape (nSamples, nDims).
 * @return A matrix of shape (nSamples + 1, nDims) containing the cumulative sums.
 */
arma::mat getCumSum(const arma::mat& X) {
    // X, shape (nSamples, nDims)
    // out, shape (nSamples + 1, nDims)
    const arma::uword nSamples = X.n_rows;
    const arma::uword nDims = X.n_cols;
    arma::mat out(nSamples + 1, nDims, arma::fill::zeros);
    arma::rowvec accum = arma::zeros<arma::rowvec>(nDims);
    for (arma::uword t = 1; t <= nSamples; ++t) {
        accum += X.row(t - 1);
        out.row(t) = accum;
    }
    return out;
}