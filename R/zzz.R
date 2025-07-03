#' @title Cost_L2 module
#' @description Load Rcpp module for Cost_L2
#' @importFrom Rcpp loadModule
#' @import methods
#' @details
#' **L2 cost function**:
#' \deqn{c_{L_2}(y_{(a+1)..b}) := \sum_{t = a+1}^{b} \| y_t - \bar{y}_{(a+1)..b} \|_2^2}
#' where \eqn{\bar{y}_{(a+1)..b}} is the empirical mean of the segment. If
#' `a+1 = b`, return 0.
#'
#' @export
#' @name Cost_L2
loadModule("Cost_L2_module", TRUE)

#' @title Cost_VAR module
#' @description Load Rcpp module for Cost_VAR
#' @importFrom Rcpp loadModule
#' @import methods
#' @details
#' **VAR cost function**:
#' \deqn{c_{\mathrm{VAR}}(y_{(a+1)..b}) := \sum_{t = a+r+1}^{b} \left\| y_t - \sum_{j=1}^r \hat A_j y_{t-j} \right\|_2^2}
#' where \eqn{\hat A_j} are the estimated VAR coefficients, commonly estimated via the OLS criterion. An approximate linear
#' solver will be used when exact `arma::solve()` fails. If no solution found, return 0. If `a-b < p*r+1` (i.e., not enough observations),
#' return 0.
#'
#' @export
#' @name Cost_VAR
loadModule("Cost_VAR_module", TRUE)

#' @title Cost_SIGMA module
#' @description Load Rcpp module for Cost_SIGMA
#' @importFrom Rcpp loadModule
#' @import methods
#' @details
#' **SIGMA cost function**:
#' \deqn{c_{\sum}(y_{(a+1)..b}) := (b - a)\log \det \hat{\Sigma}_{(a+1)..b}} where \eqn{\hat{\Sigma}_{(a+1)..b}} is
#' the empirical covariance matrix of the segment without Bessel correction. Here, if `addSmallDiag = TRUE`, a small
#' bias `epsilon` is added to the diagonal of estimated covariance matrices to improve numerical stability. If
#' \eqn{\hat{\Sigma}} is singular, return 0. If `a+1 = b`, return 0.
#'
#' @export
#' @name Cost_SIGMA
loadModule("Cost_SIGMA_module", TRUE)
