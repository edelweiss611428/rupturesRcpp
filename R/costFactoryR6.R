#' `costFactory` class
#'
#' @description An `R6` class for fast segment-cost evaluation and parameter estimation, without a detection algorithm.
#'
#' @details
#' `costFactory` builds the `C++` cost module selected by a `costFunc` object once, at initialisation, and keeps it, so
#' every query reuses the module's precomputations (e.g. cumulative sums), as `PELT`, `binSeg` and `Window` do.
#' `$eval()` returns the cost of a segment and `$get_params()` the module's `get_params()` output, unchanged:
#'
#' - `"L1"`: `median`, the coordinate-wise median.
#' - `"L2"`: `mean`.
#' - `"SIGMA"`: `mean` and `cov`, the empirical covariance (divisor `n`), plus `epsilon` on the diagonal if
#'   `addSmallDiag = TRUE`.
#' - `"VAR"`: `coef`, a \eqn{(1 + rp) \times p} matrix: intercept, then lags `1, ..., r`.
#' - `"LinearL2"` and `"LinearL1"`: `coef`, a \eqn{J \times p} matrix, intercept first (if any).
#' - `"LinearSIGMA"`: `coef` and `cov`, the residual covariance (as for `"SIGMA"`).
#' - `"Custom"`: `params`, the output of `paramFun`; an empty list if `paramFun` is `NULL`.
#'
#' Vectors and matrices are returned as computed, without names. Coefficients (and `cov` for `"LinearSIGMA"`) are `NA`
#' when the segment has fewer observations than coefficients; singular systems fall back to an approximate solve, with
#' a warning.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `costFactory` object.}
#'   \item{\code{$eval()}}{Evaluates the cost of a segment.}
#'   \item{\code{$get_params()}}{Returns the parameter estimates of a segment.}
#'   \item{\code{$clone()}}{Clones the `R6` object.}
#' }
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @docType class
#' @include costFuncR6.R
#' @importFrom R6 R6Class is.R6
#'
#' @examples
#' set.seed(1)
#' tsMat = cbind(c(rnorm(100, 0), rnorm(100, 5, 5)), c(rnorm(100, 0), rnorm(100, 5, 5)))
#' cf = costFactory$new(costFunc$new("SIGMA"), tsMat)
#' cf$eval(0, 100)
#' cf$get_params(0, 100)
#' @export
costFactory = R6Class(
  "costFactory",

  private = list(

    .spec = list(costFunc = "L2"),   # costFunc$new("L2")$pass(); captured once so later edits to the costFunc object cannot desync the module
    .module = NULL,
    .n = NULL,

    .checkAB = function(a, b) {

      if (is.null(a) || is.null(b)) {
        stop("`a` and `b` must not be NULL")
      }
      if (!is.numeric(a) || length(a) != 1L || a < 0 || a > private$.n) {
        stop("`0 <= start < nSamples` must be true!")
      }
      if (!is.numeric(b) || length(b) != 1L || b < 0 || b > private$.n) {
        stop("`0 < end <= nSamples` must be true!")
      }
      a = as.integer(a)
      b = as.integer(b)
      if (a >= b) {
        stop("a must be smaller than b!")
      }
      c(a, b)
    }
  ),

  public = list(

    #' @description Initialises a `costFactory` object and builds the `C++` cost module.
    #'
    #' @param costFunc A `R6` object of class `costFunc`. Default: `costFunc$new("L2")`.
    #' @param tsMat Numeric matrix. Time series of size \eqn{n \times p}.
    #' @param covariates Numeric matrix with `n` rows. Used by `"LinearL2"`, `"LinearSIGMA"` and `"LinearL1"`; if `NULL`, the model is force-fitted with only an intercept. Default: `NULL`.
    #'
    #' @return Invisibly returns `NULL`.
    initialize = function(costFunc, tsMat, covariates = NULL) {

      if (!missing(costFunc)) {
        if (!inherits(costFunc, "costFunc") || !is.R6(costFunc)) {
          stop("`costFunc` must be a `R6` object of class `costFunc` - can be created via costFunc$new()!")
        }
        private$.spec = costFunc$pass()
      }

      if (missing(tsMat) || is.null(tsMat)) {
        stop("`tsMat` must be provided!")
      }
      if (!is.numeric(tsMat) || !is.matrix(tsMat)) {
        stop("`tsMat` must be a numeric time series matrix!")
      }
      if (anyNA(tsMat)) {
        stop("`tsMat` contains NAs!")
      }
      private$.n = nrow(tsMat)

      spec = private$.spec

      if (spec$costFunc %in% c("LinearL2", "LinearSIGMA", "LinearL1")) {

        if (is.null(covariates)) {
          warning("No `covariates` found! Force-fitting with only an intercept!")
          covariates = matrix(1, private$.n, 1)
          intercept = FALSE
        } else {
          if (!is.numeric(covariates) || !is.matrix(covariates)) {
            stop("`covariates` must be a numeric time series matrix!")
          }
          if (anyNA(covariates)) {
            stop("`covariates` contains NAs!")
          }
          if (nrow(covariates) != private$.n) {
            stop("Numbers of observations in `covariates` and `tsMat` do not match!")
          }
          intercept = spec$intercept
        }
      }

      private$.module = switch(spec$costFunc,
        L1          = new(Cost_L1_cwMed, tsMat, FALSE),
        L2          = new(Cost_L2, tsMat, FALSE),
        SIGMA       = new(Cost_SIGMA, tsMat, spec$addSmallDiag, spec$epsilon, FALSE),
        VAR         = new(Cost_VAR, tsMat, spec$pVAR, FALSE),
        LinearL2    = new(Cost_LinearL2, tsMat, covariates, intercept, FALSE),
        LinearSIGMA = new(Cost_LinearSIGMA, tsMat, covariates, intercept, spec$addSmallDiag, spec$epsilon, FALSE),
        LinearL1    = new(Cost_LinearL1, tsMat, covariates, intercept, spec$tol, spec$maxIter, FALSE),
        Custom      = new(Cost_RFunc, tsMat, spec$evalFun, spec$paramFun, FALSE))

      invisible(NULL)
    },

    #' @description Evaluates the cost of the segment (a, b].
    #' @param a Integer. Start index (exclusive).
    #' @param b Integer. End index (inclusive). Must satisfy `a < b`.
    #' @return The segment cost.
    eval = function(a, b) {
      ab = private$.checkAB(a, b)
      private$.module$eval(ab[1L], ab[2L])
    },

    #' @description Returns the parameter estimates of the segment (a, b]: the `get_params()` output of the `C++` cost module.
    #' @param a Integer. Start index (exclusive).
    #' @param b Integer. End index (inclusive). Must satisfy `a < b`.
    #' @return A named list; fields depend on the cost function, see Details.
    get_params = function(a, b) {
      ab = private$.checkAB(a, b)
      private$.module$get_params(ab[1L], ab[2L])
    }
  )
)
