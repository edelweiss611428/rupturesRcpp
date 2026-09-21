#' `costFactory` class
#'
#' @description An `R6` class for fast segment-cost evaluation and parameter estimation, without a detection algorithm.
#'
#' @details
#' `costFactory` builds the `C++` cost module selected by a `costFunc` object and exposes `$eval()` (segment cost),
#' `$estimate()` (fitted parameters of one segment) and `$summary()` (a `segSummary` for a given segmentation, see
#' [segSummary]). Costs are identical to those used by `PELT`, `binSeg` and `Window`. Parameters are estimated in `R`
#' on the same model:
#'
#' - **L1**: coordinate-wise median \eqn{\tilde y} and scale \eqn{\frac{1}{n}\sum_t |y_t - \tilde y|}.
#' - **L2**: mean \eqn{\bar y} and variance \eqn{\frac{1}{n}\sum_t (y_t - \bar y)^2}.
#' - **SIGMA**: mean and covariance \eqn{\frac{1}{n}\sum_t (y_t - \bar y)(y_t - \bar y)^\top} (without `epsilon`).
#' - **VAR(r)**: OLS coefficients \eqn{\hat B} of \eqn{y_t} on \eqn{(1, y_{t-1}, \dots, y_{t-r})} and residual variance.
#' - **LinearL2**: OLS coefficients of \eqn{y_t} on \eqn{(1, x_t)} and residual variance. With no `covariates`,
#'   an intercept-only model is fitted, as in `$fit()` of the detection classes.
#'
#' Coefficients are `NA` when a segment is too short or the system is rank deficient.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `costFactory` object.}
#'   \item{\code{$eval()}}{Evaluates the cost of a segment.}
#'   \item{\code{$estimate()}}{Estimates the parameters of a segment.}
#'   \item{\code{$summary()}}{Summarises a segmentation.}
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
#' cf$estimate(0, 100)
#' cf$summary(c(100, 200))
#' @export
costFactory = R6Class(
  "costFactory",

  private = list(

    .costFunc = costFunc$new("L2"),   # default set here, not as an argument default (see initialize)
    .module = NULL,
    .tsMat = NULL,
    .covariates = NULL,
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
    #' @param covariates Numeric matrix with `n` rows. Required for `"LinearL2"`; if `NULL`, the model is
    #' force-fitted with only an intercept. Default: `NULL`.
    #'
    #' @return Invisibly returns `NULL`.
    initialize = function(costFunc, tsMat, covariates = NULL) {

      if (!missing(costFunc)) {
        if (!inherits(costFunc, "costFunc") || !is.R6(costFunc)) {
          stop("`costFunc` must be a `R6` object of class `costFunc` - can be created via costFunc$new()!")
        }
        private$.costFunc = costFunc
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
      private$.tsMat = tsMat
      private$.n = nrow(tsMat)

      spec = private$.costFunc$pass()

      if (spec$costFunc == "LinearL2") {

        if (is.null(covariates)) {
          warning("No `covariates` found! Force-fitting with only an intercept!")
          private$.module = new(Cost_LinearL2, tsMat, matrix(1, private$.n, 1), FALSE, FALSE)
          return(invisible(NULL))
        }
        if (!is.numeric(covariates) || !is.matrix(covariates)) {
          stop("`covariates` must be a numeric time series matrix!")
        }
        if (anyNA(covariates)) {
          stop("`covariates` contains NAs!")
        }
        if (nrow(covariates) != private$.n) {
          stop("Numbers of observations in `covariates` and `tsMat` do not match!")
        }
        private$.covariates = covariates
        private$.module = new(Cost_LinearL2, tsMat, covariates, spec$intercept, FALSE)

      } else if (spec$costFunc == "L2") {
        private$.module = new(Cost_L2, tsMat, FALSE)

      } else if (spec$costFunc == "L1") {
        private$.module = new(Cost_L1_cwMed, tsMat, FALSE)

      } else if (spec$costFunc == "SIGMA") {
        private$.module = new(Cost_SIGMA, tsMat, spec$addSmallDiag, spec$epsilon, FALSE)

      } else if (spec$costFunc == "VAR") {
        private$.module = new(Cost_VAR, tsMat, spec$pVAR, FALSE)
      }

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

    #' @description Estimates the parameters of the segment (a, b].
    #' @param a Integer. Start index (exclusive).
    #' @param b Integer. End index (inclusive). Must satisfy `a < b`.
    #' @return A named list; entries depend on the cost function, see [segSummary].
    estimate = function(a, b) {
      ab = private$.checkAB(a, b)
      .segEstimate(private$.costFunc$pass(), private$.tsMat, private$.covariates, ab[1L], ab[2L])
    },

    #' @description Summarises a segmentation: cost and parameters of every segment.
    #' @param endPts Integer vector. Sorted end points; the last element must be `n`.
    #' @return An object of class `segSummary`, see [segSummary].
    summary = function(endPts) {
      .segSummary(self$eval, private$.costFunc$pass(), private$.tsMat, private$.covariates, endPts)
    }
  )
)
