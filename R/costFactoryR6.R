#' `costFactory` class
#'
#' @description An `R6` class for fast segment-cost evaluation and parameter estimation, without a detection algorithm.
#'
#' @details
#' `costFactory` builds the `C++` cost module selected by a `costFunc` object once, at initialisation, and keeps it,
#' so every query reuses its precomputations, as `PELT`, `binSeg` and `Window` do. `$eval()` and `$get_params()` call
#' the module directly: segments are `(a, b]` with 0-based `a`, and all checks are done in `C++`.
#'
#' `$get_params()` returns:
#'
#' - `"L1"`: `median`.
#' - `"L2"`: `mean`.
#' - `"SIGMA"`: `mean` and `cov` (plus `epsilon` on the diagonal if `addSmallDiag = TRUE`).
#' - `"VAR"`, `"LinearL2"` and `"LinearL1"`: `coef`, intercept first.
#' - `"LinearSIGMA"`: `coef` and `cov`.
#' - `"Custom"`: `params`, the output of `paramFun`; an empty list if `paramFun` is `NULL`.
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
#' tsMat = cbind(c(rnorm(100, 0), rnorm(100, 5, 5)))
#' cf = costFactory$new(costFunc$new("L2"), tsMat)
#' cf$eval(0, 100)
#' cf$get_params(0, 100)
#' @export
costFactory = R6Class(
  "costFactory",

  private = list(
    .module = NULL
  ),

  public = list(

    #' @description Initialises a `costFactory` object and builds the `C++` cost module.
    #'
    #' @param costFunc A `R6` object of class `costFunc`. Default: `costFunc$new("L2")`.
    #' @param tsMat Numeric matrix. Time series of size \eqn{n \times p}.
    #' @param covariates Numeric matrix with `n` rows, used by `"LinearL2"`, `"LinearSIGMA"` and `"LinearL1"`. If `NULL`,
    #' the model is force-fitted with only an intercept. Default: `NULL`.
    #'
    #' @return Invisibly returns `NULL`.
    initialize = function(costFunc, tsMat, covariates = NULL) {

      spec = list(costFunc = "L2")
      if (!missing(costFunc)) {
        if (!inherits(costFunc, "costFunc") || !is.R6(costFunc)) {
          stop("`costFunc` must be a `R6` object of class `costFunc` - can be created via costFunc$new()!")
        }
        spec = costFunc$pass()
      }

      intercept = spec$intercept
      if (spec$costFunc %in% c("LinearL2", "LinearSIGMA", "LinearL1") && is.null(covariates)) {
        warning("No `covariates` found! Force-fitting with only an intercept!")
        covariates = matrix(1, nrow(tsMat), 1)
        intercept = FALSE
      }

      # warnOnce = FALSE: warn on every call, like $eval() of the detection classes
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

    #' @description Evaluates the cost of the segment (a, b] with the module's `eval()`.
    #' @param a Integer. Start index (exclusive, 0-based).
    #' @param b Integer. End index (inclusive).
    #' @return The segment cost.
    eval = function(a, b) {
      private$.module$eval(a, b)
    },

    #' @description Returns the parameter estimates of the segment (a, b] with the module's `get_params()`.
    #' @param a Integer. Start index (exclusive, 0-based).
    #' @param b Integer. End index (inclusive).
    #' @return A named list; see Details.
    get_params = function(a, b) {
      private$.module$get_params(a, b)
    }
  )
)
