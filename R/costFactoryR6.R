#' `costFactory` class
#'
#' @description An `R6` class for fast segment-cost evaluation and parameter estimation, without a detection algorithm.
#'
#' @details
#' `costFactory` builds the `C++` cost module selected by `$costFunc` at `$fit()`, and keeps it, so every query reuses
#' its precomputations, as `PELT`, `binSeg` and `Window` do. `$eval()` and `$get_params()` call the module directly:
#' segments are `(a, b]` with 0-based `a`, and all checks are done in `C++`.
#'
#' `$new()` only stores the `costFunc` object; `$fit()` validates and attaches the data. `$costFunc` is an active
#' binding, so it can be inspected or replaced after construction -- if data has already been supplied, replacing it
#' automatically triggers `$fit()` again.
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
#'   \item{\code{$fit()}}{Constructs the `C++` cost module.}
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
#' cf = costFactory$new(costFunc$new("L2"))
#' cf$fit(tsMat)
#' cf$eval(0, 100)
#' cf$get_params(0, 100)
#'
#' # `costFunc` is an active binding: swapping it re-fits automatically.
#' cf$costFunc = costFunc$new("SIGMA")
#' cf$eval(0, 100)
#' @export
costFactory = R6Class(
  "costFactory",

  private = list(
    .costFunc = costFunc$new("L2"), #L2 cost function
    .module = NULL,
    .tsMat = NULL,
    .covariates = NULL,
    .fitted = FALSE,
    .n = NULL,
    .p = NULL
  ),

  active = list(

    #' @field costFunc `R6` object of class `costFunc`. Can be accessed or modified via `$costFunc`.
    #' Modifying `costFunc` will automatically trigger `$fit()` if a `tsMat` has already been fitted.
    costFunc = function(costFuncObj) {

      if (missing(costFuncObj)) {
        return(private$.costFunc)
      }

      if (!inherits(costFuncObj, "costFunc") || !is.R6(costFuncObj)) {
        stop("`costFunc` must be a `R6` object of class `costFunc` - can be created via costFunc$new()!")
      }

      private$.costFunc = costFuncObj

      # If time series data exists, refit the model
      if (!is.null(private$.tsMat) & private$.fitted) {
        message("`costFunc` has been updated. Re-fitting the model.")
        self$fit()
      }
    }

  ),

  public = list(

    #' @description Initialises a `costFactory` object. Does not build the `C++` cost module; call `$fit()` for that.
    #'
    #' @param costFunc A `R6` object of class `costFunc`. Should be created via `costFunc$new()` to avoid error.
    #' Default: `costFunc$new("L2")`.
    #'
    #' @return Invisibly returns `NULL`.
    initialize = function(costFunc) {

      if (!missing(costFunc)) {
        self$costFunc = costFunc
      }

      invisible(NULL)
    },

    #' @description Validates the supplied data and constructs the `C++` cost module selected by `$costFunc`.
    #'
    #' @param tsMat Numeric matrix. Time series of size \eqn{n \times p}. If `NULL`, the method will use the
    #' previously assigned `tsMat` (i.e., from a prior `$fit(tsMat)`). Default: `NULL`.
    #' @param covariates Numeric matrix with `n` rows, used by `"LinearL2"`, `"LinearSIGMA"` and `"LinearL1"`. If
    #' `NULL` and no prior `covariates` were set, the model is force-fitted with only an intercept. Default: `NULL`.
    #'
    #' @return Invisibly returns `NULL`.
    #'
    #' @details This method constructs the `C++` cost module and sets `private$.fitted` to `TRUE`, enabling the use
    #' of `$eval()` and `$get_params()`.
    fit = function(tsMat = NULL, covariates = NULL) {

      # Only assign if explicitly called with `tsMat` argument
      if (!is.null(tsMat)) {

        if (!is.numeric(tsMat) || !is.matrix(tsMat)) {
          stop("`tsMat` must be a numeric time series matrix!")
        }

        if (any(is.na(tsMat))) {
          stop("`tsMat` contains NAs!")
        }

        private$.tsMat = tsMat
        private$.n = nrow(tsMat)
        private$.p = ncol(tsMat)

      } else {

        if (is.null(private$.tsMat)) {
          stop("No `tsMat` found! Please provide a `tsMat`!")
        } #else not null because `private$.tsMat` was set by a prior `$fit()`

      }

      costFuncObj = private$.costFunc$pass()
      intercept = costFuncObj$intercept

      if (costFuncObj$costFunc %in% c("LinearL2", "LinearSIGMA", "LinearL1")) {

        if (!is.null(covariates)) {

          if (!is.numeric(covariates) || !is.matrix(covariates)) {
            stop("`covariates` must be a numeric time series matrix!")
          }

          if (any(is.na(covariates))) {
            stop("`covariates` contains NAs!")
          }

          if (nrow(covariates) != private$.n) {
            stop("Numbers of observations in `covariates` and `tsMat` do not match!")
          }

          private$.covariates = covariates

        } else if (!is.null(private$.covariates)) {

          if (nrow(private$.covariates) != private$.n) {
            stop("Numbers of observations in `covariates` and `tsMat` do not match!")
          }

          covariates = private$.covariates

        } else {

          warning("No `covariates` found! Force-fitting with only an intercept!")
          covariates = matrix(1, private$.n, 1)
          intercept = FALSE

        }
      }

      # warnOnce = FALSE: warn on every call, like $eval() of the detection classes
      private$.module = switch(costFuncObj$costFunc,
                               L1          = new(Cost_L1_cwMed, private$.tsMat, FALSE),
                               L2          = new(Cost_L2, private$.tsMat, FALSE),
                               SIGMA       = new(Cost_SIGMA, private$.tsMat, costFuncObj$addSmallDiag, costFuncObj$epsilon, FALSE),
                               VAR         = new(Cost_VAR, private$.tsMat, costFuncObj$pVAR, FALSE),
                               LinearL2    = new(Cost_LinearL2, private$.tsMat, covariates, intercept, FALSE),
                               LinearSIGMA = new(Cost_LinearSIGMA, private$.tsMat, covariates, intercept, costFuncObj$addSmallDiag, costFuncObj$epsilon, FALSE),
                               LinearL1    = new(Cost_LinearL1, private$.tsMat, covariates, intercept, costFuncObj$tol, costFuncObj$maxIter, FALSE),
                               Custom      = new(Cost_RFunc, private$.tsMat, costFuncObj$evalFun, costFuncObj$paramFun, FALSE))

      private$.fitted = TRUE

      invisible(NULL)
    },

    #' @description Evaluates the cost of the segment (a, b] with the module's `eval()`.
    #' @param a Integer. Start index (exclusive, 0-based).
    #' @param b Integer. End index (inclusive).
    #' @return The segment cost.
    eval = function(a, b) {

      if (!private$.fitted) {
        stop("`$fit()` must be run before `$eval()`!")
      }

      private$.module$eval(a, b)
    },

    #' @description Returns the parameter estimates of the segment (a, b] with the module's `get_params()`.
    #' @param a Integer. Start index (exclusive, 0-based).
    #' @param b Integer. End index (inclusive).
    #' @return A named list; see Details.
    get_params = function(a, b) {

      if (!private$.fitted) {
        stop("`$fit()` must be run before `$get_params()`!")
      }

      private$.module$get_params(a, b)
    }
  )
)
