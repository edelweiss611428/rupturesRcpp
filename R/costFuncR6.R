#' `costFunc` class
#'
#' @description An `R6` class specifying a cost function
#'
#' @details
#' Creates an instance of `costFunc` `R6` class, used in initialisation of change-point detection modules. Currently
#' supports the following cost functions:
#'
#'
#' - **L1 cost function**:
#' \deqn{c_{L_1}(y_{(a+1)...b}) := \sum_{t = a+1}^{b} \| y_t - \tilde{y}_{(a+1)...b} \|_1}
#' where \eqn{\tilde{y}_{(a+1)...b}} is the coordinate-wise median of the segment. If \eqn{a \ge b - 1}, return 0.
#'
#' - **L2 cost function**:
#' \deqn{c_{L_2}(y_{(a+1)...b}) := \sum_{t = a+1}^{b} \| y_t - \bar{y}_{(a+1)...b} \|_2^2}
#' where \eqn{\bar{y}_{(a+1)...b}} is the empirical mean of the segment. If \eqn{a \ge b - 1}, return 0.
#'
#' - **SIGMA cost function**:
#' \deqn{c_{\sum}(y_{(a+1)...b}) := (b - a)\log \det \hat{\Sigma}_{(a+1)...b}} where \eqn{\hat{\Sigma}_{(a+1)...b}} is
#' the empirical covariance matrix of the segment without Bessel's correction. Here, if `addSmallDiag = TRUE`, a small
#' bias `epsilon` is added to the diagonal of estimated covariance matrices to improve numerical stability. \cr
#' \cr
#' By default, `addSmallDiag = TRUE` and `epsilon = 1e-6`. In case `addSmallDiag = TRUE`, if the computed determinant of covariance matrix is either 0 (singular)
#' or smaller than `p*log(epsilon)` - the lower bound, return `(b - a)*p*log(epsilon)`, otherwise, output an error message.
#'
#' - **VAR(r) cost function**:
#' \deqn{c_{\mathrm{VAR}}(y_{(a+1)...b}) := \sum_{t = a+r+1}^{b} \left\| y_t - \sum_{j=1}^r \hat A_j y_{t-j} \right\|_2^2}
#' where \eqn{\hat A_j} are the estimated VAR coefficients, commonly estimated via the OLS criterion. If system is singular,
#' \eqn{a-b < p*r+1} (i.e., not enough observations), or \eqn{a \ge n-p} (where `n` is the time series length), return 0.
#'
#' - **"LinearL2"** for piecewise linear regression process with **constant noise variance**
#' \deqn{c_{\text{LinearL2}}(y_{(a+1):b}) := \sum_{t=a+1}^b \| y_t - X_t \hat{\beta} \|_2^2} where \eqn{\hat{\beta}} are OLS estimates on segment \eqn{(a+1):b}. If segment is shorter than the minimum number of
#' points needed for OLS, return 0.
#'
#' - **"LinearSIGMA"** for piecewise linear regression process with **varying noise covariance**
#' \deqn{c_{\text{LinearSIGMA}}(y_{(a+1):b}) := (b-a)\log \det \hat\Sigma_{(a+1):b}} where \eqn{\hat\Sigma_{(a+1):b}}
#' is the empirical covariance matrix of OLS residuals \eqn{y - X\hat{\beta}} on segment \eqn{(a+1):b}, estimated
#' the same way as in the SIGMA cost function (including the `addSmallDiag`/`epsilon` stabilisation and lower-bound fallback).
#'
#' - **"LinearL1"** for piecewise linear regression process under **L1 (least absolute deviations) loss**
#' \deqn{c_{\text{LinearL1}}(y_{(a+1):b}) := \sum_{t=a+1}^b \| y_t - X_t \hat{\beta} \|_1} where \eqn{\hat{\beta}} is
#' fit column-by-column via Iteratively Reweighted Least Squares (IRLS), iterated until the change in the fit's cost
#' is within `tol` or `maxIter` iterations are reached. Unlike the other regression cost functions, this has no
#' \eqn{O(1)}-per-segment closed form, since IRLS must be re-run on each queried segment.
#'
#' - **"Custom"** for a user-defined cost function supplied from R
#' \deqn{c_{\text{Custom}}(y_{(a+1):b}) := \text{evalFun}(y_{(a+1):b}, a, b)} where `evalFun` is a
#' user-supplied function called on the raw segment matrix, plus the segment's own `(a,b]` bounds --
#' the latter let `evalFun` align the segment against any externally-captured, position-indexed data
#' (e.g. a weight vector or exogenous series `evalFun` closes over) without the package needing to know
#' that data exists. Because each call crosses back into R, this is substantially slower per call than
#' the built-in costs above; see `$evalFun` and `$paramFun`.
#'
#' If active binding `$costFunc` is modified (via assignment operator), the default parameters will be used.
#'
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `costFunc` object.}
#'   \item{\code{$pass()}}{Describes the `costFunc` object.}
#'   \item{\code{$clone()}}{Clones the `costFunc` object.}
#' }
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom utils hasName
#'
#' @examples
#'
#' ## L2 costFunc (default)
#' costFuncObj = costFunc$new()
#' costFuncObj$pass()
#' ## SIGMA costFunc
#' costFuncObj = costFunc$new(costFunc = "SIGMA")
#' costFuncObj$pass()
#' # Modify active bindings
#' costFuncObj$epsilon = 10^-5
#' costFuncObj$pass()
#' costFuncObj$costFunc = "VAR"
#' costFuncObj$pass()
#' @export
#'
#'
costFunc <- R6::R6Class(
  "costFunc",

  private = list(
    .costFunc = "L2",
    .params = list()
  ),

  active = list(

    #' @field costFunc Character. Cost function. Can be accessed or modified via `$costFunc`. If `costFunc` is modified
    #' and required parameters are missing, the default parameters are used.
    costFunc = function(charVal) {

      if (missing(charVal)) {
        return(private$.costFunc)
      }

      if(any(!charVal %in% c("L1", "L2", "SIGMA", "VAR", "LinearL2", "LinearSIGMA", "LinearL1", "Custom"))){
        stop("Cost function not supported!")
      }

      if (!is.character(charVal) | length(charVal) != 1L) {
        stop("`costFunc` must be a single character value!")
      }

      private$.costFunc = charVal

      # Set default values for required parameters if missing
      if (charVal == "VAR") {
        if (is.null(private$.params[["pVAR"]])) {
          private$.params[["pVAR"]] = 1L
        }
      }

      if (charVal == "SIGMA") {
        if (is.null(private$.params[["addSmallDiag"]])) {
          private$.params[["addSmallDiag"]] = TRUE
        }
        if (is.null(private$.params[["epsilon"]])) {
          private$.params[["epsilon"]] = 1e-6
        }
      }

      if (charVal == "LinearL2") {
        if (is.null(private$.params[["intercept"]])) {
          private$.params[["intercept"]] = TRUE
        }
      }

      if (charVal == "LinearSIGMA") {
        if (is.null(private$.params[["intercept"]])) {
          private$.params[["intercept"]] = TRUE
        }
        if (is.null(private$.params[["addSmallDiag"]])) {
          private$.params[["addSmallDiag"]] = TRUE
        }
        if (is.null(private$.params[["epsilon"]])) {
          private$.params[["epsilon"]] = 1e-6
        }
      }

      if (charVal == "LinearL1") {
        if (is.null(private$.params[["intercept"]])) {
          private$.params[["intercept"]] = TRUE
        }
        if (is.null(private$.params[["tol"]])) {
          private$.params[["tol"]] = 1e-6
        }
        if (is.null(private$.params[["maxIter"]])) {
          private$.params[["maxIter"]] = 50L
        }
      }
    },

    #' @field pVAR Integer. Vector autoregressive order. Can be accessed or modified via `$pVAR`.
    pVAR = function(intVal) {

      if (missing(intVal)) {
        return(private$.params[["pVAR"]])
      }

      if (!is.numeric(intVal) | length(intVal) != 1L | any(as.integer(intVal) < 1L)) {
        stop("`pVAR` must be a single positive integer!")

      }
      private$.params[["pVAR"]] = as.integer(intVal)

    },

    #' @field addSmallDiag Logical. Whether to add a bias value to the diagonal of estimated covariance matrices to stabilise matrix operations. Can be accessed or modified via `$addSmallDiag`.
    addSmallDiag = function(boolVal) {

      if (missing(boolVal)) {
        return(private$.params[["addSmallDiag"]])
      }

      if (!is.logical(boolVal) | length(boolVal) != 1L) {
        stop("`addSmallDiag` must be a single boolean value!")

      }
      private$.params[["addSmallDiag"]] = boolVal

    },

    #' @field epsilon Double. A bias value added to the diagonal of estimated covariance matrices to stabilise matrix operations. Can be accessed or modified via `$epsilon`.
    epsilon = function(doubleVal) {

      if (missing(doubleVal)) {
        return(private$.params[["epsilon"]])
      }

      if (!is.numeric(doubleVal) | length(doubleVal) != 1L | any(doubleVal <= 0)) {
        stop("`epsilon` must be single positive value!")

      }
      private$.params[["epsilon"]] = doubleVal

    },

    #' @field intercept Logical. Whether to include the intercept in regression problems. Can be accessed or modified via `$intercept`.
    intercept = function(boolVal) {

      if (missing(boolVal)) {
        return(private$.params[["intercept"]])
      }

      if (!is.logical(boolVal) | length(boolVal) != 1L) {
        stop("`intercept` must be a single boolean value!")

      }
      private$.params[["intercept"]] = boolVal

    },

    #' @field tol Double. IRLS convergence tolerance: iteration stops once the change in the fit's cost falls below
    #' `tol`. Can be accessed or modified via `$tol`.
    tol = function(doubleVal) {

      if (missing(doubleVal)) {
        return(private$.params[["tol"]])
      }

      if (!is.numeric(doubleVal) | length(doubleVal) != 1L | any(doubleVal <= 0)) {
        stop("`tol` must be a single positive value!")

      }
      private$.params[["tol"]] = doubleVal

    },

    #' @field maxIter Integer. Maximum number of IRLS iterations. Can be accessed or modified via `$maxIter`.
    maxIter = function(intVal) {

      if (missing(intVal)) {
        return(private$.params[["maxIter"]])
      }

      valid = is.numeric(intVal) && length(intVal) == 1L && !anyNA(intVal) &&
        intVal >= 1 && intVal == round(intVal) && intVal <= .Machine$integer.max

      if (!isTRUE(valid)) {
        stop("`maxIter` must be a single positive integer!")

      }
      private$.params[["maxIter"]] = as.integer(intVal)

    },

    #' @field evalFun Function. Required for `costFunc = "Custom"`. A user-defined cost function, called
    #' as `evalFun(segment, a, b)`, where `segment` is the numeric matrix of rows `(a+1):b` for the
    #' queried segment `(a,b]` (0-indexed, same convention as `$eval(a, b)`). `a` and `b` let `evalFun`
    #' align `segment` against externally-captured, position-indexed data it closes over (e.g.
    #' `externalSeries[(a+1):b]`), which the package itself never needs to see. Must return a single
    #' numeric value. Can be accessed or modified via `$evalFun`.
    evalFun = function(funVal) {

      if (missing(funVal)) {
        return(private$.params[["evalFun"]])
      }

      if (!is.function(funVal)) {
        stop("`evalFun` must be a function!")
      }

      private$.params[["evalFun"]] = funVal

    },

    #' @field paramFun Function or `NULL`. Optional for `costFunc = "Custom"`. A user-defined function
    #' called as `paramFun(segment, a, b)` (same convention as `evalFun`), used by `$get_params()` to
    #' report segment-level estimates. If `NULL` (default), `$get_params()` returns an empty list for
    #' `"Custom"`. Can be accessed or modified via `$paramFun`.
    paramFun = function(funVal) {

      if (missing(funVal)) {
        return(private$.params[["paramFun"]])
      }

      if (!is.null(funVal) && !is.function(funVal)) {
        stop("`paramFun` must be a function or NULL!")
      }

      private$.params[["paramFun"]] = funVal

    }

  ),

  public = list(

    #' @description Initialises a `costFunc` object.
    #'
    #' @param costFunc Character. Cost function. Supported values include \code{"L2"}, \code{"VAR"},
    #' and \code{"SIGMA"}. Default: `L2`.
    #' @param ... Optional named parameters required by specific cost functions. \cr
    #' If any required parameters are missing or null, default values will be used.
    #'
    #' For \code{"L1"} and \code{"L2"}, there is no extra parameter.
    #'
    #' For \code{"SIGMA"}, supported parameters are:
    #' \describe{
    #'   \item{`addSmallDiag`}{Logical. If \code{TRUE}, add a small value to the diagonal of estimated covariance matrices
    #'   to stabilise matrix operations. Default: `TRUE`.}
    #'   \item{`epsilon`}{Double. If `addSmallDiag = TRUE`, a small positive value added to the diagonal of estimated covariance matrices to stabilise
    #'   matrix operations. Default: `1e-6`.}
    #' }
    #'
    #' For \code{"VAR"}, \code{pVAR} is required:
    #' \describe{
    #'   \item{`pVAR`}{Integer. Vector autoregressive order. Must be a positive integer. Default: `1L`.}
    #' }
    #'
    #' For \code{"LinearL2"}, \code{intercept} is required:
    #' \describe{
    #'   \item{`intercept`}{Logical. Whether to include the intercept in regression problems. Default: `TRUE`.}
    #' }
    #'
    #' For \code{"LinearSIGMA"}, supported parameters are:
    #' \describe{
    #'   \item{`intercept`}{Logical. Whether to include the intercept in regression problems. Default: `TRUE`.}
    #'   \item{`addSmallDiag`}{Logical. If \code{TRUE}, add a small value to the diagonal of estimated residual covariance matrices
    #'   to stabilise matrix operations. Default: `TRUE`.}
    #'   \item{`epsilon`}{Double. If `addSmallDiag = TRUE`, a small positive value added to the diagonal of estimated residual covariance matrices to stabilise
    #'   matrix operations. Default: `1e-6`.}
    #' }
    #'
    #' For \code{"LinearL1"}, supported parameters are:
    #' \describe{
    #'   \item{`intercept`}{Logical. Whether to include the intercept in regression problems. Default: `TRUE`.}
    #'   \item{`tol`}{Double. IRLS convergence tolerance: iteration stops once the change in the fit's cost falls
    #'   below `tol`. Default: `1e-6`.}
    #'   \item{`maxIter`}{Integer. Maximum number of IRLS iterations. Default: `50L`.}
    #' }
    #'
    #' For \code{"Custom"}, supported parameters are:
    #' \describe{
    #'   \item{`evalFun`}{Function. Required. See `$evalFun` for details.}
    #'   \item{`paramFun`}{Function or `NULL`. Optional. See `$paramFun` for details.}
    #' }

    initialize = function(costFunc, ...) {

      if(!missing(costFunc)){
        self$costFunc = costFunc
      }

      args = list(...)

      if (private$.costFunc == "VAR") {

        if (hasName(args, "pVAR") & !is.null(args$pVAR)) {
          self$pVAR = args$pVAR

        } else {
          self$pVAR = 1L

        }
      }

      if (private$.costFunc == "SIGMA") {

        self$addSmallDiag = if (hasName(args, "addSmallDiag") & !is.null(args$addSmallDiag)) {
          args$addSmallDiag

        } else {
          TRUE

        }

        self$epsilon = if (hasName(args, "epsilon") & !is.null(args$epsilon)) {
          args$epsilon

        } else {
          1e-6

        }
      }

      if (private$.costFunc == "LinearL2") {

        self$intercept= if (hasName(args, "intercept") & !is.null(args$intercept)) {
          args$intercept

        } else {
          TRUE

        }
      }

      if (private$.costFunc == "LinearSIGMA") {

        self$intercept = if (hasName(args, "intercept") & !is.null(args$intercept)) {
          args$intercept

        } else {
          TRUE

        }

        self$addSmallDiag = if (hasName(args, "addSmallDiag") & !is.null(args$addSmallDiag)) {
          args$addSmallDiag

        } else {
          TRUE

        }

        self$epsilon = if (hasName(args, "epsilon") & !is.null(args$epsilon)) {
          args$epsilon

        } else {
          1e-6

        }
      }

      if (private$.costFunc == "LinearL1") {

        self$intercept = if (hasName(args, "intercept") & !is.null(args$intercept)) {
          args$intercept

        } else {
          TRUE

        }

        self$tol = if (hasName(args, "tol") & !is.null(args$tol)) {
          args$tol

        } else {
          1e-6

        }

        self$maxIter = if (hasName(args, "maxIter") & !is.null(args$maxIter)) {
          args$maxIter

        } else {
          50L

        }
      }

      if (private$.costFunc == "Custom") {

        if (hasName(args, "evalFun") & !is.null(args$evalFun)) {
          self$evalFun = args$evalFun
        }

        if (hasName(args, "paramFun") & !is.null(args$paramFun)) {
          self$paramFun = args$paramFun
        }
      }

    },

    #' @description Returns a list of configuration parameters to initialise `detection` modules.
    #'
    pass = function() {

      if(private$.costFunc == "L2"){
        return(list(costFunc = "L2"))

      } else if (private$.costFunc == "L1"){
        return(list(costFunc = "L1"))

      } else if(private$.costFunc == "VAR"){
        return(list(costFunc = "VAR",
                    pVAR = private$.params[["pVAR"]]))

      } else if(private$.costFunc == "SIGMA"){
        return(list(costFunc = "SIGMA",
                    addSmallDiag = private$.params[["addSmallDiag"]],
                    epsilon = private$.params[["epsilon"]]))

      } else if(private$.costFunc == "LinearL2"){
        return(list(costFunc = "LinearL2",
                    intercept = private$.params[["intercept"]]))

      } else if(private$.costFunc == "LinearSIGMA"){
        return(list(costFunc = "LinearSIGMA",
                    intercept = private$.params[["intercept"]],
                    addSmallDiag = private$.params[["addSmallDiag"]],
                    epsilon = private$.params[["epsilon"]]))

      } else if(private$.costFunc == "LinearL1"){
        return(list(costFunc = "LinearL1",
                    intercept = private$.params[["intercept"]],
                    tol = private$.params[["tol"]],
                    maxIter = private$.params[["maxIter"]]))

      } else if(private$.costFunc == "Custom"){

        if(is.null(private$.params[["evalFun"]])){
          stop("`costFunc = \"Custom\"` requires `evalFun` to be set first, e.g. via `$evalFun <- function(segment, a, b) ...`!")
        }

        return(list(costFunc = "Custom",
                    evalFun = private$.params[["evalFun"]],
                    paramFun = private$.params[["paramFun"]]))

      }
    }
  )
)
