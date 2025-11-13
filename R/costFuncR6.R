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

      if(any(!charVal %in% c("L1", "L2", "SIGMA", "VAR", "LinearL2"))){
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

      }
    }
  )
)
