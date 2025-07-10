#' `costFunc` class
#'
#' @description An `R6` class specifying a cost function
#'
#' @details
#' Creates an instance of `costFunc` `R6` class, used in initialisation of change-point detection modules. Currently
#' supports the following cost functions:
#'
#' - **"L1"** and **"L2"** for (independent) piecewise Gaussian process with **constant variance**
#'
#'  \deqn{c_{L_1}(y_{(a+1)...b}) := \sum_{t = a+1}^{b} \| y_t - \tilde{y}_{(a+1)...b} \|_1}
#' where \eqn{\tilde{y}_{(a+1)...b}} is the coordinate-wise median of the segment. If \eqn{a \ge b - 1}, return 0.
#'
#' \deqn{c_{L_2}(y_{(a+1)...b}) := \sum_{t = a+1}^{b} \| y_t - \bar{y}_{(a+1)...b} \|_2^2}
#' where \eqn{\bar{y}_{(a+1)...b}} is the empirical mean of the segment. If \eqn{a \ge b - 1}, return 0. `"L2"`
#' is faster  but less robust to contamination than `"L1"`.
#'
#' - **"SIGMA"** for (independent) piecewise Gaussian process with **varying variance**
#'
#' \deqn{c_{\sum}(y_{(a+1)...b}) := (b - a)\log \det \hat{\Sigma}_{(a+1)...b}} where \eqn{\hat{\Sigma}_{(a+1)...b}} is
#' the empirical covariance matrix of the segment without Bessel's correction. Here, if `addSmallDiag = TRUE`, a small
#' bias `epsilon` is added to the diagonal of estimated covariance matrices to improve numerical stability. \cr
#' \cr
#' By default, `addSmallDiag = TRUE` and `epsilon = 1e-6`. In case `addSmallDiag = TRUE`, if the computed determinant of covariance matrix is either 0 (singular)
#' or smaller than `p*log(epsilon)` - the lower bound, return `(b - a)*p*log(epsilon)`, otherwise, output an error message.
#'
#' - **"VAR"** for piecewise Gaussian vector-regressive process with **constant noise variance**
#'
#' \deqn{c_{\mathrm{VAR}}(y_{(a+1)...b}) := \sum_{t = a+r+1}^{b} \left\| y_t - \sum_{j=1}^r \hat A_j y_{t-j} \right\|_2^2}
#' where \eqn{\hat A_j} are the estimated VAR coefficients, commonly estimated via the OLS criterion. If system is singular,
#' \eqn{a-b < p*r+1} (i.e., not enough observations), or \eqn{a \ge n-p} (where `n` is the time series length), return 0.
#'
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `costFunc` object.}
#'   \item{\code{$describe()}}{Describes the `costFunc` object.}
#'   \item{\code{$clone()}}{Clones the `costFunc` object.}
#' }
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom utils hasName
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

    #' @field costFunc Character. Cost function. Can be accessed or modified via `$costFunc`.
    costFunc = function(charVal) {

      if (missing(charVal)) {
        return(private$.costFunc)
      }

      if(!charVal %in% c("L1", "L2", "SIGMA", "VAR")){
        stop("Cost function not supported!")
      }

      if (!is.character(charVal) | length(charVal) != 1L) {
        stop("`costFunc` must be a single character value!")
      }

      private$.costFunc = charVal
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

      }
    }
  )
)
