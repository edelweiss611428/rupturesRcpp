#' `costFunc` class
#'
#' @description An `R6` class specifying a cost function
#'
#' @details
#' Creates an instance of `costFunc` `R6` class, used in initialisation of change-point detection modules.
#'
#' Currently supports the following cost functions:
#'
#' - `"L1"` and `"L2"` for (independent) piecewise Gaussian process with **constant variance**
#' - `"SIGMA"`: for (independent) piecewise Gaussian process with **varying variance**
#' - `"VAR"`: for piecewise Gaussian vector-regressive process with **constant noise variance**
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

    #' @field coorWise Logical. If `TRUE` coordinate-wise median is computed. Otherwise, geometric median is computed via
    #' the Weiszfeld algorithm. Can be accessed or modified via `$coorWise`.
    #'
    coorWise = function(boolVal) {

      if (missing(boolVal)) {
        return(private$.params[["coorWise"]])
      }

      if (!is.logical(boolVal) | length(boolVal) != 1L) {
        stop("`coorWise` must be a single boolean value!")

      }
      private$.params[["coorWise"]] = boolVal

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
    #' #' For \code{"L1"}, supported parameters are:
    #' \describe{
    #'   \item{`coorWise`}{Logical. If \code{TRUE}, coordinate-wise median is computed. If \code{FALSE}, geometric median is computed via the
    #'   Weiszfeld algorithm. Default: `TRUE`.}
    #' }
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

      if (private$.costFunc == "L1") {

        if (hasName(args, "coorWise") & !is.null(args$coorWise)) {
          self$coorWise = args$coorWise

        } else {
          self$coorWise = TRUE

        }
      }

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
        return(list(costFunc = "L1",
                    coorWise = private$.params[["coorWise"]]))

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
