#' Create a costFunc object
#'
#' @description Constructs a cost function object of class \code{costFunc} based on the specified
#' cost function type and optional parameters.
#'
#' @param costFunc Character. Cost function. Supported values include \code{"L2"}, \code{"VAR"},
#' and \code{"SIGMA"}. Default `L2`.
#' @param ... Optional named parameters required by specific cost functions.
#'
#' For \code{"SIGMA"}, supported parameters are:
#' \describe{
#'   \item{`addSmallDiag`}{Logical. If \code{TRUE}, add a small value to the diagonal of estimated covariance matrices
#'   to improve numerical stability. Default: \code{TRUE}.}
#'   \item{`epsilon`}{Double. A small positive value added to the diagonal of estimated covariance matrices to stabilise
#'   matrix operations. Default: \code{1e-6}.}
#' }
#'
#' For \code{"VAR"}, the VAR order \code{pVAR} is required:
#' \describe{
#'   \item{`pVAR`}{Integer. Vector autoregressive order. Must be non-negative. Default: \code{1L}.}
#' }
#'
#' @return An object of class \code{costFunc} - a list containing costFunc and any additional parameters.
#'
#' @details Creates a `costFunc` object, required in segmentation methods.
#'
#' Supported cost functions include:
#'
#' - `"L2"`: for (independent) piecewise Gaussian process with **constant variance**
#' - `"SIGMA"`: for (independent) piecewise Gaussian process with **varying variance**
#' - `"VAR"`: for piecewise Gaussian vector-regressive process with **constant variance**
#'
#' @examples
#' createCostFunc("L2")
#' createCostFunc("SIGMA", addSmallDiag = TRUE, epsilon = 1e-6)
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

createCostFunc = function(costFunc = "L2", ...){

  if(is.null(costFunc) | !is.character(costFunc) | length(costFunc) != 1L){
    stop("costFunc must be a single character value!")
  }

  args = list(...)

  if (costFunc %in% c("L2")) { #Cost functions requiring no additional param
    costFuncObj = list(costFunc = costFunc)
    message(paste0(costFunc, " cost function requires no additional param!"))
  } else if (costFunc == "SIGMA") {
    message("SIGMA cost function requires `addSmallDiag` and `epsilon`!")
    costFuncObj = list(costFunc = "SIGMA")

    if (is.null(args$addSmallDiag)) {
      message("No `addSmallDiag` provided for SIGMA cost function!")
      message("Use default option `addSmallDiag = TRUE`.")
      costFuncObj[["addSmallDiag"]] = TRUE

    } else{
      if(!is.logical(args$addSmallDiag) | length(args$addSmallDiag) != 1){
        stop("`addSmallDiag` must be a single boolean (T/F) value!")

      }
      costFuncObj[["addSmallDiag"]] = args$addSmallDiag

    }

    if(is.null(args$epsilon)){

      message("No `epsilon` provided for SIGMA cost function!")
      message("Use default option `epsilon = 1e-6`.")
      costFuncObj[["epsilon"]] = 1e-6

    }  else{
      if(!is.numeric(args$epsilon) | length(args$epsilon) != 1L | any(args$epsilon <= 0)){
        stop("`epsilon` must be a single positive numeric value (ideally small)!")

      }
      costFuncObj[["epsilon"]] = args$epsilon

    }
  } else if (costFunc == "VAR"){

    message("VAR cost function requires `pVAR`!")
    costFuncObj = list(costFunc = "VAR")

    if(is.null(args$pVAR)){

      message("No `pVAR` provided for VAR cost function!")
      message("Use default option `pVAR = 1L`.")
      costFuncObj[["pVAR"]] = 1L

    } else {

      if(!is.numeric(args$pVAR) | length(args$pVAR) != 1L | any(as.integer(args$pVAR) < 1)){
        stop("`pVAR` must be a single positive integer!")

      }
      costFuncObj[["pVAR"]] = as.integer(args$pVAR)
    }

  } else{
    stop("Cost function not supported!")
  }

  class(costFuncObj) = "costFunc"

  return(costFuncObj)
}


