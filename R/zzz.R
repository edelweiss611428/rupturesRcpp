#' @title Cost_L2 module
#' @description Load Rcpp module for Cost_L2
#' @importFrom Rcpp loadModule
#' @import methods
#' @export
#' @name Cost_L2
loadModule("Cost_L2_module", TRUE)

#' @title Cost_VAR module
#' @description Load Rcpp module for Cost_VAR
#' @importFrom Rcpp loadModule
#' @import methods
#' @export
#' @name Cost_VAR
loadModule("Cost_VAR_module", TRUE)

#' @title Cost_SIGMA module
#' @description Load Rcpp module for Cost_SIGMA
#' @importFrom Rcpp loadModule
#' @import methods
#' @export
#' @name Cost_SIGMA
loadModule("Cost_SIGMA_module", TRUE)
