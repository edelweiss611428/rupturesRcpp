# nocov start
#' @importFrom methods loadMethod
NULL

.onLoad <- function(libname, pkgname) {

  #Cost modules
  Rcpp::loadModule("Cost_L2_module", TRUE)
  Rcpp::loadModule("Cost_VAR_module", TRUE)
  Rcpp::loadModule("Cost_SIGMA_module", TRUE)
  Rcpp::loadModule("Cost_L1_cwMed_module", TRUE)
  Rcpp::loadModule("Cost_LinearL2_module", TRUE)

  #Pelt
  Rcpp::loadModule("PELTCpp_L1_cwMed_module", TRUE)
  Rcpp::loadModule("PELTCpp_L2_module", TRUE)
  Rcpp::loadModule("PELTCpp_VAR_module", TRUE)
  Rcpp::loadModule("PELTCpp_SIGMA_module", TRUE)
  Rcpp::loadModule("PELTCpp_LinearL2_module", TRUE)

  #binSeg
  Rcpp::loadModule("binSegCpp_L1_cwMed_module", TRUE)
  Rcpp::loadModule("binSegCpp_L2_module", TRUE)
  Rcpp::loadModule("binSegCpp_VAR_module", TRUE)
  Rcpp::loadModule("binSegCpp_SIGMA_module", TRUE)
  Rcpp::loadModule("binSegCpp_LinearL2_module", TRUE)

  #Window
  Rcpp::loadModule("windowCpp_L1_cwMed_module", TRUE)
  Rcpp::loadModule("windowCpp_L2_module", TRUE)
  Rcpp::loadModule("windowCpp_VAR_module", TRUE)
  Rcpp::loadModule("windowCpp_SIGMA_module", TRUE)
  Rcpp::loadModule("windowCpp_LinearL2_module", TRUE)

}

.onAttach <- function(libname, pkgname) {

  packageStartupMessage(r"(

+--------------------------------------------------------------------+
|                                                                    |
|                    _                       ____                    |
|   _ __ _   _ _ __ | |_ _   _ _ __ ___  ___|  _ \ ___ _ __  _ __    |
|  | '__| | | | '_ \| __| | | | '__/ _ \/ __| |_) / __| '_ \| '_ \   |
|  | |  | |_| | |_) | |_| |_| | | |  __/\__ \  _ < (__| |_) | |_) |  |
|  |_|   \__,_| .__/ \__|\__,_|_|  \___||___/_| \_\___| .__/| .__/   |
|             |_|                                     |_|   |_|      |
|                                               version 1.0.0        |
+--------------------------------------------------------------------+

See https://github.com/edelweiss611428/rupturesRcpp/blob/main/README.md
for basic usage!
)") }
# nocov end
