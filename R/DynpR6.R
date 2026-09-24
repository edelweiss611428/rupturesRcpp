#' Exact Dynamic Programming (`Dynp`)
#'
#' @description An `R6` class implementing exact dynamic programming for offline change-point detection.
#'
#' @include costFuncR6.R
#' @docType class
#' @importFrom R6 R6Class is.R6
#' @import ggplot2
#' @import patchwork
#' @importFrom utils hasName
#' @export
#'
#' @details
#' `Dynp` finds the segmentation that globally minimises the total cost for a *specified* number of
#' change-points -- unlike `PELT` (exact for a penalty, not a change-point count) and `binSeg`
#' (greedy at each split, not guaranteed globally optimal for a fixed count), `Dynp` is exact for
#' whatever `nBkps` is requested. The price is complexity: building the full table costs
#' \eqn{O(\text{nBkpsMax} \cdot M^2)} where \eqn{M} is the number of `(minSize, jump)`-admissible
#' grid points (\eqn{M \approx n/\text{jump}} in the worst case), versus PELT's near-linear pruning
#' or binSeg's \eqn{O(n \log n)}-ish greedy search. `$fit()` computes the exact minimal cost for
#' *every* change-point count from `0` to `nBkpsMax` in one pass (see `$costPath()`), so this is
#' also the algorithm to reach for when choosing the number of change-points itself is the question
#' (the "elbow method"), not just finding them for a count you already know.
#'
#' `Dynp` requires a `R6` object of class `costFunc`, exactly as `PELT`/`binSeg`/`Window` do -- see
#' `costFunc` for the supported cost functions (`"L1"`, `"L2"`, `"SIGMA"`, `"VAR"`, `"LinearL2"`,
#' `"LinearSIGMA"`, `"LinearL1"`, `"Custom"`).
#'
#' Some examples are provided below. See the [GitHub README](https://github.com/edelweiss611428/rupturesRcpp/blob/main/README.md)
#' for detailed basic usage!
#'
#' @examples
#'
#' ## L2 example
#' set.seed(1121)
#' signals = as.matrix(c(rnorm(100,0,1),
#'                      rnorm(100,5,1)))
#' # Default L2 cost function; nBkpsMax left NULL -> resolved to min(20, maximum feasible value)
#' DynpObj = Dynp$new(minSize = 1L, jump = 1L)
#' DynpObj$fit(signals)
#' DynpObj$predict(nBkps = 1)
#' DynpObj$plot()
#'
#' # Same result, reached via a penalty instead of an explicit count
#' DynpObj$predict(pen = 100)
#'
#' # The exact cost of every change-point count from 0 to nBkpsMax -- elbow-method data
#' DynpObj$getHistory()
#' DynpObj$plotElbow()
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `Dynp` object.}
#'   \item{\code{$describe()}}{Describes the `Dynp` object.}
#'   \item{\code{$fit()}}{Constructs a `Dynp` module in `C++`.}
#'   \item{\code{$eval()}}{Evaluates the cost of a segment.}
#'   \item{\code{$predict()}}{Returns the exact optimal breakpoints, for a given `pen` or `nBkps`.}
#'   \item{\code{$costPath()}}{Returns the exact minimal cost for every change-point count up to `nBkpsMax`.}
#'   \item{\code{$getHistory()}}{Same data as `$costPath()`, as a `k`/`cost` `data.frame` (same shape as `binSeg`/`Window`'s `$getHistory()`, minus `added_bkp` -- see its own docs for why).}
#'   \item{\code{$plotElbow()}}{Plots the elbow curve (Total Cost vs. Number of Change-Points).}
#'   \item{\code{$plot()}}{Plots change-point segmentation in `ggplot` style.}
#'   \item{\code{$clone()}}{Clones the `R6` object.}
#' }
#'
#' @references
#' Truong, C., Oudre, L., & Vayatis, N. (2020). Selective review of offline change point detection methods.
#' Signal Processing, 167, 107299.
#'
#' @author
#' Minh Long Nguyen \email{edelweiss611428@gmail.com} \cr
#' Toby Dylan Hocking \email{toby.hocking@r-project.org} \cr
#' Charles Truong \email{ctruong@ens-paris-saclay.fr}
#'
#' @export

Dynp = R6Class(
  "Dynp",

  private = list(

    .minSize = 1L,
    .jump = 1L,
    .nBkpsMax = NULL,           # user setting: NULL (auto) or a positive integer
    .resolvedNBkpsMax = NULL,   # value actually used by the last $fit()
    .costFunc = costFunc$new("L2"), #L2 cost function
    .DynpModule = NULL,
    .tsMat = NULL,
    .covariates = NULL,
    .fitted = FALSE,
    .n = NULL,
    .p = NULL,
    .tmpEndPts = NULL, #Temporary end points
    .tmpPen = NULL, #Temporary penalty value
    .tmpNBkps = NULL #Temporary nBkps value

  ),

  active = list(

    #' @field minSize Integer. Minimum allowed segment length. Can be accessed or modified via `$minSize`.
    #' Modifying `minSize` will automatically trigger `$fit()`.
    minSize = function(intVal) {

      if(missing(intVal)){
        return(private$.minSize)
      }

      if(is.null(intVal)){
        stop("'minSize' must not be null!")
      }

      if (!is.numeric(intVal) | any(intVal < 1) | length(intVal) != 1) {
        stop("'minSize' must be a single positive integer!")
      }

      private$.minSize = as.integer(intVal)

      if (!is.null(private$.tsMat) & private$.fitted) {
        message("`minSize` has been updated. Re-fitting the model.")
        self$fit()
      }

    },


    #' @field jump Integer. Search grid step size. Can be accessed or modified via `$jump`.
    #' Modifying `jump` will automatically trigger `$fit()`.
    jump = function(intVal) {

      if(missing(intVal)){
        return(private$.jump)
      }

      if(is.null(intVal)){
        stop("'jump' must not be null!")
      }

      if (!is.numeric(intVal) | any(intVal < 1) | length(intVal) != 1) {
        stop("'jump' must be a single positive integer!")
      }

      private$.jump = as.integer(intVal)

      if (!is.null(private$.tsMat) & private$.fitted) {
        message("`jump` has been updated. Re-fitting the model.")
        self$fit()
      }

    },

    #' @field nBkpsMax Integer or `NULL`. Upper bound on the number of change-points the exact DP
    #' table is built for; `$predict(nBkps = ...)` treats this as a cap, not a requirement -- see
    #' `$predict()`. If `NULL` (default), `$fit()` resolves it to
    #' \code{min(20, floor(n / minSize) - 1)}, i.e. the maximum feasible value given `minSize`, capped
    #' at 20. Set explicitly to go beyond 20 (up to the maximum feasible value) or to lower it; the
    #' table costs \eqn{O(\text{nBkpsMax} \cdot M^2)} work and \eqn{O(\text{nBkpsMax} \cdot M)} memory.
    #' Can be accessed or modified via `$nBkpsMax`; modifying it will automatically trigger `$fit()`.
    nBkpsMax = function(intVal) {

      if(missing(intVal)){
        return(private$.nBkpsMax)
      }

      if(!is.null(intVal)){
        if (!is.numeric(intVal) | any(intVal < 0) | length(intVal) != 1) {
          stop("'nBkpsMax' must be NULL or a single non-negative integer!")
        }
        intVal = as.integer(intVal)
      }

      private$.nBkpsMax = intVal

      if (!is.null(private$.tsMat) & private$.fitted) {
        message("`nBkpsMax` has been updated. Re-fitting the model.")
        self$fit()
      }

    },

    #' @field costFunc `R6` object of class `costFunc`. Can be accessed or modified via `$costFunc`.
    #' Modifying `costFunc` will automatically trigger `$fit()`.
    costFunc = function(costFuncObj) {

      if(missing(costFuncObj)){
        return(private$.costFunc)
      }

      if (!inherits(costFuncObj, "costFunc") | !is.R6(costFuncObj)) {
        stop("`costFunc` must be a `R6` object of class `costFunc` - can be created via costFunc$new()!")
      }

      private$.costFunc = costFuncObj

      if (!is.null(private$.tsMat) & private$.fitted) {
        message("`costFunc` has been updated. Re-fitting the model.")
        self$fit()
      }
    },

    #' @field tsMat Numeric matrix. Input time series matrix of size \eqn{n \times p}. Can be accessed or modified via `$tsMat`.
    #' Modifying `tsMat` will automatically trigger `$fit()`.
    #'
    tsMat = function(numMat) {

      if(missing(numMat)){
        return(private$.tsMat)
      }

      if(is.null(numMat)){
        stop("`tsMat` must not be null!")
      }

      if (!is.numeric(numMat) | !is.matrix(numMat)) {
        stop("`tsMat` must be a numeric time series matrix!")
      }

      if(any(is.na(numMat))){
        stop("`tsMat` contains NAs!")
      }

      private$.tsMat = numMat
      private$.n = nrow(numMat)
      private$.p = ncol(numMat)

      if (private$.fitted) {
        self$fit()
      }

    },

    #' @field covariates Numeric matrix. Input time series matrix having a similar number of observations as `tsMat`. Can be accessed or modified via `$covariates`.
    #' Modifying `covariates` will automatically trigger `$fit()`.
    #'
    covariates = function(numMat) {

      if(missing(numMat)){
        return(private$.covariates)
      }

      if(is.null(numMat)){
        stop("`covariates` must not be null!")
      }

      if (!is.numeric(numMat) | !is.matrix(numMat)) {
        stop("`covariates` must be a numeric time series matrix!")
      }

      if(any(is.na(numMat))){
        stop("`covariates` contains NAs!")
      }

      private$.covariates = numMat

      if(private$.costFunc$pass()[["costFunc"]] %in% c("LinearL2", "LinearSIGMA", "LinearL1")){
        if (!is.null(private$.tsMat) & private$.fitted) {
          self$fit()
        }
      }

    }


  ),

  public = list(

    #' @description Initialises a `Dynp` object.
    #'
    #' @param minSize Integer. Minimum allowed segment length. Default: `1L`.
    #' @param jump Integer. Search grid step size: only positions in \{k, 2k, ...\} are considered. Default: `1L`.
    #' @param nBkpsMax Integer or `NULL`. Upper bound on the number of change-points to build the exact
    #' DP table for. Default: `NULL` (resolved to \code{min(20, floor(n / minSize) - 1)} at `$fit()` time).
    #' @param costFunc A `R6` object of class `costFunc`. Should be created via `costFunc$new()` to avoid error.
    #' Default: `costFunc$new("L2")`.
    #'
    #' @return Invisibly returns `NULL`.

    initialize = function(minSize, jump, nBkpsMax, costFunc) {

      if(!missing(minSize)){
        self$minSize = minSize
      }

      if(!missing(jump)){
        self$jump = jump
      }

      if(!missing(nBkpsMax)){
        self$nBkpsMax = nBkpsMax
      }

      if(!missing(costFunc)){
        self$costFunc = costFunc
      }

      invisible(NULL)
    },

    #' @description Describes a `Dynp` object.
    #'
    #' @param printConfig Logical. Whether to print object configurations. Default: `FALSE`.
    #'
    #' @return Invisibly returns a list storing at least the following fields:
    #'
    #' \describe{
    #'   \item{\code{minSize}}{Minimum allowed segment length.}
    #'   \item{\code{jump}}{Search grid step size.}
    #'   \item{\code{nBkpsMax}}{The user-set `nBkpsMax` (possibly `NULL`).}
    #'   \item{\code{resolvedNBkpsMax}}{The `nBkpsMax` actually used by the last `$fit()` (`NULL` if not fitted).}
    #'   \item{\code{costFunc}}{The `costFunc` object.}
    #'   \item{\code{fitted}}{Whether or not `$fit()` has been run.}
    #'   \item{\code{tsMat}}{Time series matrix.}
    #'   \item{\code{covariates}}{Covariate matrix (if exists).}
    #'   \item{\code{n}}{Number of observations.}
    #'   \item{\code{p}}{Number of features.}
    #' }

    describe = function(printConfig = FALSE) {

      if(is.null(printConfig)){
        stop("`printConfig` is null!")
      } else{
        if(!is.logical(printConfig) | length(printConfig) != 1L){
          stop("`printConfig` must be a single boolean value!")
        }
      }

      params = list(minSize = private$.minSize,
                    jump = private$.jump,
                    nBkpsMax = private$.nBkpsMax,
                    resolvedNBkpsMax = private$.resolvedNBkpsMax,
                    costFunc = private$.costFunc,
                    fitted = private$.fitted,
                    tsMat = private$.tsMat,
                    covariates = private$.covariates,
                    n = private$.n,
                    p = private$.p)

      if(printConfig){

        cat(sprintf("Exact Dynamic Programming (Dynp) \n"))
        cat(sprintf("minSize          : %sL\n", private$.minSize))
        cat(sprintf("jump             : %sL\n", private$.jump))
        cat(sprintf("nBkpsMax         : %s\n", if(is.null(private$.nBkpsMax)) "NULL (auto)" else paste0(private$.nBkpsMax, "L")))
        cat(sprintf("resolvedNBkpsMax : %s\n", if(is.null(private$.resolvedNBkpsMax)) "NULL" else paste0(private$.resolvedNBkpsMax, "L")))
        cat(sprintf("costFunc         : \"%s\"\n", private$.costFunc$pass()[["costFunc"]]))

      }


      if(private$.costFunc$pass()[["costFunc"]] == "SIGMA"){

        if(printConfig){

          cat(sprintf("addSmallDiag     : %s\n", private$.costFunc$pass()[["addSmallDiag"]]))
          cat(sprintf("epsilon          : %s\n", private$.costFunc$pass()[["epsilon"]]))

        }

        params[["addSmallDiag"]] = private$.costFunc$pass()[["addSmallDiag"]]
        params[["epsilon"]] = private$.costFunc$pass()[["epsilon"]]

      }

      if(private$.costFunc$pass()[["costFunc"]] == "VAR"){

        if(printConfig){

          cat(sprintf("pVAR             : %sL\n", private$.costFunc$pass()[["pVAR"]]))

        }

        params[["pVAR"]] = private$.costFunc$pass()[["pVAR"]]

      }

      if(private$.costFunc$pass()[["costFunc"]] == "LinearL2"){

        if(printConfig){

          cat(sprintf("intercept        : %sL\n", private$.costFunc$pass()[["intercept"]]))

        }

        params[["intercept"]] = private$.costFunc$pass()[["intercept"]]

      }

      if(private$.costFunc$pass()[["costFunc"]] == "LinearSIGMA"){

        if(printConfig){

          cat(sprintf("intercept        : %sL\n", private$.costFunc$pass()[["intercept"]]))
          cat(sprintf("addSmallDiag     : %s\n", private$.costFunc$pass()[["addSmallDiag"]]))
          cat(sprintf("epsilon          : %s\n", private$.costFunc$pass()[["epsilon"]]))

        }

        params[["intercept"]] = private$.costFunc$pass()[["intercept"]]
        params[["addSmallDiag"]] = private$.costFunc$pass()[["addSmallDiag"]]
        params[["epsilon"]] = private$.costFunc$pass()[["epsilon"]]

      }

      if(private$.costFunc$pass()[["costFunc"]] == "LinearL1"){

        if(printConfig){

          cat(sprintf("intercept        : %sL\n", private$.costFunc$pass()[["intercept"]]))
          cat(sprintf("tol              : %s\n", private$.costFunc$pass()[["tol"]]))
          cat(sprintf("maxIter          : %sL\n", private$.costFunc$pass()[["maxIter"]]))

        }

        params[["intercept"]] = private$.costFunc$pass()[["intercept"]]
        params[["tol"]] = private$.costFunc$pass()[["tol"]]
        params[["maxIter"]] = private$.costFunc$pass()[["maxIter"]]

      }

      if(private$.costFunc$pass()[["costFunc"]] == "Custom"){

        if(printConfig){

          cat(sprintf("evalFun          : <function>\n"))
          cat(sprintf("paramFun         : %s\n", if(is.null(private$.costFunc$pass()[["paramFun"]])) "NULL" else "<function>"))

        }

        params[["evalFun"]] = private$.costFunc$pass()[["evalFun"]]
        params[["paramFun"]] = private$.costFunc$pass()[["paramFun"]]

      }

      if(printConfig){

        cat(sprintf("fitted           : %s\n", private$.fitted))
        cat(sprintf("n                : %sL\n", private$.n))
        cat(sprintf("p                : %sL\n", private$.p))

      }

      invisible(params)

    },

    #' @description Constructs a `C++` module for `Dynp` and builds the exact DP table.
    #'
    #' @param tsMat Numeric matrix. A time series matrix of size \eqn{n \times p} whose rows are observations ordered in time.
    #' If `tsMat = NULL`, the method will use the previously assigned `tsMat` (e.g., set via the active binding `$tsMat`
    #' or from a prior `$fit(tsMat)`). Default: `NULL`.
    #'
    #' @param covariates Numeric matrix. A time series matrix having a similar number of observations as `tsMat`.
    #' Required for models involving both dependent and independent variables.
    #' If `covariates = NULL` and no prior covariates were set (i.e., `$covariates` is still `NULL`),
    #' the model is force-fitted with only an intercept. Default: `NULL`.
    #'
    #' @return Invisibly returns `NULL`.
    #'
    #' @details This method constructs a `C++` `Dynp` module and sets `private$.fitted` to `TRUE`, enabling the
    #' use of `$predict()`, `$costPath()` and `$eval()`. If `$nBkpsMax` is `NULL`, it is resolved here to
    #' \code{min(20, floor(n / minSize) - 1)} (with a message); if set higher than \code{floor(n / minSize) - 1},
    #' it is capped to that (with a warning).

    fit = function(tsMat = NULL, covariates = NULL) {

      if (!is.null(tsMat)) {

        if (!is.numeric(tsMat) | !is.matrix(tsMat)) {
          stop("`tsMat` must be a numeric time series matrix!")
        }

        if(any(is.na(tsMat))){
          stop("`tsMat` contains NAs!")
        }

        private$.tsMat = tsMat
        private$.n = nrow(tsMat)
        private$.p = ncol(tsMat)
      } else{

        if (is.null(private$.tsMat)) {
          stop("No `tsMat` found! Please provide a `tsMat`!")
        }

      }

      naturalMax = private$.n %/% private$.minSize - 1L

      if(is.null(private$.nBkpsMax)){
        resolvedNBkpsMax = min(naturalMax, 20L)
        message(sprintf("`nBkpsMax` not set; using %dL.", resolvedNBkpsMax))
      } else if(private$.nBkpsMax > naturalMax){
        resolvedNBkpsMax = naturalMax
        warning(sprintf(
          "`nBkpsMax` (%dL) exceeds the maximum feasible value given `minSize` (%dL); using %dL instead. This may be slow/memory-intensive for large series.",
          private$.nBkpsMax, naturalMax, naturalMax))
      } else{
        resolvedNBkpsMax = private$.nBkpsMax
      }

      private$.resolvedNBkpsMax = resolvedNBkpsMax

      if(private$.costFunc$pass()[["costFunc"]] %in% c("LinearL2", "LinearSIGMA", "LinearL1")){

        if(!is.null(covariates)){

          if (!is.numeric(covariates) | !is.matrix(covariates)) {
            stop("`covariates` must be a numeric time series matrix!")
          }

          if(any(is.na(covariates))){
            stop("`covariates` contains NAs!")
          }

          if(nrow(covariates) != private$.n){
            stop("Numbers of observations in `covariates` and `tsMat` do not match!")
          }

          private$.covariates = covariates

        } else{

          if(is.null(private$.covariates)){
            warning("No `covariates` found! Force-fitting with only an intercept!")

            if(private$.costFunc$pass()[["costFunc"]] == "LinearL2"){

              private$.DynpModule = new(DynpCpp_LinearL2, private$.tsMat, matrix(1, nrow = private$.n, ncol = 1),
                                        FALSE, #no intercept
                                        private$.minSize, private$.jump, resolvedNBkpsMax)

            } else if(private$.costFunc$pass()[["costFunc"]] == "LinearSIGMA"){

              private$.DynpModule = new(DynpCpp_LinearSIGMA, private$.tsMat, matrix(1, nrow = private$.n, ncol = 1),
                                        FALSE, #no intercept
                                        private$.costFunc$pass()[["addSmallDiag"]],
                                        private$.costFunc$pass()[["epsilon"]],
                                        private$.minSize, private$.jump, resolvedNBkpsMax)

            } else if(private$.costFunc$pass()[["costFunc"]] == "LinearL1"){

              private$.DynpModule = new(DynpCpp_LinearL1, private$.tsMat, matrix(1, nrow = private$.n, ncol = 1),
                                        FALSE, #no intercept
                                        private$.costFunc$pass()[["tol"]],
                                        private$.costFunc$pass()[["maxIter"]],
                                        private$.minSize, private$.jump, resolvedNBkpsMax)

            }

            private$.fitted = TRUE
            private$.DynpModule$fit()

            return(invisible(NULL))


          } else{

            if(nrow(private$.covariates) != private$.n){
              stop("Numbers of observations in `covariates` and `tsMat` do not match!")
            }

          }
        }
      }


      if(private$.costFunc$pass()[["costFunc"]] == "L2"){
        private$.DynpModule = new(DynpCpp_L2, private$.tsMat, private$.minSize, private$.jump, resolvedNBkpsMax)

      } else if(private$.costFunc$pass()[["costFunc"]] == "VAR"){
        private$.DynpModule = new(DynpCpp_VAR, private$.tsMat, private$.costFunc$pass()[["pVAR"]],
                                  private$.minSize, private$.jump, resolvedNBkpsMax)

      } else if(private$.costFunc$pass()[["costFunc"]] == "SIGMA"){
        private$.DynpModule = new(DynpCpp_SIGMA, private$.tsMat,
                                  private$.costFunc$pass()[["addSmallDiag"]],
                                  private$.costFunc$pass()[["epsilon"]],
                                  private$.minSize, private$.jump, resolvedNBkpsMax)

      } else if(private$.costFunc$pass()[["costFunc"]] == "L1"){

        private$.DynpModule = new(DynpCpp_L1_cwMed, private$.tsMat,
                                  private$.minSize, private$.jump, resolvedNBkpsMax)

      } else if(private$.costFunc$pass()[["costFunc"]] == "LinearL2"){

        private$.DynpModule = new(DynpCpp_LinearL2, private$.tsMat, private$.covariates,
                                  private$.costFunc$pass()[["intercept"]],
                                  private$.minSize, private$.jump, resolvedNBkpsMax)

      } else if(private$.costFunc$pass()[["costFunc"]] == "LinearSIGMA"){

        private$.DynpModule = new(DynpCpp_LinearSIGMA, private$.tsMat, private$.covariates,
                                  private$.costFunc$pass()[["intercept"]],
                                  private$.costFunc$pass()[["addSmallDiag"]],
                                  private$.costFunc$pass()[["epsilon"]],
                                  private$.minSize, private$.jump, resolvedNBkpsMax)

      } else if(private$.costFunc$pass()[["costFunc"]] == "LinearL1"){

        private$.DynpModule = new(DynpCpp_LinearL1, private$.tsMat, private$.covariates,
                                  private$.costFunc$pass()[["intercept"]],
                                  private$.costFunc$pass()[["tol"]],
                                  private$.costFunc$pass()[["maxIter"]],
                                  private$.minSize, private$.jump, resolvedNBkpsMax)

      } else if(private$.costFunc$pass()[["costFunc"]] == "Custom"){

        private$.DynpModule = new(DynpCpp_RFunc, private$.tsMat,
                                  private$.costFunc$pass()[["evalFun"]],
                                  private$.costFunc$pass()[["paramFun"]],
                                  private$.minSize, private$.jump, resolvedNBkpsMax)

      } else{
        # nocov start
        stop("Cost function not supported!")
        # nocov end
      }

      private$.fitted = TRUE

      private$.DynpModule$fit()

      invisible(NULL)
    },

    #' @description Evaluate the cost of the segment (a,b]
    #'
    #' @param a Integer. Start index of the segment (exclusive). Must satisfy \code{start < end}.
    #' @param b Integer. End index of the segment (inclusive).
    #'
    #' @return The segment cost. See `costFunc` for the cost formulas.

    eval = function(a, b){

      if(!private$.fitted){
        stop("$fit() must be run before $eval()!")
      }

      if(is.null(a) | is.null(b)){
        stop("`a` and `b` must not be NULL")
      }

      if (!is.numeric(a) | any(a < 0) | length(a) != 1 | any(a > private$.n)) {
        stop("`0 <= start < nSamples` must be true!")
      }

      if (!is.numeric(b) | any(b < 0) | length(b) != 1 | any(b>private$.n)) {
        stop("`0 < end <= nSamples` must be true!")
      }

      a = as.integer(a)
      b = as.integer(b)

      if(a >= b){
        stop("a must be smaller than b!")
      }

      return(private$.DynpModule$eval(a, b))

    },

    #' @description Returns the exact optimal breakpoints, either under a linear penalty or for a
    #' specified number of change-points.
    #'
    #' @param pen Numeric. Penalty per change-point; the change-point count is chosen by minimising
    #' \eqn{\text{cost} + \text{pen} \cdot k} over \eqn{k = 0, \dots, \text{nBkpsMax}} (same convention
    #' as `PELT`/`binSeg`/`Window`). Ignored if `nBkps` is supplied. Default: `0`.
    #' @param nBkps Integer. If supplied, takes precedence over `pen`: returns the exact optimal
    #' segmentation for exactly `nBkps` change-points. Treated as an upper bound, not a strict
    #' requirement, in the same spirit as `binSeg`/`Window`: if `nBkps` exceeds `nBkpsMax`, it is
    #' capped to `nBkpsMax` and a message reports the shortfall. Unlike `binSeg`/`Window`, the result
    #' is still exact for whatever count is actually used -- capping only ever happens because the DP
    #' table wasn't built that far (a configuration choice via `$nBkpsMax`), not because `Dynp` ran
    #' out of candidates the way a greedy search can. If no valid segmentation exists at all for the
    #' (possibly capped) count -- only possible when `jump` is large relative to `minSize` -- this
    #' still errors, since there is no well-defined smaller-but-still-what-you-asked-for answer to
    #' fall back to. Default: `NULL`.
    #'
    #' @return An integer vector of regime end-points. By design, the last element is the
    #' number of observations.
    #'
    #' @details
    #' With `nBkps = k`, this returns the segmentation into exactly `k` change-points that globally
    #' minimises the total cost -- exact, unlike `binSeg$predict()` for the same `k`. With `pen`,
    #' this instead minimises the penalised cost over every change-point count the DP table covers,
    #' matching how `PELT`/`binSeg`/`Window` already select a count from a penalty.
    #'
    #' Both the DP table (via `$fit()`) and the traceback here only depend on `minSize`/`jump`
    #' admissible positions, so `$predict()` itself is cheap (`O(nBkps)`) regardless of which mode
    #' is used -- the cost was already paid in `$fit()`.
    #'
    #' Temporary segment end-points are saved to `private$.tmpEndPoints` after `$predict()`, enabling
    #' users to call `$plot()` without specifying endpoints manually.

    predict = function(pen = 0, nBkps = NULL){

      if(!private$.fitted){
        stop("`$fit()` must be run before `$predict()`!")
      }

      if(!is.null(nBkps)){

        if(!is.numeric(nBkps) | length(nBkps) != 1){
          stop("`nBkps` must be a single non-negative integer!")
        }

        if(any(nBkps < 0) | any(nBkps != round(nBkps))){
          stop("`nBkps` must be a single non-negative integer!")
        }

        nBkps = as.integer(nBkps)
        kUse = min(nBkps, private$.resolvedNBkpsMax)

        if(kUse < nBkps){
          message(sprintf(
            "`nBkps` (%dL) exceeds `nBkpsMax` (%dL); using %dL instead.",
            nBkps, private$.resolvedNBkpsMax, kUse))
        }

        endPts = c(private$.DynpModule$predictK(kUse))
        private$.tmpNBkps = nBkps
        private$.tmpPen = NULL

      } else{

        if(is.null(pen)){
          stop("`pen` must not be null!")
        }

        if(!is.numeric(pen) | length(pen) != 1 | any(pen < 0)){
          stop("`pen` must be a single non-negative value!")
        }

        endPts = c(private$.DynpModule$predictPen(pen))
        private$.tmpPen = pen
        private$.tmpNBkps = NULL

      }

      private$.tmpEndPts = endPts

      return(endPts)

    },

    #' @description Returns the exact minimal cost for every change-point count the DP table covers.
    #'
    #' @return A numeric vector of length `resolvedNBkpsMax + 1`: element `k+1` is the exact minimal
    #' total cost of segmenting the series into exactly `k` change-points, for `k = 0, ..., resolvedNBkpsMax`.
    #' `Inf` at position `k+1` means no valid segmentation with exactly `k` change-points exists given
    #' `minSize`/`jump` (only possible when `jump` is large relative to `minSize`).
    #'
    #' @details This is the data behind the "elbow method" for choosing the number of change-points:
    #' plot it (see `$getHistory()`/`$plotElbow()`) and look for where the marginal decrease in cost
    #' flattens out. It comes for free out of `$fit()` -- no extra computation is triggered here.

    costPath = function(){

      if(!private$.fitted){
        stop("`$fit()` must be run before `$costPath()`!")
      }

      return(private$.DynpModule$costPath())

    },

    #' @description Retrieves the exact minimal cost for every change-point count the DP table covers,
    #' as a `data.frame` (same shape as `binSeg`/`Window`'s `$getHistory()`, minus `added_bkp`).
    #'
    #' @return A `data.frame` with two columns:
    #' \describe{
    #'   \item{\code{k}}{The number of change-points, from `0` to the resolved `nBkpsMax`.}
    #'   \item{\code{cost}}{The exact minimal total cost of segmenting the series into exactly `k`
    #'   change-points (see `$costPath()`).}
    #' }
    #'
    #' @details Unlike `binSeg`/`Window`, there is no `added_bkp` column here: `binSeg` and `Window`
    #' each build a single *nested* sequence of change-points, where every `k`'s answer is the
    #' previous `k-1`'s answer plus one more point, so "the breakpoint added at this step" is
    #' well-defined. `Dynp`'s per-`k` solutions are each independently exact and need not be nested at
    #' all -- the optimal segmentation for `k` change-points can differ completely from the one for
    #' `k-1`. Use `$predict(nBkps = k)` to get the full breakpoint set for a given `k`.

    getHistory = function(){

      cost_vec = self$costPath()

      return(data.frame(k = 0:(length(cost_vec) - 1), cost = cost_vec))

    },

    #' @description Plots the elbow curve (Total Cost vs. Number of Change-Points).
    #'
    #' @param maxK Integer. The maximum number of change-points to display on the plot.
    #' If `NULL`, displays the full history. Default: `NULL`.
    #'
    #' @return A `ggplot` object.

    plotElbow = function(maxK = NULL) {

      if(!private$.fitted){
        stop("`$fit()` must be run before `$plotElbow()`!")
      }

      hist_df = self$getHistory()

      if(!is.null(maxK)) {
        if(!is.numeric(maxK) || maxK < 1) stop("`maxK` must be a positive integer!")
        hist_df = hist_df[hist_df$k <= maxK, ]
      }

      p = ggplot(hist_df, aes(x = k, y = cost)) +
        geom_line(color = "#5B9BD5", linewidth = 0.8) +
        geom_point(color = "#5B9BD5", size = 2) +
        scale_x_continuous(breaks = hist_df$k) +
        theme_minimal() +
        theme(panel.grid.minor.x = element_blank()) +
        labs(
          title = "Exact Dynamic Programming Elbow Plot",
          subtitle = paste("Cost Function:", private$.costFunc$pass()[["costFunc"]]),
          x = "Number of Change-Points (k)",
          y = "Total Cost"
        )

      return(p)
    },


    #' @description Plots change-point segmentation
    #'
    #' @param d Integer vector. Dimensions to plot. Default: `1L`.
    #' @param endPts Integer vector. End points. Default: latest temporary changepoints obtained via `$predict()`.
    #' @param dimNames Character vector. Feature names matching length of `d`. Defaults to `"X1", "X2", ...`.
    #' @param main Character. Main title. Defaults to `"Dynp: d = ..."`.
    #' @param xlab Character. X-axis label. Default: `"Time"`.
    #' @param tsWidth Numeric. Line width for time series and segments. Default: `0.25`.
    #' @param tsCol Character. Time series color. Default: `"#5B9BD5"`.
    #' @param bgCol Character vector. Segment colors, recycled to length of `endPts`. Default: `c("#A3C4F3", "#FBB1BD")`.
    #' @param bgAlpha Numeric. Background transparency. Default: `0.5`.
    #' @param ncol Integer. Number of columns in facet layout. Default: `1L`.
    #'
    #' @details Plots change-point segmentation results. Based on `ggplot2`. Multiple plots can easily be
    #' horizontally and vertically stacked using `patchwork`'s operators `/` and `|`, respectively.
    #'
    #' @return An object of classes `gg` and `ggplot`.

    plot = function(d = 1L, endPts, dimNames, main, xlab, tsWidth = 0.25,
                    tsCol = "#5B9BD5",
                    bgCol = c("#A3C4F3", "#FBB1BD"),
                    bgAlpha = 0.5,
                    ncol = 1L){


      if(missing(main)){
        main = paste0("Dynp: ", paste0("d = (", toString(d), ")"))

      } else {
        if(!is.character(main) | length(main )!= 1L){
          stop("`main` must be a single character!")
        }
      }

      if(missing(xlab)){
        xlab = "Time"

      } else {
        if(!is.character(xlab) | length(xlab)!= 1L){
          stop("`xlab` must be a single character!")
        }
      }

      if(missing(endPts)){
        message("`endPts` is missing. Proceed to use the temporary `endPts`!")

        if(is.null(private$.tmpEndPts)){
          stop("Temporary `endPts` is null. Must run `$predict()` to initialise this!")
        } else{
          endPts = private$.tmpEndPts
        }

      } else{
        if(!is.numeric(endPts)){
          stop("`endPts` must be an integer vector specifying endpoints! Could be obtained via `$predict()`.")
        } else{
          endPts = as.integer(sort(endPts))

          if(min(endPts) < 1){
            stop("`min(endPts)` must not be less than 1!")
          }

          if(max(endPts) != private$.n){
            stop("By construction, `max(endPts)` must be `n`! Can use `$predict()` to obtain `endPts`.")
          }

          if(length(unique(endPts)) != length(endPts)){
            stop("`as.integer(endPts)` contains duplicated elements!")
          }
        }
      }

      if (!is.numeric(d)) {
        stop("`d` must be a numeric/integer vector specifying dimensions! e.g., `1`, or `c(1,2)`.")
      }

      if(missing(dimNames)){
        message("`dimNames` is missing. Proceed to use the default `dimNames`! e.g., `paste0('X', d))`.")
        dimNames = paste0("X", d)
      } else {
        if (!is.character(dimNames)) {
          stop("`dimNames` must be a character vector specifying feature names!")
        }
        if(length(dimNames) != length(d)){
          stop("length(dimNames) != length(d)!")
        }
      }

      d = as.integer(d)
      if(any(d < 1) | any(d > private$.p)){
        stop("Dimensions must be between [1,p]!")
      }

      tsList = lapply(seq_along(d), function(i) {
        data.frame(
          time = 1:private$.n,
          value = private$.tsMat[, d[i]],
          dimension = dimNames[i]
        )
      })

      allTsDf = do.call(rbind, tsList)

      segInt = data.frame(
        xmin = c(1, endPts[-length(endPts)]),
        xmax = endPts,
        fill = rep(bgCol, length.out = length(endPts))
      )

      ggplot(allTsDf, aes(x = time, y = value)) +
        scale_fill_identity() +
        geom_rect(
          data = segInt,
          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
          inherit.aes = FALSE,
          alpha = bgAlpha
        ) +
        geom_line(color = tsCol, linewidth = tsWidth) +
        geom_vline(xintercept = endPts[-length(endPts)], linetype = "dashed", color = "black", linewidth = tsWidth) +
        facet_wrap(~ dimension, scales = "free_y", ncol = ncol) +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold")
        ) +
        labs(x = xlab, y = NULL, title = main)
    }

  )

)
