#' Binary Segmentation (`binSeg`)
#'
#' @description An `R6` class implementing binary segmentation for offline change-point detection.
#'
#' @include costFuncR6.R
#' @docType class
#' @importFrom R6 R6Class is.R6
#' @importFrom ggplot2 aes ggplot geom_rect geom_line scale_fill_identity theme_minimal theme geom_vline labs element_blank element_text facet_wrap
#' @import patchwork
#' @importFrom utils hasName
#' @export
#'
#' @details
#' Binary segmentation is a classic algorithm for change-point detection that recursively
#' splits the data at locations that minimise the cost function.
#'
#' `binSeg` requires  a `R6` object of class `costFunc`, which can be created via `costFunc$new()`. Currently, the following cost functions are supported:
#'
#' - `"L1"` and `"L2"` for (independent) piecewise Gaussian process with **constant variance**
#' - `"SIGMA"`: for (independent) piecewise Gaussian process with **varying variance**
#' - `"VAR"`: for piecewise Gaussian vector-regressive process with **constant noise variance**
#' - `"LinearL2"`: for piecewise linear regression process with **constant noise variance**
#'
#' See `$eval()` method for more details on computation of cost.
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
#' # Default L2 cost function
#' binSegObj = binSeg$new(minSize = 1L, jump = 1L)
#' binSegObj$fit(signals)
#' binSegObj$predict(pen = 100)
#' binSegObj$plot()
#'
#' ## SIGMA example
#' set.seed(111)
#' signals = as.matrix(c(rnorm(100,-5,1),
#'                       rnorm(100,-5,10),
#'                       rnorm(100,-5,1)))
#' # L2 cost function
#' binSegObj = binSeg$new(minSize = 1L, jump = 1L)
#' binSegObj$fit(signals)
#' # We choose pen = 50.
#' binSegObj$predict(pen = 50)
#' binSegObj$plot()
#'
#' # The standard L2 cost function is not suitable.
#' # Use the SIGMA cost function.
#' binSegObj$costFunc = costFunc$new(costFunc = "SIGMA")
#' binSegObj$predict(pen = 50)
#' binSegObj$plot()
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `binSeg` object.}
#'   \item{\code{$describe()}}{Describes the `binSeg` object.}
#'   \item{\code{$fit()}}{Constructs a `binSeg` module in `C++`.}
#'   \item{\code{$eval()}}{Evaluates the cost of a segment.}
#'   \item{\code{$predict()}}{Performs `binSeg` given a linear penalty value.}
#'   \item{\code{$plot()}}{Plots change-point segmentation in `ggplot` style.}
#'   \item{\code{$clone()}}{Clones the `R6` object.}
#' }
#'
#' @references
#' Truong, C., Oudre, L., & Vayatis, N. (2020). Selective review of offline change point detection methods.
#' Signal Processing, 167, 107299.
#'
#' Hocking, T. D. (2024). Finite Sample Complexity Analysis of Binary Segmentation. arXiv preprint arXiv:2410.08654.
#'
#' @author
#' Minh Long Nguyen \email{edelweiss611428@gmail.com} \cr
#' Toby Dylan Hocking \email{toby.hocking@r-project.org} \cr
#' Charles Truong \email{ctruong@ens-paris-saclay.fr}
#'
#' @export

binSeg = R6Class(
  "binSeg",

  private = list(

    .minSize = 1L,
    .jump = 1L,
    .costFunc = costFunc$new("L2"), #L2 cost function
    .binSegModule = NULL,
    .tsMat = NULL,
    .covariates = NULL,
    .fitted = FALSE,
    .n = NULL,
    .p = NULL,
    .tmpEndPts = NULL, #Temporary end points
    .tmpPen = NULL #Temporary penalty value

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

      # If time series data exists, refit the model
      if (!is.null(private$.tsMat) & private$.fitted) {
        message("`costFunc` has been updated. Re-fitting the model.")
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

      # If time series data exists, refit the model
      if (!is.null(private$.tsMat) & private$.fitted) {
        message("`costFunc` has been updated. Re-fitting the model.")
        self$fit()
      }

    },

    #' @field costFunc `R6` object of class `costFunc`. Search grid step size. Can be accessed or modified via `$costFunc`.
    #' Modifying `costFunc` will automatically trigger `$fit()`.
    costFunc = function(costFuncObj) {

      if(missing(costFuncObj)){
        return(private$.costFunc)
      }

      if (!inherits(costFuncObj, "costFunc") | !is.R6(costFuncObj)) {
        stop("`costFunc` must be a `R6` object of class `costFunc` - can be created via costFunc$new()!")
      }

      private$.costFunc = costFuncObj

      # If time series data exists, refit the model

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

      if(private$.costFunc$pass()[["costFunc"]] %in% c("LinearL2")){
        if (!is.null(private$.tsMat) & private$.fitted) {
          self$fit()
        }
      }

    }


  ),
  public = list(

    #' @description Initialises a `binSeg` object.
    #'
    #' @param minSize Integer. Minimum allowed segment length. Default: `1L`.
    #' @param jump Integer. Search grid step size: only positions in \{k, 2k, ...\} are considered. Default: `1L`.
    #' @param costFunc A `R6` object of class `costFunc`. Should be created via `costFunc$new()` to avoid error.
    #' Default: `costFunc$new("L2")`.
    #'
    #' @return Invisibly returns `NULL`.

    initialize = function(minSize, jump, costFunc) {

      if(!missing(minSize)){
        self$minSize = minSize
      }

      if(!missing(jump)){
        self$jump = jump
      }

      if(!missing(costFunc)){
        self$costFunc = costFunc
      }

      invisible(NULL)
    },

    #' @description Describes a `binSeg` object.
    #'
    #' @param printConfig Logical. Whether to print object configurations. Default: `FALSE`.
    #'
    #' @return Invisibly returns a list storing at least the following fields:
    #'
    #' \describe{
    #'   \item{\code{minSize}}{Minimum allowed segment length.}
    #'   \item{\code{jump}}{Search grid step size.}
    #'   \item{\code{costFunc}}{The `costFun` object.}
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
                    costFunc = private$.costFunc,
                    fitted = private$.fitted,
                    tsMat = private$.tsMat,
                    covariates = private$.covariates,
                    n = private$.n,
                    p = private$.p)


      if(printConfig){

        cat(sprintf("Binary Segmentation (binSeg) \n"))
        cat(sprintf("minSize      : %sL\n", private$.minSize))
        cat(sprintf("jump         : %sL\n", private$.jump))
        cat(sprintf("costFunc     : \"%s\"\n", private$.costFunc$pass()[["costFunc"]]))

      }


      if(private$.costFunc$pass()[["costFunc"]] == "SIGMA"){

        if(printConfig){

          cat(sprintf("addSmallDiag : %s\n", private$.costFunc$pass()[["addSmallDiag"]]))
          cat(sprintf("epsilon      : %s\n", private$.costFunc$pass()[["epsilon"]]))

        }

        params[["addSmallDiag"]] = private$.costFunc$pass()[["addSmallDiag"]]
        params[["epsilon"]] = private$.costFunc$pass()[["epsilon"]]

      }

      if(private$.costFunc$pass()[["costFunc"]] == "VAR"){

        if(printConfig){

          cat(sprintf("pVAR         : %sL\n", private$.costFunc$pass()[["pVAR"]]))

        }

        params[["pVAR"]] = private$.costFunc$pass()[["pVAR"]]

      }

      if(private$.costFunc$pass()[["costFunc"]] == "LinearL2"){

        if(printConfig){

          cat(sprintf("intercept    : %sL\n", private$.costFunc$pass()[["intercept"]]))

        }

        params[["intercept"]] = private$.costFunc$pass()[["intercept"]]

      }

      if(printConfig){

        cat(sprintf("fitted       : %s\n", private$.fitted))
        cat(sprintf("n            : %sL\n", private$.n))
        cat(sprintf("p            : %sL\n", private$.p))

      }

      invisible(params)

    },

    #' @description Constructs a `C++` module for binary segmentation.
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
    #' @details This method constructs a `C++` `binSeg` module and sets `private$.fitted` to `TRUE`,
    #' enabling the use of `$predict()` and `$eval()`. Some precomputations are performed to allow
    #' `$predict()` to run in linear time with respect to the number of data points
    #' (see `$predict()` for more details).

    fit = function(tsMat = NULL, covariates = NULL) {

      # Only assign if explicitly called with `tsMat` argument

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
        } #else not null because `private$.tsMat` is set via active binding

      }

      if(private$.costFunc$pass()[["costFunc"]] %in% c("LinearL2")){

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

        } else{ #if covariates is null

          if(is.null(private$.covariates)){
            warning("No `covariates` found! Force-fitting with only an intercept!")

            if(private$.costFunc$pass()[["costFunc"]] == "LinearL2"){

              private$.binSegModule = new(binSegCpp_LinearL2, private$.tsMat, matrix(1, nrow = private$.n, ncol = 1),
                                          FALSE, #no intercept
                                          private$.minSize, private$.jump)

            }

            private$.fitted = TRUE
            private$.binSegModule$fit()

            return(invisible(NULL))


          } else{

            if(nrow(private$.covariates) != private$.n){
              stop("Numbers of observations in `covariates` and `tsMat` do not match!")
            }

          }
        }
      }


      if(private$.costFunc$pass()[["costFunc"]] == "L2"){
        private$.binSegModule = new(binSegCpp_L2, private$.tsMat, private$.minSize, private$.jump)

      } else if(private$.costFunc$pass()[["costFunc"]] == "VAR"){
        private$.binSegModule = new(binSegCpp_VAR, private$.tsMat, private$.costFunc$pass()[["pVAR"]],
                                    private$.minSize, private$.jump)

      } else if(private$.costFunc$pass()[["costFunc"]] == "SIGMA"){
        private$.binSegModule = new(binSegCpp_SIGMA, private$.tsMat,
                                    private$.costFunc$pass()[["addSmallDiag"]],
                                    private$.costFunc$pass()[["epsilon"]],
                                    private$.minSize, private$.jump)

      } else if(private$.costFunc$pass()[["costFunc"]] == "L1"){

        private$.binSegModule = new(binSegCpp_L1_cwMed, private$.tsMat,
                                    private$.minSize, private$.jump)

      } else if(private$.costFunc$pass()[["costFunc"]] == "LinearL2"){

        private$.binSegModule = new(binSegCpp_LinearL2, private$.tsMat, private$.covariates,
                                    private$.costFunc$pass()[["intercept"]],
                                    private$.minSize, private$.jump)

      } else{
        # nocov start
        stop("Cost function not supported!")
        # nocov end
      }

      private$.fitted = TRUE

      private$.binSegModule$fit()

      invisible(NULL)
    },

    #' @description Evaluate the cost of the segment (a,b]
    #'
    #' @param a Integer. Start index of the segment (exclusive). Must satisfy \code{start < end}.
    #' @param b Integer. End index of the segment (inclusive).
    #'
    #' @return The segment cost.
    #'
    #' @details
    #' The segment cost is evaluated as follows:
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

      return(private$.binSegModule$eval(a, b))

    },

    #' @description Performs `binSeg` given a linear penalty value.
    #'
    #' @param pen Numeric. Penalty per change-point. Default: `0`.
    #'
    #' @return An integer vector of regime end-points. By design, the last element is the
    #' number of observations.
    #'
    #' @details The algorithm recursively partitions a time series to detect multiple change-points.
    #' At each step, the algorithm identifies the segment that, if split, would result in the greatest reduction in total cost.
    #' This process continues until no further splits are possible (e.g., each segment is of minimal length or each breakpoint
    #' corresponds to a single data point).
    #'
    #' Then, the algorithm selects the "optimal" set of break-points given the linear penalty threshold. Let \eqn{[c_1, \dots, c_k, c_{k+1}]} denote the set of segment end-points with \eqn{c_1 < c_2 < \dots < c_k < c_{k+1} = n},
    #' where \eqn{k} is the number of detected change-points and \eqn{n} is the total number of data points.
    #' and \eqn{k} is the number of change-points. Let \eqn{c_{(c_i, c_{i+1}]}} be the cost of segment \eqn{(c_i, c_{i+1}]}.
    #' The total penalised cost is then
    #' \deqn{
    #'   \text{TotalCost} = \sum_{i=1}^{k+1} c_{(c_i, c_{i+1}]} + \lambda \cdot k,
    #' }
    #' where \eqn{\lambda} is a linear penalty applied per change-point. We then optimise over
    #' \eqn{k} to minimise the penalised cost function.
    #'
    #' This approach allows detecting multiple change-points in a time series while controlling
    #' model complexity through the linear penalty threshold.
    #'
    #' In our implementation, the recursive step is carried out during `$fit()`.
    #' Therefore, `$predict()` runs in linear time with respect to the number of data points.
    #'
    #' Temporary segment end-points are saved to `private$.tmpEndPoints` after `$predict()`, enabling users to call `$plot()` without
    #' specifying endpoints manually.

    predict = function(pen = 0){

      if(!private$.fitted){
        stop("`$fit()` must be run before `$predict()`!")
      }

      if(is.null(pen)){
        stop("`pen` must not be null!")
      }

      if(!is.numeric(pen) | length(pen)!= 1 |  any(pen < 0)){
        stop("`pen` must be a single non-negative value!")
      }

      endPts = sort(c(private$.binSegModule$predict(pen), private$.n))

      private$.tmpEndPts = endPts
      private$.tmpPen = pen

      return(endPts)

    },


    #' @description Plots change-point segmentation
    #'
    #' @param d Integer vector. Dimensions to plot. Default: `1L`.
    #' @param endPts Integer vector. End points. Default: latest temporary changepoints obtained via `$predict()`.
    #' @param dimNames Character vector. Feature names matching length of `d`. Defaults to `"X1", "X2", ...`.
    #' @param main Character. Main title. Defaults to `"binSeg: d = ..."`.
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
        main = paste0("binSeg: ", paste0("d = (", toString(d), ")"))

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

      # Build long-format dataframe for all selected dimensions
      tsList = lapply(seq_along(d), function(i) {
        data.frame(
          time = 1:private$.n,
          value = private$.tsMat[, d[i]],
          dimension = dimNames[i]
        )
      })

      allTsDf = do.call(rbind, tsList)

      # Create segment info
      segInt = data.frame(
        xmin = c(1, endPts[-length(endPts)]),
        xmax = endPts,
        fill = rep(bgCol, length.out = length(endPts))
      )

      # Build plot with facet_wrap
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

