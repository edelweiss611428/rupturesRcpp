#' Binary Segmentation (binSeg)
#'
#' @description An R6 class implementing binary segmentation for offline change point detection.
#'
#' @include createCostFunc.R
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom ggplot2 aes ggplot geom_rect geom_line scale_fill_identity theme_minimal theme geom_vline labs element_blank element_text facet_wrap
#' @import patchwork
#' @export
#'
#' @details
#' Binary segmentation is a classic algorithm for change point detection that recursively
#' splits the data at locations that minimise the cost function.
#'
#' Currently supports the following cost functions:
#'
#' - `"L2"`: for (independent) piecewise Gaussian process with **constant variance**
#' - `"SIGMA"`: for (independent) piecewise Gaussian process with **varying variance**
#' - `"VAR"`: for piecewise Gaussian vector-regressive process with **constant variance**
#'
#' `binSeg` requires  a `costFunc` object, which can be created via `createCostFunc()`.
#'
#'  Basic usage: <https://github.com/edelweiss611428/rupturesRcpp/tree/main/README.md>
#'
#'  See `PELT$eval()` method for more details on computation of cost.
#'
#' @examples
#'
#' # 2-regime simulated data example
#' set.seed(1)
#' tsMat = cbind(c(rnorm(100,0), rnorm(100,5,5)),
#'               c(rnorm(100,0), rnorm(100,5,5)))
#' # Create a `"SIGMA` cost object.
#' SIGMAObj = createCostFunc(costFunc = "SIGMA")
#' # Initialise a `binSeg` object
#' binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFuncObj = SIGMAObj)
#' # Input the time series matrix and perform binary segmentation for the maximum number of regimes
#' binSegObj$fit(tsMat)
#' # Describe the `binSeg` object
#' binSegObj$describe()
#' # Perform binSeg with `pen = 100`
#' binSegObj$predict(pen = 100)
#' # Plot the segmentation results
#' binSegObj$plot(d = 1:2, main = "method: binSeg; costFunc: SIGMA; pen: 100")
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `binSeg` object.}
#'   \item{\code{$describe()}}{Describes the `binSeg` object.}
#'   \item{\code{$fit()}}{Takes a time series matrix as input and perform `binSeg` for the
#' maximum number of change points.}
#'   \item{\code{$predict()}}{Performs binSeg given a linear penalty value.}
#'   \item{\code{$plot()}}{Plots change point segmentation in ggplot style.}
#'   \item{\code{$clone()}}{Clones the `binSeg` object.}
#' }
#'
#' @references
#' Truong, C., Oudre, L., & Vayatis, N. (2020). Selective review of offline change point detection methods.
#' Signal Processing, 167, 107299.
#'
#' Hocking, T. D. (2024). Finite Sample Complexity Analysis of Binary Segmentation. arXiv preprint arXiv:2410.08654.
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @export

binSeg = R6Class(
  "binSeg",

  private = list(

    .minSize = 1L,
    .jump = 1L,
    .costFuncObj = createCostFunc(), #L2 cost function
    .costModule = NULL,
    .tsMat = NULL,
    .fitted = FALSE,
    .n = NULL,
    .p = NULL,
    .bkps = NULL,
    .cost = NULL,
    .tmpEndPts = NULL, #Temporary end points
    .tmpPen = NULL #Temporary penalty value

  ),

  active = list(

    #' @field minSize Active binding. Sets the internal variable \code{.minSize} but should not be called directly.
    minSize = function(intVal) {
      if (!is.numeric(intVal) || any(intVal < 1) || length(intVal) != 1) {
        stop("minSize must be a single positive integer!")
      }
      private$.minSize = as.integer(intVal)
    },

    #' @field jump Active binding. Sets the internal variable \code{.jump} but should not be called directly.
    jump = function(intVal) {
      if (!is.numeric(intVal) || any(intVal < 1) || length(intVal) != 1) {
        stop("jump must be a single positive integer!")
      }
      private$.jump = as.integer(intVal)
    },

    #' @field costFuncObj Active binding. Sets the internal variable \code{.costFuncObj} but should not be called directly.
    costFuncObj = function(Obj) {

      if (!inherits(Obj, "costFunc") | !is.list(Obj)) {
        stop("costFuncObj must be a costFunc object! See createCostObj()!")
      }

      private$.costFuncObj = Obj
    },

    #' @field tsMat Active binding. Sets the internal variable \code{.tsMat} but should not be called directly.
    tsMat = function(numMat) {
      if (!is.numeric(numMat) || !is.matrix(numMat)) {
        stop("tsMat must be a numeric time series matrix!")
      }
      private$.tsMat = numMat
    }
  ),

  public = list(

    #' @description Initialises a `binSeg` object.
    #'
    #' @param minSize Integer. Minimum allowed segment length. Default: `1L`.
    #' @param jump Integer. Search grid step size: only positions in \{1, k+1, 2k+1, ...\} are considered. Default: `1L`.
    #' @param costFuncObj List of class `costFunc`. Should be created via `createFuncObj()` to avoid unintended error.
    #' Default, `createFuncObj("L2")`.
    #' @return Invisibly returns `NULL`.

    initialize = function(minSize, jump, costFuncObj) {

      if(!missing(minSize)){
        self$minSize = minSize
      }

      if(!missing(jump)){
        self$jump = jump
      }

      if(!missing(costFuncObj)){
        self$costFuncObj = costFuncObj
      }

      print("You have created a binSeg object!")

      invisible(NULL)
    },

    #' @description Describes a `binSeg` object.
    #'
    #' @return Invisibly returns a list containing some of the following fields:
    #' \describe{
    #'   \item{\code{minSize}}{Minimum allowed segment length.}
    #'   \item{\code{jump}}{Search grid step size.}
    #'   \item{\code{costFuncObj}}{The `costFun` object.}
    #'   \item{\code{fitted}}{Whether or not `$fit()` has been run.}
    #'   \item{\code{tsMat}}{Input time series matrix.}
    #'   \item{\code{n}}{Number of observations.}
    #'   \item{\code{p}}{Number of features.}
    #'   \item{\code{bkps}}{Vector of unordered breakpoint positions.}
    #'   \item{\code{cost}}{From the second element, split costs corresponding to \code{bkps}; first element is total cost without splits.}
    #' }

    describe = function() {

      params = list(minSize = private$.minSize,
                    jump = private$.minSize,
                    costFuncObj = private$.costFuncObj,
                    fitted = private$.fitted,
                    tsMat = private$.tsMat,
                    n = private$.n,
                    p = private$.p,
                    bkps = private$.bkps,
                    cost = private$.cost)

      cat(sprintf("Binary Segmentation (binSeg) \n"))
      cat(sprintf("minSize      : %sL\n", private$.minSize))
      cat(sprintf("jump         : %sL\n", private$.jump))
      cat(sprintf("costFunc     : \"%s\"\n", private$.costFuncObj$costFunc))

      if(private$.costFuncObj$costFunc == "SIGMA"){
        cat(sprintf("addSmallDiag : %s\n", private$.costFuncObj$addSmallDiag))
        cat(sprintf("epsilon      : %s\n", private$.costFuncObj$epsilon))
        params[["addSmallDiag"]] = private$.costFuncObj$addSmallDiag
        params[["epsilon"]] = private$.costFuncObj$epsilon
      }

      if(private$.costFuncObj$costFunc == "VAR"){
        cat(sprintf("pVAR         : %s\n", private$.costFuncObj$pVAR))
        params[["pVAR"]] = private$.costFuncObj$pVAR
      }

      cat(sprintf("fitted       : %s\n", private$.fitted))
      cat(sprintf("n            : %sL\n", private$.n))
      cat(sprintf("p            : %sL\n", private$.p))

      invisible(params)

    },

    #' @description Takes a time series matrix as input and performs binSeg for the
    #' maximum number of change points.
    #'
    #' @param tsMat Numeric matrix. A time series matrix of size \eqn{n \times p} whose rows are observations ordered in time.
    #'
    #' @return Invisibly returns `NULL`.
    #'
    #' @details
    #' This method does the following:
    #'- Initialises `private$.tsMat`, `private$.n`, and `private$.p`.
    #'- Sets `private$.fitted` to `TRUE`, enabling the use of `$predict()` and `$eval()`.
    #'- Initialises a cost module corresponding to `tsMat` and the `costFunc`
    #'  object, enabling the use of `$eval()`.
    #'- Performs binSeg for the maximum number of change points, making `$predict()`
    #'  more efficient.


    fit = function(tsMat) {

      self$tsMat = tsMat
      private$.n = nrow(tsMat)
      private$.p = ncol(tsMat)
      private$.fitted = TRUE
      detection = binSegCpp(private$.tsMat, private$.minSize, private$.jump,
                            costFuncObj = private$.costFuncObj)
      private$.cost = detection$cost
      private$.bkps = detection$bkps

      if(private$.costFuncObj$costFunc == "L2"){
        private$.costModule = new(rupturesRcpp::Cost_L2, tsMat)

      } else if(private$.costFuncObj$costFunc == "SIGMA"){
        private$.costModule = new(rupturesRcpp::Cost_SIGMA, tsMat,
                                  private$.costFuncObj$addSmallDiag, private$.costFuncObj$epsilon)

      } else if(private$.costFuncObj$costFunc == "VAR"){
        private$.costModule = new(rupturesRcpp::Cost_VAR, tsMat,
                                  private$.costFuncObj$pVAR)

      } else {
        stop("Cost function not supported!")

      }

      invisible(NULL)
    },

    #' @description Evaluate the cost of the segment (a,b]
    #'
    #' @param a Integer. Start index of the segment (exclusive). Must satisfy \code{start < end}.
    #' @param b Integer. End index of the segment (inclusive).
    #'
    #' @return The segment cost
    #'
    #' @details
    #' The segment cost is evaluated as follows:
    #'
    #' - **L2 cost function**:
    #' \deqn{c_{L_2}(y_{(a+1)..b}) := \sum_{t = a+1}^{b} \| y_t - \bar{y}_{(a+1)..b} \|_2^2}
    #' where \eqn{\bar{y}_{(a+1)..b}} is the empirical mean of the segment. If
    #' `a+1 = b`, return 0.
    #'
    #' - **SIGMA cost function**:
    #' \deqn{c_{\sum}(y_{(a+1)..b}) := (b - a)\log \det \hat{\Sigma}_{(a+1)..b}} where \eqn{\hat{\Sigma}_{(a+1)..b}} is
    #' the empirical covariance matrix of the segment without Bessel correction. Here, if `addSmallDiag = TRUE`, a small
    #' bias `epsilon` is added to the diagonal of estimated covariance matrices to improve numerical stability. If
    #' \eqn{\hat{\Sigma}} is singular, return 0. If `a+1 = b`, return 0.
    #'
    #' - **VAR(r) cost function**:
    #' \deqn{c_{\mathrm{VAR}}(y_{(a+1)..b}) := \sum_{t = a+r+1}^{b} \left\| y_t - \sum_{j=1}^r \hat A_j y_{t-j} \right\|_2^2}
    #' where \eqn{\hat A_j} are the estimated VAR coefficients, commonly estimated via the OLS criterion. An approximate linear
    #' solver will be used when exact `arma::solve()` fails. If no solution found, return 0. If `a-b < p*r+1` (i.e., not enough observations),
    #' return 0.
    eval = function(a, b){

      if(!private$.fitted){
        stop("$fit() must be run before $eval()!")
      }

      if (!is.numeric(start) || any(start < 0) || length(start) != 1 || any(start > private$.n)) {
        stop("`0 <= start <= n` must be true!")
      }

      if (!is.numeric(end) || any(end < 0) || length(end) != 1 || any(end>private$.n)) {
        stop("`0 <= end <= n` must be true!")
      }

      start = as.integer(start)
      end = as.integer(end)
      if(start >= end){
        stop("start must be smaller than end!")
      }

      return(private$.costModule$eval(start, end))

    },

    #' @description Performs binSeg given a linear penalty value.
    #'
    #' @param pen Numeric. Penalty per change point. Default: `0`.
    #'
    #' @return An integer vector of regime end points. By design, the last element is the
    #' number of observations.
    #'
    #' @details Performs binSeg given a linear penalty value. Temporary end points are saved
    #' to `private$.tmpEndPoints`, allowing users to use `$plot()` without specifying
    #' end points.

    predict = function(pen = 0){

      if(!private$.fitted){
        stop("$fit() must be run before $predict()!")
      }

      if(!is.numeric(pen) || length(pen)!= 1 ||  any(pen < 0)){
        stop("pen must be a single non-negative numeric value!")
      }

      endPts = c(sort(binSegPredCpp(private$.bkps, private$.cost, pen)), private$.n)

      private$.tmpEndPts = endPts
      private$.tmpPen = pen

      return(endPts)
    },


    #' @description Plots change point segmentation
    #'
    #' @param d Integer vector. Dimensions to plot. Default: `1L`.
    #' @param endPts Integer vector. End points; defaults to latest temporary changepoints from `$predict()`.
    #' @param dimNames Character vector. Feature names matching length of `d`. Defaults to `"X1", "X2", ...`.
    #' @param main Character. Main title. Defaults to `"binSeg: d = ..."`.
    #' @param xlab Character. X-axis label. Default: `"Time"`.
    #' @param tsWidth Numeric. Line width for time series and segments. Default: `0.25`.
    #' @param tsCol Character. Time series color. Default: `"#5B9BD5"`.
    #' @param bgCol Character vector. Segment colors, recycled to length of `endPts`. Default: `c("#A3C4F3", "#FBB1BD")`.
    #' @param bgAlpha Numeric. Background transparency. Default: `0.5`.
    #' @param ncol Integer. Number of columns in facet layout. Default: `1L`.
    #'
    #' @details Plots change point segmentation results. Based on `ggplot2`. Multiple plots can easily be
    #' horizontally and vertically stacked using `patchwork`'s operators `/` and `|`, respectively.
    #'
    #' @return An object of classes `gg` and `ggplot`.

    plot = function(d = 1L, endPts, dimNames, main, xlab, tsWidth = 0.25,
                    tsCol = "#5B9BD5",
                    bgCol = c("#A3C4F3", "#FBB1BD"),
                    bgAlpha = 0.5,
                    ncol = 1L){


      if(missing(main)){
        main = paste0("binSeg: ", title_text = paste0("d = (", toString(d), ")"))

      } else {
        if(!is.character(main) || length(main )!= 1L){
          stop("main must be a single character!")
        }
      }

      if(missing(xlab)){
        xlab = "Time"

      } else {
        if(!is.character(xlab) || length(xlab)!= 1L){
          stop("xlab must be a single character!")
        }
      }

      if(missing(endPts)){
        message("endPts is missing. Proceed to use the temporary endPts!")

        if(is.null(private$.tmpEndPts)){
          stop("Temporary endPts is NULL. Must run $predict() to obtain this!")
        }
      } else{
        if(!is.numeric(endPts)){
          stop("endPts must be a numeric/integer vector specifying endpoints! Could be obtained via $predict().")
        } else{
          endPts = as.integer(sort(endPts))

          if(min(endPts) < 1){
            stop("min(endPts) >= 1!")
          }

          if(max(endPts) != private$.n){
            stop("By construction, max(endPts) must be n! Can use $predict() to obtain endPts.")
          }

          if(length(unique(endPts)) != length(endPts)){
            stop("as.integer(endPts) contains duplicated elements!")
          }
        }
      }

      if (!is.numeric(d)) {
        stop("d must be a numeric/integer vector specifying dimensions! e.g., 1, or c(1,2).")
      }

      if(missing(dimNames)){
        message("dimNames is missing. Proceed to use the default dimNames! e.g., paste0('X', d)).")
        dimNames = paste0("X", d)
      } else {
        if (!is.character(dimNames)) {
          stop("dimNames must be a character vector specifying feature names!")
          if(length(dimNames) != length(d)){
            stop("length(dimNames) != length(d)!")
          }
        }
      }

      d = as.integer(d)
      if(any(d < 1) | any(d > private$.p)){
        stop("Dimensions must be between [1,p]!")
      }

      endPts = private$.tmpEndPts

      # Build long-format dataframe for all selected dimensions
      tsList <- lapply(seq_along(d), function(i) {
        data.frame(
          time = 1:private$.n,
          value = private$.tsMat[, d[i]],
          dimension = dimNames[i]
        )
      })

      allTsDf <- do.call(rbind, tsList)

      # Create segment info
      segInt <- data.frame(
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

