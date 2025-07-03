#' Pruned Exact Linear Time (PELT)
#'
#' @description An R6 class implementing the PELT algorithm for offline change point detection.
#'
#' @include createCostFunc.R
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom ggplot2 aes ggplot geom_rect geom_line scale_fill_identity theme_minimal theme geom_vline labs element_blank element_text facet_wrap
#' @import patchwork
#'
#' @export
#'
#' @details
#' PELT (Pruned Exact Linear Time) is an efficient algorithm for change point detection
#' that prunes the search space to achieve optimal segmentation in linear time under certain conditions.
#'
#' Currently support the following cost functions:
#'
#' - `"L2"`: for (independent) piecewise Gaussian process with **constant covariance**
#' - `"SIGMA"`: for (independent) piecewise Gaussian process with **varying covariance**
#' - `"VAR"`: for piecewise Gaussian vector-regressive process with **constant covariance**
#'
#' Cost function objects can be created via `createCostFunc()`.
#'
#' See **Methods** section for more details.
#'
#' @examples
#' # Toy example
#' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,0, 10)))
#' # Initialise a `PELT` object and fit `PELT` to `tsMat`
#' PELTObj = PELT$new(costFuncObj = createCostFunc("SIGMA"))
#' PELTObj$fit(tsMat)
#' # Perform PELT for a specific linear penalty threshold
#' PELTObj$predict(pen = 50)
#' # Plot the latest segmentation solution
#' PELTObj$plot(main = "PELT:SIGMA:pen=50", ncol = 1)
#' # Describe the `PELT` object (and invisibly return the object's fields)
#' PELTObj$describe()
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `PELT` object.}
#'   \item{\code{$describe()}}{Describes the `PELT` object.}
#'   \item{\code{$fit()}}{Takes a time series matrix as input.}
#'   \item{\code{$predict()}}{Performs PELT given a linear penalty value.}
#'   \item{\code{$plot()}}{Plots change point segmentation in ggplot style.}
#'   \item{\code{$clone()}}{Clone the `PELT` object.}
#' }
#'
#' @references
#' Truong, C., Oudre, L., & Vayatis, N. (2020). Selective review of offline change point detection methods.
#' Signal Processing, 167, 107299.
#'
#' Killick, R., Fearnhead, P., & Eckley, I. A. (2012). Optimal detection of change points with a linear computational cost.
#' Journal of the American Statistical Association, 107(500), 1590-1598.
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @docType class
#'
#' @importFrom R6 R6Class
#' @export

PELT = R6Class(
  "PELT",

  private = list(

    .minSize = 1L,
    .jump = 1L,
    .tsMat = NULL,
    .fitted = FALSE,
    .n = NULL,
    .p = NULL,
    .costFuncObj = createCostFunc(), #L2 cost function
    .tmpEndPts = NULL, #Temporary end points - obtained after running $predict()
    .tmpPen = NULL #Temporary penalty value - obtained after running $predict()
  ),

  active = list(

    #' @field minSize Active binding. Sets the internal variable \code{.minSize} but should not be called directly.
    minSize = function(intVal) {
      if (is.null(intVal) || !is.numeric(intVal) || (as.integer(intVal) < 1) || length(intVal) != 1) {
        stop("minSize must be a single positive integer!")
      }
      private$.minSize = as.integer(intVal)
    },

    #' @field jump Active binding. Sets the internal variable \code{.jump} but should not be called directly.
    jump = function(intVal) {
      if (is.null(intVal) || !is.numeric(intVal) || (as.integer(intVal) < 1) || length(intVal) != 1) {
        stop("jump must be a single positive integer!")
      }
      private$.jump = as.integer(intVal)
    },

    #' @field costFuncObj Active binding. Sets the internal variable \code{.costFuncObj} but should not be called directly.
    costFuncObj = function(Obj) {
      if (!inherits(Obj, "costFunc") | !is.list(Obj)) {
        stop("costFuncObj must be a costFunc object! See createCostObj()")
      }
      private$.costFuncObj = Obj
    },

    #' @field tsMat Active binding. Sets the internal variable \code{.tsMat} but should not be called directly.
    tsMat = function(numMat) {
      if (is.null(numMat) || !is.numeric(numMat) || !is.matrix(numMat)) {
        stop("tsMat must be a numeric time series matrix!")
      }
      private$.tsMat = numMat
    }
  ),

  public = list(

    #' @description Initialises a `PELT` object.
    #'
    #' @param minSize Integer. Minimum allowed segment length. Default: `1L`.
    #' @param jump Integer. Search grid step size: only positions in \{1, k+1, 2k+1, ...\} are considered. Default: `1L`.
    #' @param costFuncObj List of class `costFunc`. Created via `costFuncObj()` function. Default, `costFuncObj("L2")`.
    #' @return Invisibly returns `NULL`.
    #'
    #' @examples
    #' L2Obj = createCostFunc("L2")
    #' PELTObj = PELT$new(minSize = 1L, jump = 1L, costFuncObj = L2Obj)

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

      print("You have created a PELT object!")

      invisible(NULL)
    },

    #' @description Describes a `PELT` object.
    #'
    #' @return Invisibly returns a list containing some of the following fields:
    #' \describe{
    #'   \item{\code{minSize}}{Minimum allowed segment length.}
    #'   \item{\code{jump}}{Search grid step size.}
    #'   \item{\code{costFuncObj}}{The `costFunc` object.}
    #'   \item{\code{fitted}}{Whether or not `$fit()` has been run.}
    #'   \item{\code{tsMat}}{Input time series matrix.}
    #'   \item{\code{n}}{Number of observations.}
    #'   \item{\code{p}}{Number of features.}
    #' }
    #' @examples
    #' L2Obj = createCostFunc("L2")
    #' PELTObj = PELT$new(minSize = 1L, jump = 1L, costFuncObj = L2Obj)
    #' PELTObj$describe()

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

      cat(sprintf("Pruned Exact Linear Time (PELT) \n"))
      cat(sprintf("minSize      : %sL\n", private$.minSize))
      cat(sprintf("jump         : %sL\n", private$.jump))
      cat(sprintf("costFunc.    : \"%s\"\n", private$.costFuncObj$costFunc))

      if(private$.costFuncObj$costFunc == "SIGMA"){
        cat(sprintf("addSmallDiag : %s\n", private$.costFuncObj$addSmallDiag))
        cat(sprintf("epsilon      : %s\n", private$.costFuncObj$epsilon))
        params[["addSmallDiag"]] = private$.costFuncObj$addSmallDiag
        params[["epsilon"]] = private$.costFuncObj$epsilon
      }

      if(private$.costFuncObj$costFunc == "VAR"){
        cat(sprintf("pVAR.        : %s\n", private$.costFuncObj$pVAR))
        params[["pVAR"]] = private$.costFuncObj$pVAR
      }

      cat(sprintf("fitted       : %s\n", private$.fitted))
      cat(sprintf("n            : %sL\n", private$.n))
      cat(sprintf("p            : %sL\n", private$.p))

      invisible(params)

    },

    #' @description Takes a time series matrix as input.
    #'
    #' @param tsMat Numeric matrix. A time series matrix of size \eqn{n \times p} whose rows are observations ordered in time.
    #'
    #' @return Invisibly returns `NULL`.
    #'
    #' @details Initialises `private$.tsMat`, `private$.n`, and `private$.p`, and sets private$.fitted to TRUE,
    #' enabling the use of `$predict()`. Run `$describe()` for detailed configurations.
    #'
    #' @examples
    #' L2Obj = createCostFunc("L2")
    #' PELTObj = PELT$new(minSize = 1L, jump = 1L, costFuncObj = L2Obj)
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' PELTObj$fit(tsMat)

    fit = function(tsMat) {
      self$tsMat = tsMat
      private$.n = nrow(tsMat)
      private$.p = ncol(tsMat)
      private$.fitted = TRUE
      invisible(NULL)
    },

    #' @description Performs PELT given a linear penalty value.
    #' @param pen Numeric. Penalty per change point. Default: `0`.
    #'
    #' @return An integer vector of regime end points. By design, the last element is the
    #' number of observations.
    #'
    #' @details Performs PELT given a linear penalty value. Temporary end points are saved
    #' to `private$.tmpEndPoints`, allowing users to use `$plot()` without specifying
    #' end points.
    #'
    #' @examples
    #' L2Obj = createCostFunc("L2")
    #' PELTObj = PELT$new(minSize = 1L, jump = 1L, costFuncObj = L2Obj)
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' PELTObj$fit(tsMat)
    #' PELTObj$predict()

    predict = function(pen = 0){

      if(!private$.fitted){
        stop("$fit() must be run before $predict()!")
      }

      if(!is.numeric(pen) || length(pen)!= 1 || any(pen < 0)){
        stop("pen must be a single non-negative numeric value!")
      }

      endPts = PELTCpp(private$.tsMat, pen, private$.minSize, private$.jump,
                         costFuncObj = private$.costFuncObj)

      private$.tmpEndPts = endPts
      private$.tmpPen = pen

      return(endPts)
    },

    #' @description Plots change point segmentation
    #'
    #' @param d Integer vector. Dimensions to plot. Default: `1L`.
    #' @param endPts Integer vector. End points; defaults to latest temporary changepoints from `$predict()`.
    #' @param dimNames Character vector. Feature names matching length of `d`. Defaults to `"X1", "X2", ...`.
    #' @param main Character. Main title. Defaults to `"PELT: d = ..."`.
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
    #'
    #' @examples
    #' L2Obj = createCostFunc("L2")
    #' PELTObj = PELT$new(minSize = 1L, jump = 1L, costFuncObj = L2Obj)
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' PELTObj$fit(tsMat)
    #' PELTObj$predict(pen = 1)
    #' pen1 = PELTObj$plot(main = "PELT: pen = 1")
    #' PELTObj$predict(pen = 25)
    #' pen25 = PELTObj$plot(main = "PELT: pen = 25")
    #' pen1 | pen25

    plot = function(d = 1L, endPts, dimNames, main, xlab, tsWidth = 0.25,
                    tsCol = "#5B9BD5",
                    bgCol = c("#A3C4F3", "#FBB1BD"),
                    bgAlpha = 0.5,
                    ncol = 1L){


      if(missing(main)){
        main = paste0("PELT: ", title_text = paste0("d = (", toString(d), ")"))

      } else {
        if(is.null(main) || !is.character(main) || length(main )!= 1L){
          stop("main must be a single character!")
        }
      }

      if(missing(xlab)){
        xlab = "Time"

      } else {
        if(is.null(xlab) || !is.character(xlab) || length(xlab)!= 1L){
          stop("xlab must be a single character!")
        }
      }

      if(missing(endPts)){
        message("endPts is missing. Proceed to use the temporary endPts!")

        if(is.null(private$.tmpEndPts)){
          stop("Temporary endPts is NULL. Must run $predict() to obtain this!")
        }
      } else{
        if(is.null(endPts) | !is.numeric(endPts)){
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

      if (is.null(d) || !is.numeric(d)) {
        stop("d must be a numeric/integer vector specifying dimensions! e.g., 1, or c(1,2).")
      }

      if(missing(dimNames)){
        message("dimNames is missing. Proceed to use the default dimNames! e.g., paste0('X', d)).")
        dimNames = paste0("X", d)
      } else {
        if (is.null(dimNames) || !is.character(dimNames)) {
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

