#' Binary Segmentation (binSeg)
#'
#' @description An R6 class implementing binary segmentation for offline change point detection.
#'
#'
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
#' Currently supported cost functions:
#'
#' - `"L2"`: for (independent) piecewise Gaussian process with **constant covariance**
#' - `"SIGMA"`: for (independent) piecewise Gaussian process with **varying covariance**
#' - `"VAR"`: for piecewise Gaussian vector-regressive process with **constant covariance**
#'
#' See **Methods** and section for more details.
#'
#' @examples
#' # Toy example
#' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,0, 10)))
#' # Initialise a binSeg object and fit the method to tsMat
#' binSegObj = binSeg$new(costFunc = "SIGMA")
#' binSegObj$fit(tsMat) #Need to run this before running $predict()
#' # Perform binSeg for a specific linear penalty threshold
#' binSegObj$predict(pen = 50)
#' # Plot the latest segmentation solution
#' binSegObj$plot(main = "binSeg:SIGMA:pen=50", ncol = 1)
#' # Describe the binSeg object (and invisibly return the object's fields)
#' binSegObj$describe()
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a binSeg object.}
#'   \item{\code{$describe()}}{Describes the binSeg object.}
#'   \item{\code{$fit()}}{Takes a time series matrix as input and perform binSeg for the
#' maximum number of change points.}
#'   \item{\code{$predict()}}{Performs binSeg given a linear penalty value.}
#'   \item{\code{$plot()}}{Plots change point segmentation in ggplot style.}
#'   \item{\code{$clone()}}{Clone the binSeg object.}
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
    .costFunc = "L2",
    .tsMat = NULL,
    .fitted = FALSE,
    .n = NULL,
    .p = NULL,
    .addSmallDiag = TRUE,
    .epsilon = 10^-6,
    .pVAR = 1L,
    .bkps = NULL,
    .cost = NULL,
    .tmpEndPts = NULL, #Temporary end points - obtained after running $predict()
    .tmpPen = NULL #Temporary penalty value - obtained after running

  ),

  active = list(

    #' @field minSize Active binding. Sets the internal variable \code{.minSize} but should not be called directly.
    minSize = function(intVal) {
      if (is.null(intVal) || !is.numeric(intVal) || as.integer(intVal) < 1 || length(intVal) != 1) {
        stop("minSize must be a single positive integer!")
      }
      private$.minSize = as.integer(intVal)
    },

    #' @field jump Active binding. Sets the internal variable \code{.jump} but should not be called directly.
    jump = function(intVal) {
      if (is.null(intVal) || !is.numeric(intVal) || as.integer(intVal) < 1 || length(intVal) != 1) {
        stop("jump must be a single positive integer!")
      }
      private$.jump = as.integer(intVal)
    },

    #' @field costFunc Active binding. Sets the internal variable \code{.costFunc} but should not be called directly.
    costFunc = function(charVal) {
      if (is.null(charVal) || !is.character(charVal) || length(charVal) != 1) {
        stop("costFunc must be a single character!")
      } else{
        if(!charVal %in% c("L2", "SIGMA", "VAR")){
          stop("costFunc is not support!")
        }
      }
      private$.costFunc = charVal
    },

    #' @field tsMat Active binding. Sets the internal variable \code{.tsMat} but should not be called directly.
    tsMat = function(numMat) {
      if (is.null(numMat) || !is.numeric(numMat) || !is.matrix(numMat)) {
        stop("tsMat must be a numeric time series matrix!")
      }
      private$.tsMat = numMat
    },

    #' @field addSmallDiag Active binding. Sets the internal variable \code{.addSmallDiag} but should not be called directly.
    addSmallDiag = function(boolVal) {
      if (is.null(boolVal) || !is.logical(boolVal) || length(boolVal) != 1) {
        stop("addSmallDiag must be a single boolean value!")
      }
      private$.addSmallDiag = boolVal
    },

    #' @field epsilon Active binding. Sets the internal variable \code{.epsilon} but should not be called directly.
    epsilon = function(doubleVal) {
      if (is.null(doubleVal) || !is.numeric(doubleVal) || as.integer(doubleVal) < 0 || length(doubleVal) != 1) {
        stop("epsilon must be a single positive double!")
      }
      private$.epsilon = doubleVal
    },

    #' @field pVAR Active binding. Sets the internal variable \code{.pVAR} but should not be called directly.
    pVAR = function(intVal) {
      if (is.null(intVal) || !is.numeric(intVal) || as.integer(intVal) < 1 || length(intVal) != 1) {
        stop("pVAR must be a single positive integer!")
      }
      private$.pVAR = as.integer(intVal)
    }

  ),

  public = list(

    #' @description Initialises a binSeg object.
    #'
    #' @param minSize Integer. Minimum allowed segment length. Default: 1L.
    #' @param jump Integer. Search grid step size: only positions in \{1, k+1, 2k+1, ...\} are considered. Default: 1L.
    #' @param costFunc Character. Cost function to use: one of `"L2"`, `"SIGMA"`, or `"VAR"`. Default: `"L2"`.
    #' @param addSmallDiag Logical. (SIGMA) If `TRUE`, add a small value to the diagonal of estimated covariance matrices
    #' to improve numerical stability. Default: `TRUE`.
    #' @param epsilon Double. (SIGMA) A small positive value used to the diagonal of estimated covariance matrices to stabilise
    #' matrix operations. Default: `1e-6`.
    #' @param pVAR Integer (VAR). Order of the vector autoregressive (VAR) model. Must be non-negative. Default: `1L`.
    #'
    #' @return Invisibly returns NULL.
    #'
    #' @examples
    #' peltObj = PELT$new(minSize = 1L, jump = 1L, costFunc = "L2")

    initialize = function(minSize, jump, costFunc, addSmallDiag, epsilon, pVAR) {

      if(!missing(minSize)){
        self$minSize = minSize
      }

      if(!missing(jump)){
        self$jump = jump
      }

      if(!missing(costFunc)){
        self$costFunc = costFunc
      }

      if(!missing(epsilon)){
        self$epsilon = epsilon
      }

      if (!missing(addSmallDiag)){
        self$addSmallDiag = addSmallDiag
      }

      if (!missing(pVAR)){
        self$pVAR = pVAR
      }

      print("You have created a binSeg object!")

      invisible(NULL)
    },

    #' @description Describes a binSeg object.
    #'
    #' @return Invisibly returns a list containing some of the following fields:
    #' \describe{
    #'   \item{\code{minSize}}{Minimum allowed segment length.}
    #'   \item{\code{jump}}{Search grid step size.}
    #'   \item{\code{costFunc}}{The cost function.}
    #'   \item{\code{addSmallDiag}}{(SIGMA only) Whether to add bias for numerical stability.}
    #'   \item{\code{epsilon}}{(SIGMA, VAR) Bias added to diagonal entries.}
    #'   \item{\code{pVAR}}{(VAR) VAR order.}
    #'   \item{\code{fitted}}{Whether or not `$fit()` has been run.}
    #'   \item{\code{tsMat}}{Input time series matrix.}
    #'   \item{\code{n}}{Number of observations.}
    #'   \item{\code{p}}{Number of features.}
    #'   \item{\code{bkps}}{Vector of unordered breakpoint positions.}
    #'   \item{\code{cost}}{From the second element, split costs corresponding to \code{bkps}; first element is total cost without splits.}
    #' }
    #'
    #' @examples
    #' binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' binSegObj$describe()
    #'
    describe = function() {

      params = list(minSize = private$.minSize,
                    jump = private$.minSize,
                    costFunc = private$.costFunc,
                    fitted = private$.fitted,
                    tsMat = private$.tsMat,
                    n = private$.n,
                    p = private$.p,
                    bkps = private$.bkps,
                    cost = private$.cost)

      cat(sprintf("Binary Segmentation (binSeg) \n"))
      cat(sprintf("minSize      : %sL\n", private$.minSize))
      cat(sprintf("jump         : %sL\n", private$.jump))
      cat(sprintf("costFunc.    : \"%s\"\n", private$.costFunc))

      if(private$.costFunc == "SIGMA"){
        cat(sprintf("addSmallDiag : %s\n", private$.addSmallDiag))
        cat(sprintf("epsilon      : %s\n", private$.epsilon))
        params[["addSmallDiag"]] = private$.addSmallDiag
        params[["epsilon"]] = private$.epsilon
      }

      if(private$.costFunc == "VAR"){
        cat(sprintf("pVAR.        : %s\n", private$.epsilon))
        params[["pVAR"]] = private$.pVAR
      }

      cat(sprintf("fitted       : %s\n", private$.fitted))
      cat(sprintf("n            : %sL\n", private$.n))
      cat(sprintf("p            : %sL\n", private$.p))

      invisible(params)

    },

    #' @description Takes a time series matrix as input and perform binSeg for the
    #' maximum number of change points.
    #'
    #' @param tsMat Numeric matrix. A time series matrix of size \eqn{n \times p} whose rows are observations ordered in time.
    #'
    #' @return Invisibly returns NULL.
    #'
    #' @details Initialises `private$.tsMat`, `private$.n`, `private$.p`, `private$.bkps`, and `private$.cost`. Sets private$.fitted to TRUE,
    #' enabling the use of `$predict()`. Run `$describe()` for detailed configurations.
    #'
    #' @examples
    #' binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' binSegObj$fit(tsMat)

    fit = function(tsMat) {

      self$tsMat = tsMat
      private$.n = nrow(tsMat)
      private$.p = ncol(tsMat)
      private$.fitted = TRUE
      detection = binSegCpp(private$.tsMat, private$.minSize, private$.jump,
                            costFunc = private$.costFunc,
                            addSmallDiag = private$.addSmallDiag,
                            epsilon = private$.epsilon,
                            pVAR = private$.pVAR)
      private$.cost = detection$cost
      private$.bkps = detection$bkps
      invisible(NULL)
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
    #'
    #' @examples
    #' binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' binSegObj$fit(tsMat)
    #' binSegObj$predict()

    predict = function(pen = 0){

      if(!private$.fitted){
        stop("$fit() must be run before $predict()!")
      }

      if(!is.numeric(pen) || length(pen)!= 1 || pen < 0){
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
    #' @details Plots change point segmentation results. Based on ggplot2. Multiple plots can easily be
    #' horizontally and vertically stacked using patchwork's operators `/` and `|`, respectively.
    #'
    #' @return An object of classes "gg" and "ggplot".
    #'
    #' @examples
    #' binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' binSegObj$fit(tsMat)
    #' binSegObj$predict(pen = 1)
    #' pen1 = binSegObj$plot(main = "binSeg: pen = 1")
    #' binSegObj$predict(pen = 25)
    #' pen25 = binSegObj$plot(main = "binSeg: pen = 25")
    #' pen1 | pen25

    plot = function(d = 1L, endPts, dimNames, main, xlab, tsWidth = 0.25,
                    tsCol = "#5B9BD5",
                    bgCol = c("#A3C4F3", "#FBB1BD"),
                    bgAlpha = 0.5,
                    ncol = 1L){


      if(missing(main)){
        main = paste0("binSeg: ", title_text = paste0("d = (", toString(d), ")"))

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

