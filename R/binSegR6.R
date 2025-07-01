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

    #' @field minSize An active binding. Sets the internal variable \code{.minSize} but should not be called directly.
    minSize = function(intVal) {
      if (is.null(intVal) || !is.numeric(intVal) || as.integer(intVal) < 1 || length(intVal) != 1) {
        stop("minSize must be a single positive integer!")
      }
      private$.minSize = as.integer(intVal)
    },

    #' @field jump An active binding. Sets the internal variable \code{.jump} but should not be called directly.
    jump = function(intVal) {
      if (is.null(intVal) || !is.numeric(intVal) || as.integer(intVal) < 1 || length(intVal) != 1) {
        stop("jump must be a single positive integer!")
      }
      private$.jump = as.integer(intVal)
    },

    #' @field costFunc An active binding. Sets the internal variable \code{.costFunc} but should not be called directly.
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

    #' @field tsMat An active binding. Sets the internal variable \code{.tsMat} but should not be called directly.
    tsMat = function(numMat) {
      if (is.null(numMat) || !is.numeric(numMat) || !is.matrix(numMat)) {
        stop("tsMat must be a numeric time series matrix!")
      }
      private$.tsMat = numMat
    },

    #' @field addSmallDiag An active binding. Sets the internal variable \code{.addSmallDiag} but should not be called directly.
    addSmallDiag = function(boolVal) {
      if (is.null(boolVal) || !is.logical(boolVal) || length(boolVal) != 1) {
        stop("addSmallDiag must be a single boolean value!")
      }
      private$.addSmallDiag = boolVal
    },

    #' @field epsilon An active binding. Sets the internal variable \code{.epsilon} but should not be called directly.
    epsilon = function(doubleVal) {
      if (is.null(doubleVal) || !is.numeric(doubleVal) || as.integer(doubleVal) < 0 || length(doubleVal) != 1) {
        stop("epsilon must be a single positive double!")
      }
      private$.epsilon = doubleVal
    },

    #' @field pVAR An active binding. Sets the internal variable \code{.pVAR} but should not be called directly.
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
    #' @param minSize An integer specifying the minimum segment size. By default, minSize = 1L.
    #' @param jump An integer k defining the search grid - only candidate change points in \{1,k+1,2k+1,...\}
    #' will be considered. By default, jump = 1L.
    #' @param costFunc A character specifying a cost function. Currently, only "L2", "SIGMA", and "VAR" are supported. By default,
    #' costFunc = "L2".
    #' @param addSmallDiag (SIGMA) An boolean value indicating whether or not to add a small bias to the diagonal entries
    #' of estimated covariance matrices (for costFunc "SIGMA"). This improves numerical stability in near-singularity
    #' scenarios. By default, addSmallDiag = TRUE.
    #' @param epsilon (SIGMA) A double value specifying a bias value to be added to the diagonal entries
    #' of estimated covariance matrices. By default, epsilon = 10^-6.
    #' @param pVAR (VAR) A non-negative integer specify vector-autoregressive order. By default, pVAR = 1.
    #'
    #' @return Invisibly returns NULL. Creates a PELT object with params minSize, jump, and costFunc.
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
    #' @return Invisibly returns a list containing the following fields of the binSeg object:
    #' \describe{
    #'   \item{\code{minSize}}{The minimum segment size.}
    #'   \item{\code{jump}}{The integer k defining the search grid \{1,k+1,2k+1,...\}.}
    #'   \item{\code{costFunc}}{The cost function.}
    #'   \item{\code{addSmallDiag}}{A boolean value indicating whether to add a bias value to diagonal entries of estimated covariance matrices.}
    #'   \item{\code{epsilon}}{The bias value to be added to diagonal entries of estimated covariance matrices.}
    #'   \item{\code{fitted}}{A boolean indicating whether or not $fit() has been run.}
    #'   \item{\code{tsMat}}{The input time series matrix.}
    #'   \item{\code{n}}{The number of observations in tsMat.}
    #'   \item{\code{p}}{The number of features in tsMat.}
    #'   \item{\code{bkps}}{A vector containing unordered breakpoint positions, given the maximum number of allowable breakpoints.}
    #'   \item{\code{cost}}{From the second elemement, the corresponding costs of splits in bkps. By default, the
    #'   first element is the total cost without any split}
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
    #' @param tsMat tsMat A time series matrix of size \eqn{n \times p} whose rows are observations ordered in time.
    #'
    #' @return Invisibly returns NULL. Initialises .tsMat, .n, and .p, and sets private$.fitted to TRUE,
    #' enabling the use of $predict().
    #'
    #' @examples
    #' binSegObj = binSeg$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' binSegObj$fit(tsMat)

    fit = function(tsMat) {

      self$tsMat = tsMat
      private$.n = nrow(tsMat)
      private$.p = ncol(tsMat)
      private$.fitted = TRUE #Needed for the $predict() method.
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
    #' @param pen A single non-negative numeric value specifying a penalty for each additional change point. By default,
    #' pen = 0.
    #'
    #' @return A vector of indexes corresponding to the end point of each regime. By design, the last element
    #' of the vector is the number of observations. Temporary end points are saved to private$.tmpEndPoints,
    #' which allows users to usethe $plot() method without specifying end points.
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
    #' @param d An integer vector specifying dimensions to plot. By default, d = 1L.
    #' @param endPts An integer vector specifying end points! Could be obtained via $predict().").
    #' By default, endPts is missing, in which, the method will proceed to use the (latest) temporary
    #' changepoints obtained by running $predict().
    #' @param dimNames A character vector specifying feature names!", whose length must match that of d if not missing.
    #' By default, dimNames is missing, which forces dimNames = ("X1", "X2",...).
    #' @param main A character specifying the main title of the plot. By default, main is missing, in which the method
    #' will use the default title "binSeg: d = ...".
    #' @param xlab A character specifying the x-axis label. By default, xlab is missing, which force xlab = "Time".
    #' @param tsWidth A numeric value specifying the linewidth of the time series and also that of the segments' dashed lines.
    #' By default, tsWidth = 0.25.
    #' @param tsCol A character color of the plotted time series. By default, tsCol = "#5B9BD5".
    #' @param bgCol A character vector specifying segment colors. This will be repeated up to the length of
    #' endPts. By defaults, bgCol = c("#A3C4F3", "#FBB1BD").
    #' @param bgAlpha A numeric value specifying the degree of transparency of the background.
    #' By default, bgAlpha = 0.5.
    #' @param ncol An integer specifying the number of columns to be used in the facet layout. By default, ncol  = 1L.
    #'
    #' @details Plots change point segmentation results. Based on ggplot2. Multiple plots can easily be
    #' combined using / (vertically) and | (horizontally) operator via patchwork.
    #'
    #' @return An object of classes "gg"/"ggplot".
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

