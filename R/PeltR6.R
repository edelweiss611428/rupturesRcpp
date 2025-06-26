#' Pruned Exact Linear Time (PELT)
#'
#' @description An R6 class implementing the PELT algorithm for offline change point detection.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom ggplot2 aes ggplot geom_rect geom_line scale_fill_identity theme_minimal theme geom_vline labs element_blank
#' @import patchwork
#'
#' @export
#'
#' @details
#' PELT (Pruned Exact Linear Time) is an efficient algorithm for change point detection
#' that prunes the search space to achieve optimal segmentation in linear time under certain conditions.
#'
#' This implementation currently only supports the L2 cost function.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a PELT object.}
#'   \item{\code{$describe()}}{Describes the PELT object.}
#'   \item{\code{$fit()}}{Takes a time series matrix as input.}
#'   \item{\code{$predict()}}{Performs PELT given a linear penalty value.}
#'   \item{\code{$plot()}}{Plots change point segmentation in ggplot/patchwork style.}
#'   \item{\code{$clone()}}{Clone the PELT object.}
#' }
#'
#' @references
#' Truong, C., Oudre, L., & Vayatis, N. (2020). Selective review of offline change point detection methods.
#' Signal Processing, 167, 107299.
#'
#' Killick, R., Fearnhead, P., & Eckley, I. A. (2012). Optimal detection of change points with a linear computational cost.
#' Journal of the American Statistical Association, 107(500), 1590-1598.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @export

PELT = R6Class(
  "PELT",

  private = list(

    .minSize = 1L,
    .jump = 1L,
    .costFunc = "L2",
    .tsMat = NULL,
    .fitted = FALSE,
    .n = NULL,
    .p = NULL,
    .tmpEndPts = NULL, #Temporary end points - obtained after running $predict()
    .tmpPen = NULL #Temporary penalty value - obtained after running $predict()
  ),

  active = list(

    #' @field minSize An active binding. Sets the internal variable \code{.minSize} but should not be called directly.
    minSize = function(intVal) {
      if (missing(intVal)) return(private$.minSize)
      if (is.null(intVal) || !is.numeric(intVal) || as.integer(intVal) < 1 || length(intVal) != 1) {
        stop("minSize must be a single positive integer!")
      }
      private$.minSize = as.integer(intVal)
    },

    #' @field jump An active binding. Sets the internal variable \code{.jump} but should not be called directly.
    jump = function(intVal) {
      if (missing(intVal)) return(private$.jump)
      if (is.null(intVal) || !is.numeric(intVal) || as.integer(intVal) < 1 || length(intVal) != 1) {
        stop("jump must be a single positive integer!")
      }
      private$.jump = as.integer(intVal)
    },

    #' @field costFunc An active binding. Sets the internal variable \code{.costFunc} but should not be called directly.
    costFunc = function(charVal) {
      if (missing(charVal)) return(private$.costFunc)
      if (is.null(charVal) || !is.character(charVal) || length(charVal) != 1) {
        stop("costFunc must be a single character!")
      } else{
        if(!charVal %in% c("L2")){
          stop("costFunc is not support!")
        }
      }
      private$.costFunc = charVal
    },

    #' @field tsMat An active binding. Sets the internal variable \code{.tsMat} but should not be called directly.
    tsMat = function(numMat) {
      if (is.null(tsMat) || !is.numeric(numMat) || !is.matrix(numMat)) {
        stop("tsMat must be a numeric time series matrix!")
      }
      private$.tsMat = numMat
    }

  ),

  public = list(

    #' @description Initialises a PELT object.
    #'
    #' @param minSize An integer specifying the minimum segment size. By default, minSize = 1L.
    #' @param jump An integer k defining the search grid - only candidate change points in \{1,k+1,2k+1,...\}
    #' will be considered. By default, jump = 1L.
    #' @param costFunc A character specifying a cost function. Currently, only "L2" is supported.
    #'
    #' @return Invisibly returns NULL. Creates a PELT object with params minSize, jump, and costFunc.
    #'
    #' @examples
    #' peltObj = PELT$new(minSize = 1L, jump = 1L, costFunc = "L2")

    initialize = function(minSize, jump, costFunc) {
      self$minSize = minSize
      self$jump = jump
      self$costFunc = costFunc
      print("You have created a PELT object!")

      invisible(NULL)
    },

    #' @description Describes a PELT object.
    #'
    #' @return Invisibly returns a list containing the following fields of the PELT object:
    #' \describe{
    #'   \item{\code{minSize}}{The minimum segment size.}
    #'   \item{\code{jump}}{The integer k defining the search grid \{1,k+1,2k+1,...\}.}
    #'   \item{\code{costFunc}}{The cost function.}
    #'   \item{\code{fitted}}{A boolean indicating whether or not $fit() has been run.}
    #'   \item{\code{tsMat}}{The input time series matrix.}
    #'   \item{\code{n}}{The number of observations in tsMat.}
    #'   \item{\code{p}}{The number of features in tsMat.}
    #' }
    #'
    #' @examples
    #' peltObj = PELT$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' peltObj$describe()

    describe = function() {

      cat(sprintf("Pruned Exact Linear Time (PELT) \n"))
      cat(sprintf("minSize  : %sL\n", private$.minSize))
      cat(sprintf("jump     : %sL\n", private$.jump))
      cat(sprintf("costFunc : \"%s\"\n", private$.costFunc))
      cat(sprintf("fitted   : %s\n", private$.fitted))
      cat(sprintf("n        : %sL\n", private$.n))
      cat(sprintf("p        : %sL\n", private$.p))


      params = list(minSize = private$.minSize,
                    jump = private$.minSize,
                    costFunc = private$.costFunc,
                    fitted = private$.fitted,
                    tsMat = private$.tsMat,
                    n = private$.n,
                    p = private$.p)

      invisible(params)

    },

    #' @description Takes a time series matrix as input.
    #'
    #' @param tsMat tsMat A time series matrix of size \eqn{n \times p}whose rows are observations ordered in time.
    #'
    #' @return Invisibly returns NULL. Initialises .tsMat, .n, and .p, and sets private$.fitted to TRUE,
    #' enabling the use of $predict().
    #'
    #' @examples
    #' peltObj = PELT$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' peltObj$fit(tsMat)

    fit = function(tsMat) {
      self$tsMat = tsMat
      private$.n = nrow(tsMat)
      private$.p = ncol(tsMat)
      private$.fitted = TRUE #Needed for the $predict() method.
      invisible(NULL)
    },

    #' @description Performs PELT given a linear penalty value.
    #'
    #' @param pen A single non-negative numeric value specifying a penalty for each additional change point. By default,
    #' pen is missing, which forces pen = log(n).
    #'
    #' @return A vector of indexes corresponding to the end point of each regime. By design, the last element
    #' of the vector is the number of observations. Temporary end points are saved to private$.tmpEndPoints,
    #' which allows users to usethe $plot() method without specifying end points.
    #'
    #' @examples
    #' peltObj = PELT$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' peltObj$fit(tsMat)
    #' peltObj$predict()

    predict = function(pen){

      if(missing(pen)){
        pen = log(private$.n)
      }

      if(!private$.fitted){
        stop("$fit() must be run before $predict()!")
      }

      if(!is.numeric(pen) || length(pen)!= 1 || pen < 0){
        stop("pen must be a single non-negative numeric value!")
      }

      endPts = sort(peltL2(private$.tsMat,
                           pen,
                           minSize = private$.minSize,
                           jump = private$.jump))

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
    #' will use the default title "PELT: d = ...".
    #' @param xlab A character specifying the x-axis label. By default, xlab is missing, which force xlab = "Time".
    #' @param tsWidth A numeric value specifying the linewidth of the time series and also that of the segments' dashed lines.
    #' By default, tsWidth = 0.25.
    #' @param tsCol A character color of the plotted time series. By default, tsCol = "#5B9BD5".
    #' @param bgCol A character vector specifying segment colors. This will be repeated up to the length of
    #' endPts. By defaults, bgCol = c("#A3C4F3", "#FBB1BD").
    #' @param bgAlpha A numeric value specifying the degree of transparency of the background.
    #' By default, bgAlpha = 0.5.
    #'
    #' @details Plots change point segmentation results. Based on ggplot2 and patchwork. Multiple plots can easily be
    #' combined using / (vertically) and | (horizontally) operator.
    #'
    #' @return An object of classes "patchwork"/"gg"/"ggplot".
    #'
    #' @examples
    #' peltObj = PELT$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' peltObj$fit(tsMat)
    #' peltObj$predict(pen = 1)
    #' pen1 = peltObj$plot(main = "PELT: pen = 1")
    #' peltObj$predict(pen = 25)
    #' pen25 = peltObj$plot(main = "PELT: pen = 25")
    #' pen1 | pen25

    plot = function(d = 1L, endPts, dimNames, main, xlab, tsWidth = 0.25,
                    tsCol = "#5B9BD5",
                    bgCol = c("#A3C4F3", "#FBB1BD"),
                    bgAlpha = 0.5){


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
        warning("endPts is missing. Proceed to use the temporary endPts!")

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

          if(max(endPts) != n){
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
        warning("dimNames is missing. Proceed to use the default dimNames! e.g., paste0('X', d)).")
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
      segInt = data.frame(
        xmin = c(1, endPts[-length(endPts)]),
        xmax = endPts,
        fill = rep(bgCol, length.out = length(endPts))
      )

      nPlots = length(d)

      if(nPlots == 1L){

        df = data.frame(time = 1:private$.n, value = private$.tsMat[,d])

        ggplot() +
          scale_fill_identity()+
          geom_rect(data = segInt, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                    alpha = bgAlpha, inherit.aes = FALSE) +  # translucent background
          geom_line(data = df, aes(x = time, y = value), color = tsCol, linewidth = tsWidth)+
          theme_minimal()+
          theme(panel.grid = element_blank())+
          geom_vline(xintercept = endPts[-length(endPts)], linetype = "dashed", color = "black", linewidth = tsWidth) +
          labs(x = xlab, y = dimNames, title = main) -> plotObj

        return(plotObj)
      }


      for(i in 1:nPlots){

        df = data.frame(time = 1:private$.n, value = private$.tsMat[,d[i]])

        if(i == 1L){

          ggplot() +
            scale_fill_identity()+
            geom_rect(data = segInt, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                      alpha = bgAlpha, inherit.aes = FALSE) +  # translucent background
            geom_line(data = df, aes(x = time, y = value), color = tsCol, linewidth = tsWidth)+
            theme_minimal()+
            theme(panel.grid = element_blank())+
            geom_vline(xintercept = endPts[-length(endPts)], linetype = "dashed", color = "black", linewidth = tsWidth) +
            theme(
              axis.title.x = element_blank(),   # Remove x-axis title
              axis.text.x  = element_blank(),   # Remove x-axis tick labels
            ) +
            labs(y = dimNames[i], title = main) -> plotObj

        } else if (i == nPlots){

          plotObj = plotObj / ggplot() +
            scale_fill_identity()+
            geom_rect(data = segInt, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                      alpha = bgAlpha, inherit.aes = FALSE) +  # translucent background
            geom_line(data = df, aes(x = time, y = value), color = tsCol, linewidth = tsWidth)+
            theme_minimal()+
            theme(panel.grid = element_blank())+
            geom_vline(xintercept = endPts[-length(endPts)], linetype = "dashed", color = "black", linewidth = tsWidth) +
            labs(x = xlab, y = dimNames[i])

        } else {

          plotObj = plotObj/ggplot() +
            scale_fill_identity()+
            geom_rect(data = segInt, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                      alpha = bgAlpha, inherit.aes = FALSE) +  # translucent background
            geom_line(data = df, aes(x = time, y = value), color = tsCol, linewidth = tsWidth)+
            theme_minimal()+
            theme(panel.grid = element_blank())+
            geom_vline(xintercept = endPts[-length(endPts)], linetype = "dashed", color = "black", linewidth = tsWidth) +
            theme(
              axis.title.x = element_blank(),   # Remove x-axis title
              axis.text.x  = element_blank(),   # Remove x-axis tick labels
            ) +
            labs(y = dimNames[i])

        }

      }

      return(plotObj)

    }
  )

)

