#' PeltL2
#'
#' An R6 class implements Pelt for offline changepoint detection,
#' currently only supports the L2 cost function.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @export
#'
#' @field tsMat A time series matrix of size $n \times p$ where
#' where each row is an observation ordered in time.
#' @field minSize An integer specifying the minimum size of a segment; by default, minSize = 1.
#' @field jump  An integer value defines the step size for the search grid during changepoint
#' detection -  only candidates changepoints \{1,k+1,2k+1,...\}  will be considered.
#' @field fitted  A boolean value specifies whether or not PeltL2 has been fitted on a dataset
#'- otherwise methods like $plot() and $predict() will return an error; by default, fitted = FALSE.
#' @field n The number of observations.
#' @field p The number of dimensions.

PeltL2 <- R6Class(
  "PeltL2",

  public = list(

    tsMat = NULL,
    minSize = NULL,
    jump = NULL,
    fitted = FALSE,
    n = NULL,
    p = NULL,

    #' Initialise a PeltL2 object
    #' @param minSize An integer specifying the minimum regime size - is 1 by default.
    #' @param jump <int> An integer determines the search space - only changepoints
    #'  \{1,k+1,2k+1,...\}  will be considered; by default jump = 1.
    #' @return Return no value but initialise a PeltL2 object with minSize and jump
    #' parameters.
    #
    initialize = function(minSize = 1, jump = 1) {
      self$minSize = minSize
      self$jump = jump
      print("You have created a PeltL2 object!")
    },

    #' Input a time series matrix for PeltL2 changepoint detection
    #' @param tsMat A time series matrix of size $n \times p$ where
    #' where each row is an observation ordered in time.
    #' @return Return no value but input tsMat, n, and p, and fitted = TRUE.
    #'
    fit = function(tsMat) { #Add error handling for tsMat
      self$tsMat = tsMat
      self$n = nrow(tsMat)
      self$p = ncol(tsMat)
      self$fitted = TRUE
    },

    #' Carry out change-point detection for a specified penalty - must be run after fitted!
    #' @param pen A non-negative double specifying a linear penaltu value for each additional
    #' change-point.
    #' @return Return no value but input tsMat, n, and p, and fitted = TRUE.
    #'
    predict = function(pen = 0){
      if(!self$fitted){
        stop("Must run $fit() before running $predict()!!!")
      }

      endPts = sort(peltL2(self$tsMat, pen, #Change the function's name to PeltCpp
                           minSize = self$minSize,
                           jump = self$jump))
      nEndPts = length(endPts)

      return(endPts[-nEndPts]) #Remove n
    },

    #' Plot PeltL2 segmentation results for a specified dimension - must be run after fitted
    #' @param pen A non-negative double specifying a linear penaltu value for each additional
    #' change-point.
    #'
    #' @return Return no value but input tsMat, n, and p, and fitted = TRUE.
    #'
    plot = function(pen, d = 1, xlab = "Iteration", ylab = "Value", #Add error handling for (pen, d)
                    main = paste0("peltSeg with pen = ", round(pen,2))){

      if(!self$fitted){
        stop("Must run $fit() before running $plot()!!!")
      }

      if(d > self$p){
        stop("d exceeds the number of dimensions in tsMat!!!")
      }

      optBkps = self$predict(pen) #optimal set of bkps
      nOptBkps = length(optBkps)
      color = rainbow(nOptBkps+1)

      ts.plot(self$tsMat[,d], xlab = xlab, ylab = ylab, main = main, col = color[1])

      for(i in 1:nOptBkps){

        lines(self$tsMat[1:optBkps[nOptBkps-i+1],d], col = color[i+1])

      }

    }

  ),

)

