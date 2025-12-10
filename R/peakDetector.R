#' peakDetector (`binSeg`)
#'
#' @description A prototype `R6` class implementing peak detector for local maxima
#' @docType class
#' @importFrom R6 R6Class is.R6
#' @importFrom ggplot2 aes ggplot geom_rect geom_line scale_fill_identity theme_minimal theme geom_vline labs element_blank element_text facet_wrap
#' @import patchwork
#' @importFrom utils hasName
#' @export
#'
#' @details
#' Lorem ipsum
#'
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new()}}{Initialises a `peakDetector` object.}
#'   \item{\code{$fit()}}{Constructs a `peakDetector` module in `R`.}
#'   \item{\code{$predict()}}{Performs `peakDetector` given a linear penalty value.}
#'   \item{\code{$plot()}}{Plots peak detection in `ggplot` style.}
#'   \item{\code{$clone()}}{Clones the `R6` object.}
#' }
#'
#' @author
#' Minh Long Nguyen \email{edelweiss611428@gmail.com}
#'
#' @export

peakDetector = R6Class(
  "peakDetector",

  public = list(

    minSize = 1L,
    window = NULL,
    tsVec = NULL,
    detectedPeaks = NULL,
    costVec = NULL,
    n = NULL,
    tempPeaks = NULL,

    #' @description Initialises a `peakDetector` object.
    #'
    #' @param minSize Integer. Minimum allowed segment length. Default: `1L`.
    #' @param radius Integer. Radius of each sliding window. Default: `NULL`, resulting in `radius = minSize`. Otherwise, `max(minSize, radius)`.
    #' @return Invisibly returns `NULL`.

    initialize = function(minSize = 1L, radius = NULL) {

      self$minSize = minSize

      if(is.null(radius)){
        self$radius = minSize
      } else{
        self$radius = max(radius, minSize)
      }

      invisible(NULL)
    },


    #' @description Constructs an `R` module for peak detection.
    #'
    #' @param tsVec Numeric vector. A time series vector size \eqn{n}
    #' @return Invisibly returns `NULL`.
    #'
    #' @details Lorem ipsum.

    fit = function(tsVec) {

      self$tsVec = tsVec
      self$n = length(tsVec)
      detectedPeaks = integer()
      costVec = numeric()
      skipID = -Inf
      for(i in self$radius:(n-self$radius)){
        if(i <= skipID) next

        ci = i
        ri = (i+1):(i + self$radius)
        li = (i - self$radius):(i-1)

        if(self$tsVec[ci] >= max(self$tsVec[ri]) & self$tsVec[ci] >= max(self$tsVec[li])){
          detectedPeaks = c(detectedPeaks, i)
          costVec = c(costVec, self$tsVec[ci])
          skipID = i + self$radius
        }

      }

      self$detectedPeaks = detectedPeaks
      self$costVec = costVec

      invisible(NULL)
    },

    #' @description Performs `peakDetection` given a linear penalty value.
    #'
    #' @param pen Numeric. Penalty per change-point. Default: `0`.
    #'
    #' @return An integer vector of detected peaks.

    predict = function(pen = 0){

      orderID = order(self$costVec, decreasing = T)
      orderedCosts = self$costVec[orderID]
      orderedPeaks = self$detectedPeaks[orderID]
      nPeaks = length(self$detectedPeaks)
      penCosts = -pen*1:nPeaks + cumsum(self$costVec)
      w.max = which.max(penCosts)
      self$tempPeaks = orderedPeaks[1:w.max]

      return(self$tempPeaks)


    },


    #' @description Plots change-point segmentation
    #'
    #' @param lipsum
    #' @return An object of classes `gg` and `ggplot`.

    plot = function(){

      ts.plot(self$tsVec)
      points(self$tempPeaks, self$tsVec[self$tempPeaks], col = "red")
    }

  )

)

