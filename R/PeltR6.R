#' Pruned Exact Linear Time (PELT)
#'
#' @description An R6 class implementing the PELT algorithm for offline change point detection.
#'
#' @docType class
#' @importFrom R6 R6Class
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
#'   \item{\code{$describe()}}{Describes a PELT object.}
#'   \item{\code{$fit()}}{Takes a time series matrix as input.}
#'   \item{\code{$predict()}}{Performs PELT given a linear penalty value.}
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
    .p = NULL

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
        stop("costFunc must be a single string!")
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
    #' @param costFunc A string specifying a cost function. Currently, only "L2" is supported.
    #'
    #' @return Invisibly returns NULL. Creates a PELT object with params minSize, jump, and costFunc.
    #'
    #' @examples
    #' peltObj = PELT$new(minSize = 1L, jump = 1L, costFunc = "L2")

    initialize = function(minSize = 1L, jump = 1L, costFunc = "L2") {
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
    #' pen = NULL, which forces pen = 2*log(n).
    #'
    #' @return A vector of indexes corresponding to the end point of each regime. By design, the last element
    #' of the vector is the number of observations.
    #'
    #' @examples
    #' peltObj = PELT$new(minSize = 1L, jump = 1L, costFunc = "L2")
    #' tsMat = as.matrix(c(rnorm(100,0), rnorm(100,5)))
    #' peltObj$fit(tsMat)
    #' peltObj$predict(pen = NULL)

    predict = function(pen = NULL){

      if(is.null(pen)){
         pen = 2*log(private$.n)
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

      return(endPts)
    }
  )

)

