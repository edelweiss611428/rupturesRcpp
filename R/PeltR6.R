#' @importFrom R6 R6Class
#' @export

library("R6")

PeltL2 <- R6Class(
  "PeltL2",

  public = list(

    tsMat = NULL,
    minSize = NULL,
    jump = NULL,
    fitted = FALSE,
    n = NULL,
    p = NULL,

    initialize = function(minSize = 1, jump = 1) {
      self$minSize = minSize
      self$jump = jump
      print("You have created a PeltL2 object!")
    },

    fit = function(tsMat) { #Add error handling for tsMat
      self$tsMat = tsMat
      self$n = nrow(tsMat)
      self$p = ncol(tsMat)
      self$fitted = TRUE
    },

    predict = function(pen){
      if(!self$fitted){
        stop("Must run $fit() before running $predict()!!!")
      }

      endPts = sort(peltL2(self$tsMat, pen, #Change the function's name to PeltCpp
                           minSize = self$minSize,
                           jump = self$jump))
      nEndPts = length(endPts)

      return(endPts[-nEndPts]) #Remove n
    },

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

