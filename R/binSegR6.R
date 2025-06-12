#' @importFrom R6 R6Class
#' @export

library("R6")
#Current only supports L2 metric
BinSegL2 <- R6Class(
  "BinSegL2",

  public = list(

    tsMat = NULL,
    minSize = NULL,
    jump = NULL,
    fitted = FALSE,
    bkps = NULL,
    cost = NULL,
    n = NULL,
    p = NULL,

    initialize = function(minSize = 1, jump = 1) {
      self$minSize = minSize
      self$jump = jump
      print("You have created a BinSegL2 object!")
    },

    fit = function(tsMat) { #Add error handling for tsMat
      self$tsMat = tsMat
      self$n = nrow(tsMat)
      self$p = ncol(tsMat)
      detection = binSegCpp(self$tsMat, self$minSize, self$jump)
      self$bkps = detection$bkps
      self$cost = detection$cost
      self$fitted = TRUE
    },

    predict = function(pen){
      if(!self$fitted){
        stop("Must run $fit() before running $predict()!!!")
      }
      return(sort(binSegPredCpp(self$bkps, self$cost, pen)))
    },

    plot = function(pen, d = 1, xlab = "Iteration", ylab = "Value", #Add error handling for (pen, d)
                    main = paste0("Optimal binSeg with pen = ", round(pen,2))){

      if(!self$fitted){
        stop("Must run $fit() before running $plot()!!!")
      }

      if(d > self$p){
        stop("d exceeds the number of dimensions in tsMat!!!")
      }

      optBkps = self$predict(pen) #Optimal regimes
      nOptBkps = length(optBkps)
      color = rainbow(nOptBkps+1)

      ts.plot(self$tsMat[,d], xlab = xlab, ylab = ylab, main = main, col = color[1])

      for(i in 1:nOptBkps){

        lines(self$tsMat[1:optBkps[nOptBkps-i+1],d], col = color[i+1])

      }




    }

  ),
)

