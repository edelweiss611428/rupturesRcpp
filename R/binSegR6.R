#' @importFrom R6 R6Class
#' @export

BinSegL2 <- R6Class(
  "BinSegL2",

  public = list(
    nr = NULL,
    cp = NULL,
    tsMat = NULL,
    k = NULL,
    fitted = FALSE,
    output = NULL,
    cost = NULL,

    initialize = function(tsMat, k) {
      self$tsMat = tsMat
      self$nr = nrow(tsMat)
      self$k = k
      print("You have created a BinSegL2 object!")
    },

    fit = function() {
      self$output = binSegCpp(self$tsMat, self$k)
      self$cp = self$output$changePoints
      self$cost = self$output$cost
      self$fitted = TRUE
    },

    plot = function(what = NULL,
                    nRegimes = NULL, #number of changepoints
                    main = NULL) {

      if(is.null(main)){
        main = "binSeg clustering"
      }

      if(is.null(what)){

        if(ncol(self$tsMat == 1)){ #Current version only plot the k-partition

          if(!self$fitted){
            ts.plot(self$tsMat, ylab = "X", main = main)
            warning("Should have run the fit() method first!")
          } else if (is.null(nRegimes)){

            ts.plot(self$tsMat, ylab = "X",
                    main = main)

            sortedRegimes = c(sort(self$cp), self$nr)
            colors = rainbow(self$k)
            for(i in self$k:1){

              lines(self$tsMat[1:sortedRegimes[i]], col = colors[i])

            }

          } else if (is.numeric(nRegimes) &
                     length(nRegimes) == 1){

            if(nRegimes == 1){
              ts.plot(self$tsMat, ylab = "X",
                      main = main)
            } else if(nRegimes <= self$k){

              tempCPts = self$cp[1:(nRegimes-1)]
              ts.plot(self$tsMat, ylab = "X",
                      main = main)

              sortedRegimes = c(sort(tempCPts), self$nr)
              colors = rainbow(nRegimes)
              for(i in nRegimes:1){

                lines(self$tsMat[1:sortedRegimes[i]], col = colors[i])

              }

            } else{
              print("Invalid nRegimes!")
            }
          } else {
            print("Invalid or unsupported method!")
          }

        } else {
          print("Currently does not support high-dimensional plots!")
        }

      } else if (what == "cost"){

        if(!self$fitted){
          stop("Should have run the fit() method first!")
        } else {
          ts.plot(self$cost, xlab = "Number of regimes", ylab = "Cost")
        }

      }

    }

  ),
)

