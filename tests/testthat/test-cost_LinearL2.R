#LinearL2 module

# ========================================================
#               (R) LinearL2 cost function
# ========================================================

R_LinearL2eval = function(Y, X, start, end){

  subX = X[(start+1):end,]
  subY = Y[(start+1):end,]

  return(sum(lm(subY~subX)$residuals^2))
}

# ========================================================
#                   Simulated datasets
# ========================================================

set.seed(12345)
X = cbind(rnorm(100), rnorm(100))
Y = X*cbind(rep(c(1,5), each = 50), rep(c(1,5), each = 50)) + cbind(rnorm(100),rnorm(100))
nr = nrow(X)

nCases = 10
idx1 = sample.int(nr-4, nCases)
idx2 = integer(nCases)

for(i in 1:nCases){
  repeat{
    idx2[i] = sample((idx1[i]+1):nr, 1)

    if(idx2[i] > idx1[i] + 2){
      break #Ensure there is no singularity
    }
  }
}


test_that("[No singularity] Expect C++ .eval() method in LinearL2 cost module gives the correct results", {


  LinearL2module = new(rupturesRcpp:::Cost_LinearL2, Y, X, TRUE, TRUE) # (intercept_, warnOnce) = (TRUE, TRUE)

  for(i in 1:nCases){
    gs = R_LinearL2eval(Y,X, idx1[i],idx2[i])
    expect_equal(LinearL2module$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(LinearL2module$eval(0,nr+1),
               regexp = "out of bounds") #arma::mat indexing error

  expect_error(LinearL2module$eval(-1,nr),
               regexp = "out of bounds") #arma::mat indexing error

  #Give 0 if end - start = 0 or 1
  expect_equal(LinearL2module$eval(0,0), 0)
  expect_equal(LinearL2module$eval(0,1), 0)

  expect_no_error(LinearL2module$resetWarning(FALSE)) #Generally does nothing here

})



test_that("[No singularity] Expect $eval() method in PELT_LinearL2 gives the correct results/error message", {

  set.seed(12345)

  PELTObj = PELT$new(costFunc = costFunc$new("LinearL2"))
  PELTObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearL2eval(Y,X, idx1[i],idx2[i])
    expect_equal(PELTObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(PELTObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(PELTObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(PELTObj$eval(0,0),
               regexp = "a must be smaller than b")

  #Give 0 if end - start = 1
  expect_equal(PELTObj$eval(0,1), 0)


})


test_that("[No singularity] Expect .eval() method in C++ PELT_LinearL2 class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, int, int>()
  PELTCppObj = new(PELTCpp_LinearL2, Y, X, TRUE, 1L, 1L)

  for(i in 1:nCases){
    gs = R_LinearL2eval(Y,X, idx1[i],idx2[i])
    expect_equal(PELTCppObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(PELTCppObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(PELTCppObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(PELTCppObj$eval(0,0),
               regexp = "start < end` must be true!")

  #Give 0 if end - start = 1
  expect_equal(PELTCppObj$eval(0,1), 0)


})


test_that("[No singularity] Expect $eval() method in binSeg_LinearL2 gives the correct results/error message", {

  binSegObj = binSeg$new(costFunc = costFunc$new("LinearL2"))
  binSegObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearL2eval(Y,X, idx1[i],idx2[i])
    expect_equal(binSegObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(binSegObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(binSegObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(binSegObj$eval(0,0),
               regexp = "a must be smaller than b")

  #Give 0 if end - start = 1
  expect_equal(binSegObj$eval(0,1), 0)


})



test_that("[No singularity] Expect .eval() method in C++ binSeg_LinearL2 class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, int, int>()
  binSegCppObj = new(binSegCpp_LinearL2, Y, X, TRUE, 1L, 1L)

  for(i in 1:nCases){
    gs = R_LinearL2eval(Y,X, idx1[i],idx2[i])
    expect_equal(binSegCppObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(binSegCppObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(binSegCppObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(binSegCppObj$eval(0,0),
               regexp = "start < end` must be true!")

  #Give 0 if end - start = 1
  expect_equal(binSegCppObj$eval(0,1), 0)


})



test_that("[No singularity] Expect $eval() method in window_LinearL2 gives the correct results/error message", {

  windowObj = Window$new(costFunc = costFunc$new("LinearL2"))
  windowObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearL2eval(Y,X, idx1[i],idx2[i])
    expect_equal(windowObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(windowObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(windowObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(windowObj$eval(0,0),
               regexp = "a must be smaller than b")

  #Give 0 if end - start = 1
  expect_equal(windowObj$eval(0,1), 0)


})



test_that("[No singularity] Expect .eval() method in C++ window_LinearL2 class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, int, int, int>()
  windowCppObj = new(windowCpp_LinearL2, Y, X, TRUE, 1L, 1L, 10L)
  idx1 = sample.int(nr-3, nCases)
  idx2 = integer(nCases)

  for(i in 1:nCases){
    repeat{
      idx2[i] = sample((idx1[i]+1):nr, 1)

      if(idx2[i] > idx1[i] + 1){
        break #Generally, ensure there is no singularity
      }
    }
  }

  for(i in 1:10){
    gs = R_LinearL2eval(Y,X, idx1[i],idx2[i])
    expect_equal(windowCppObj$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(windowCppObj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")

  expect_error(windowCppObj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

  expect_error(windowCppObj$eval(0,0),
               regexp = "start < end` must be true!")

  #Give 0 if end - start = 1
  expect_equal(windowCppObj$eval(0,1), 0)


})



test_that("Expect correct error/warning messages when initialising C++ Cost_LinearL2 module", {

  expect_error(new(rupturesRcpp:::Cost_LinearL2, Y, as.matrix(X[-1]), TRUE, TRUE), "Number of observations in response and covariate matrices must match!")

  X2 = matrix(1:5, nrow = 1)
  Y2 = matrix(1, nrow = 1)

  expect_error(new(rupturesRcpp:::Cost_LinearL2, Y2, X2, TRUE, TRUE), "The full dataset contains not enough observations to fit a linear regression model!")

  Xconst = matrix(rep(1,100))
  Y3 = matrix(rnorm(100))

  LinearL2module = new(rupturesRcpp:::Cost_LinearL2, Y3, Xconst, TRUE, TRUE)
  #warnOnce_ = TRUE -> keepWarning = FALSE (expected behavior when running eval() inside segmentation methods)
  expect_warning(LinearL2module$eval(0, 100), "System is singular")

  #warnOnce_ = FALSE -> keepWarning = TRUE (expected behavior when running eval() outside segmentation methods)
  LinearL2module = new(rupturesRcpp:::Cost_LinearL2, Y3, Xconst, TRUE, FALSE)
  expect_warning(LinearL2module$eval(0, 100), "Singular system encountered")

})







