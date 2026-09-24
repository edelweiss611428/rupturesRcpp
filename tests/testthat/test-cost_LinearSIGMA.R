#LinearSIGMA module

# ========================================================
#             (R) LinearSIGMA cost function
# ========================================================

R_LinearSIGMAeval = function(Y, X, start, end, addSmallDiag = TRUE, epsilon = 1e-6){

  subX = X[(start+1):end, , drop = FALSE]
  subY = Y[(start+1):end, , drop = FALSE]
  n = end - start

  resid = as.matrix(residuals(lm(subY~subX)))
  covMat = crossprod(resid)/n

  if(addSmallDiag){
    covMat = covMat + diag(epsilon, ncol(resid))
  }

  return(n*log(det(covMat)))
}

# ========================================================
#                   Simulated datasets
# ========================================================

set.seed(12345)
X = cbind(rnorm(100), rnorm(100))
Y = X*cbind(rep(c(1,5), each = 50), rep(c(1,5), each = 50)) + cbind(rnorm(100),rnorm(100))
nr = nrow(X)
J = ncol(X) + 1L #with intercept

nCases = 10
idx1 = sample.int(nr-6, nCases)
idx2 = integer(nCases)

for(i in 1:nCases){
  repeat{
    idx2[i] = sample((idx1[i]+1):nr, 1)

    if(idx2[i] > idx1[i] + J + 1){
      break #Ensure there are enough observations to fit both the regression and covariance
    }
  }
}


test_that("Expect C++ .eval() method in LinearSIGMA cost module gives the correct results", {

  LinearSIGMAmodule = new(rupturesRcpp:::Cost_LinearSIGMA, Y, X, TRUE, TRUE, 1e-6, TRUE) # (intercept, addSmallDiag, epsilon, warnOnce)

  for(i in 1:nCases){
    gs = R_LinearSIGMAeval(Y,X, idx1[i],idx2[i])
    expect_equal(LinearSIGMAmodule$eval(idx1[i],idx2[i]), gs)
  }

  expect_error(LinearSIGMAmodule$eval(0,nr+1),
               regexp = "out of bounds") #arma::mat indexing error

  expect_error(LinearSIGMAmodule$eval(-1,nr),
               regexp = "out of bounds") #arma::mat indexing error

  #Give 0 if end - start = 0, 1, or fewer than J observations
  expect_equal(LinearSIGMAmodule$eval(0,0), 0)
  expect_equal(LinearSIGMAmodule$eval(0,1), 0)
  expect_equal(LinearSIGMAmodule$eval(0,J-1), 0)

  expect_no_error(LinearSIGMAmodule$resetWarning(FALSE)) #Generally does nothing here

})


test_that("Expect $get_params() method in LinearSIGMA cost module gives the correct results", {

  LinearSIGMAmodule = new(rupturesRcpp:::Cost_LinearSIGMA, Y, X, TRUE, TRUE, 1e-6, TRUE)

  for(i in 1:nCases){

    fit = lm(Y[(idx1[i]+1):idx2[i],]~X[(idx1[i]+1):idx2[i],])
    gsCoef = unname(coef(fit))
    resid = as.matrix(residuals(fit))
    gsCov = crossprod(resid)/(idx2[i]-idx1[i]) + diag(1e-6, ncol(resid))

    params = LinearSIGMAmodule$get_params(idx1[i],idx2[i])
    expect_equal(unname(params$coef), gsCoef)
    expect_equal(unname(params$cov), unname(gsCov))
  }

  #NA coef/cov if segment has fewer usable observations than coefficients
  params = LinearSIGMAmodule$get_params(0, J-1)
  expect_true(all(is.na(params$coef)))
  expect_true(all(is.na(params$cov)))

})


test_that("Expect $eval() method in PELT_LinearSIGMA gives the correct results/error message", {

  set.seed(12345)

  PELTObj = PELT$new(costFunc = costFunc$new("LinearSIGMA"))
  PELTObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearSIGMAeval(Y,X, idx1[i],idx2[i])
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


test_that("Expect .eval() method in C++ PELT_LinearSIGMA class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, bool, double, int, int>()
  PELTCppObj = new(PELTCpp_LinearSIGMA, Y, X, TRUE, TRUE, 1e-6, 1L, 1L)

  for(i in 1:nCases){
    gs = R_LinearSIGMAeval(Y,X, idx1[i],idx2[i])
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


test_that("Expect $eval() method in binSeg_LinearSIGMA gives the correct results/error message", {

  binSegObj = binSeg$new(costFunc = costFunc$new("LinearSIGMA"))
  binSegObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearSIGMAeval(Y,X, idx1[i],idx2[i])
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


test_that("Expect .eval() method in C++ binSeg_LinearSIGMA class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, bool, double, int, int>()
  binSegCppObj = new(binSegCpp_LinearSIGMA, Y, X, TRUE, TRUE, 1e-6, 1L, 1L)

  for(i in 1:nCases){
    gs = R_LinearSIGMAeval(Y,X, idx1[i],idx2[i])
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


test_that("Expect $eval() method in window_LinearSIGMA gives the correct results/error message", {

  windowObj = Window$new(costFunc = costFunc$new("LinearSIGMA"))
  windowObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearSIGMAeval(Y,X, idx1[i],idx2[i])
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


test_that("Expect .eval() method in C++ window_LinearSIGMA class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, bool, double, Rcpp::IntegerVector>() -- last arg is c(minSize, jump, radius)
  windowCppObj = new(windowCpp_LinearSIGMA, Y, X, TRUE, TRUE, 1e-6, c(1L, 1L, 10L))

  idx1b = sample.int(nr-6, nCases)
  idx2b = integer(nCases)

  for(i in 1:nCases){
    repeat{
      idx2b[i] = sample((idx1b[i]+1):nr, 1)

      if(idx2b[i] > idx1b[i] + J + 1){
        break
      }
    }
  }

  for(i in 1:nCases){
    gs = R_LinearSIGMAeval(Y,X, idx1b[i],idx2b[i])
    expect_equal(windowCppObj$eval(idx1b[i],idx2b[i]), gs)
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


test_that("Expect correct error/warning messages when initialising C++ Cost_LinearSIGMA module", {

  expect_error(new(rupturesRcpp:::Cost_LinearSIGMA, Y, as.matrix(X[-1]), TRUE, TRUE, 1e-6, TRUE),
               "Number of observations in response and covariate matrices must match!")

  X2 = matrix(1:5, nrow = 1)
  Y2 = matrix(1, nrow = 1)

  expect_error(new(rupturesRcpp:::Cost_LinearSIGMA, Y2, X2, TRUE, TRUE, 1e-6, TRUE),
               "The full dataset contains not enough observations to fit a linear regression model!")

  #Singular design matrix (X duplicates the intercept column) -> shares RegressionCost::solveSegment
  #fallback/warning logic with Cost_LinearL2 and Cost_VAR
  Xconst = matrix(rep(1,100))
  Y3 = cbind(rnorm(100), rnorm(100))

  LinearSIGMAmodule = new(rupturesRcpp:::Cost_LinearSIGMA, Y3, Xconst, TRUE, TRUE, 1e-6, TRUE)
  #warnOnce_ = TRUE -> keepWarning = FALSE (expected behavior when running eval() inside segmentation methods)
  expect_warning(LinearSIGMAmodule$eval(0, 100), "System is singular")

  #warnOnce_ = FALSE -> keepWarning = TRUE (expected behavior when running eval() outside segmentation methods)
  LinearSIGMAmodule = new(rupturesRcpp:::Cost_LinearSIGMA, Y3, Xconst, TRUE, TRUE, 1e-6, FALSE)
  expect_warning(LinearSIGMAmodule$eval(0, 100), "Singular system encountered")

})


test_that("Active binding `costFunc` (LinearSIGMA) works as intended for PELT/binSeg/Window", {

  set.seed(12345)
  costFuncObj = costFunc$new("LinearSIGMA")

  expect_equal(costFuncObj$intercept, TRUE)
  expect_equal(costFuncObj$addSmallDiag, TRUE)
  expect_equal(costFuncObj$epsilon, 1e-6)

  costFuncObj$addSmallDiag = FALSE
  costFuncObj$epsilon = 1e-4
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "intercept", "addSmallDiag", "epsilon")))
  expect_equal(args$addSmallDiag, FALSE)
  expect_equal(args$epsilon, 1e-4)

  #No covariates found -> force-fit with intercept-only, with a warning, for all three algorithms
  X_constantSeg = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))

  PELTObj = PELT$new(costFunc = costFunc$new("LinearSIGMA"))
  expect_warning(PELTObj$fit(X_constantSeg), "No `covariates` found!")
  expect_no_error(PELTObj$predict(0.1))

  binSegObj = binSeg$new(costFunc = costFunc$new("LinearSIGMA"))
  expect_warning(binSegObj$fit(X_constantSeg), "No `covariates` found!")
  expect_no_error(binSegObj$predict(0.1))

  windowObj = Window$new(costFunc = costFunc$new("LinearSIGMA"))
  expect_warning(windowObj$fit(X_constantSeg), "No `covariates` found!")
  expect_no_error(windowObj$predict(0.1))

})
