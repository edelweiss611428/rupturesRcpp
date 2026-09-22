#LinearL1 module (IRLS)

# ========================================================
#          (R) LinearL1 cost function, via IRLS
# ========================================================
# Mirrors Cost_LinearL1::fitColumnIRLS() in src/costs.cpp step-by-step (same
# initial OLS start, same 1e-6 weight floor, same convergence check), so
# results can be compared to the C++ implementation with tight tolerance.

R_fitColumnIRLS = function(y, Z, tol = 1e-6, maxIter = 50L){

  #`%*%`/solve() always return a matrix, even from vector inputs; as.vector()
  #keeps y/B/resid plain vectors throughout so elementwise ops below (which mix
  #an n-length vector with the n x J matrix Z) don't hit R's stricter
  #array-vs-array conformability check that a stray n x 1 matrix would trigger.
  y = as.vector(y)
  B = as.vector(solve(t(Z)%*%Z, t(Z)%*%y))
  resid = y - as.vector(Z%*%B)
  prevCost = sum(abs(resid))

  for(iter in seq_len(maxIter)){

    w = 1/pmax(abs(resid), 1e-6)
    sqrtw = sqrt(w)
    Zw = sqrtw*Z
    yw = sqrtw*y

    B = as.vector(solve(t(Zw)%*%Zw, t(Zw)%*%yw))
    resid = y - as.vector(Z%*%B)
    cost = sum(abs(resid))

    if(abs(prevCost - cost) <= tol*(1+prevCost)){
      prevCost = cost
      break
    }
    prevCost = cost
  }

  list(coef = as.vector(B), cost = prevCost)
}

R_LinearL1eval = function(Y, X, start, end, intercept = TRUE, tol = 1e-6, maxIter = 50L){

  Z = X[(start+1):end, , drop = FALSE]
  if(intercept){
    Z = cbind(1, Z)
  }
  Ysub = Y[(start+1):end, , drop = FALSE]

  total = 0
  for(j in seq_len(ncol(Ysub))){
    total = total + R_fitColumnIRLS(Ysub[,j], Z, tol, maxIter)$cost
  }
  total
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
      break #Ensure there are enough observations to fit the regression
    }
  }
}


test_that("Expect C++ .eval() method in LinearL1 cost module gives the correct results", {

  LinearL1module = new(rupturesRcpp:::Cost_LinearL1, Y, X, TRUE, 1e-6, 50L, TRUE) # (intercept, tol, maxIter, warnOnce)

  for(i in 1:nCases){
    gs = R_LinearL1eval(Y, X, idx1[i], idx2[i])
    expect_equal(LinearL1module$eval(idx1[i], idx2[i]), gs, tolerance = 1e-6)
  }

  expect_error(LinearL1module$eval(0,nr+1),
               regexp = "out of bounds") #arma::mat indexing error

  expect_error(LinearL1module$eval(-1,nr),
               regexp = "out of bounds") #arma::mat indexing error

  #Give 0 if end - start = 0, 1, or fewer than J observations
  expect_equal(LinearL1module$eval(0,0), 0)
  expect_equal(LinearL1module$eval(0,1), 0)
  expect_equal(LinearL1module$eval(0,J-1), 0)

  expect_no_error(LinearL1module$resetWarning(FALSE)) #Generally does nothing here

})


test_that("Expect $get_params() method in LinearL1 cost module gives the correct results", {

  LinearL1module = new(rupturesRcpp:::Cost_LinearL1, Y, X, TRUE, 1e-6, 50L, TRUE)

  for(i in 1:nCases){

    Zseg = cbind(1, X[(idx1[i]+1):idx2[i],])
    gsCoef = cbind(R_fitColumnIRLS(Y[(idx1[i]+1):idx2[i],1], Zseg)$coef,
                    R_fitColumnIRLS(Y[(idx1[i]+1):idx2[i],2], Zseg)$coef)

    params = LinearL1module$get_params(idx1[i],idx2[i])
    expect_equal(unname(params$coef), unname(gsCoef), tolerance = 1e-6)
  }

  #NA coef if segment has fewer usable observations than coefficients
  params = LinearL1module$get_params(0, J-1)
  expect_true(all(is.na(params$coef)))

})


test_that("Expect $eval() method in PELT_LinearL1 gives the correct results/error message", {

  set.seed(12345)

  PELTObj = PELT$new(costFunc = costFunc$new("LinearL1"))
  PELTObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearL1eval(Y, X, idx1[i], idx2[i])
    expect_equal(PELTObj$eval(idx1[i],idx2[i]), gs, tolerance = 1e-6)
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


test_that("Expect .eval() method in C++ PELT_LinearL1 class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, double, int, int, int>()
  PELTCppObj = new(PELTCpp_LinearL1, Y, X, TRUE, 1e-6, 50L, 1L, 1L)

  for(i in 1:nCases){
    gs = R_LinearL1eval(Y, X, idx1[i], idx2[i])
    expect_equal(PELTCppObj$eval(idx1[i],idx2[i]), gs, tolerance = 1e-6)
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


test_that("Expect $eval() method in binSeg_LinearL1 gives the correct results/error message", {

  binSegObj = binSeg$new(costFunc = costFunc$new("LinearL1"))
  binSegObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearL1eval(Y, X, idx1[i], idx2[i])
    expect_equal(binSegObj$eval(idx1[i],idx2[i]), gs, tolerance = 1e-6)
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


test_that("Expect .eval() method in C++ binSeg_LinearL1 class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, double, int, int, int>()
  binSegCppObj = new(binSegCpp_LinearL1, Y, X, TRUE, 1e-6, 50L, 1L, 1L)

  for(i in 1:nCases){
    gs = R_LinearL1eval(Y, X, idx1[i], idx2[i])
    expect_equal(binSegCppObj$eval(idx1[i],idx2[i]), gs, tolerance = 1e-6)
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


test_that("Expect $eval() method in window_LinearL1 gives the correct results/error message", {

  windowObj = Window$new(costFunc = costFunc$new("LinearL1"))
  windowObj$fit(Y, X)

  for(i in 1:nCases){
    gs = R_LinearL1eval(Y, X, idx1[i], idx2[i])
    expect_equal(windowObj$eval(idx1[i],idx2[i]), gs, tolerance = 1e-6)
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


test_that("Expect .eval() method in C++ window_LinearL1 class gives the correct results/error message", {

  #.constructor<arma::mat, arma::mat, bool, double, int, int, int, int>()
  windowCppObj = new(windowCpp_LinearL1, Y, X, TRUE, 1e-6, 50L, 1L, 1L, 10L)

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
    gs = R_LinearL1eval(Y, X, idx1b[i], idx2b[i])
    expect_equal(windowCppObj$eval(idx1b[i],idx2b[i]), gs, tolerance = 1e-6)
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


test_that("Expect correct error/warning messages when initialising C++ Cost_LinearL1 module", {

  expect_error(new(rupturesRcpp:::Cost_LinearL1, Y, as.matrix(X[-1]), TRUE, 1e-6, 50L, TRUE),
               "Number of observations in response and covariate matrices must match!")

  X2 = matrix(1:5, nrow = 1)
  Y2 = matrix(1, nrow = 1)

  expect_error(new(rupturesRcpp:::Cost_LinearL1, Y2, X2, TRUE, 1e-6, 50L, TRUE),
               "The full dataset contains not enough observations to fit a linear regression model!")

  expect_error(new(rupturesRcpp:::Cost_LinearL1, Y, X, TRUE, 0, 50L, TRUE),
               "`tol` must be a single positive value!")

  expect_error(new(rupturesRcpp:::Cost_LinearL1, Y, X, TRUE, -1, 50L, TRUE),
               "`tol` must be a single positive value!")

  expect_error(new(rupturesRcpp:::Cost_LinearL1, Y, X, TRUE, 1e-6, 0L, TRUE),
               "`maxIter` must be at least 1!")

  #Singular design matrix (X duplicates the intercept column) -> initial OLS solve is singular
  Xconst = matrix(rep(1,100))
  Y3 = matrix(rnorm(100))

  LinearL1module = new(rupturesRcpp:::Cost_LinearL1, Y3, Xconst, TRUE, 1e-6, 50L, TRUE)
  #warnOnce_ = TRUE -> keepWarning = FALSE (expected behavior when running eval() inside segmentation methods)
  expect_warning(LinearL1module$eval(0, 100), "System is singular")

  #warnOnce_ = FALSE -> keepWarning = TRUE (expected behavior when running eval() outside segmentation methods)
  LinearL1module = new(rupturesRcpp:::Cost_LinearL1, Y3, Xconst, TRUE, 1e-6, 50L, FALSE)
  expect_warning(LinearL1module$eval(0, 100), "Singular system encountered")

  #maxIter = 1L is (essentially) never enough to converge -> non-convergence warning
  Y4 = matrix(rnorm(100))
  X4 = matrix(rnorm(100))

  LinearL1module = new(rupturesRcpp:::Cost_LinearL1, Y4, X4, TRUE, 1e-12, 1L, TRUE)
  expect_warning(LinearL1module$eval(0, 100), "IRLS did not converge")

  LinearL1module = new(rupturesRcpp:::Cost_LinearL1, Y4, X4, TRUE, 1e-12, 1L, FALSE)
  expect_warning(LinearL1module$eval(0, 100), "IRLS did not converge")

})


test_that("Active binding `costFunc` (LinearL1) works as intended for PELT/binSeg/Window", {

  set.seed(12345)
  costFuncObj = costFunc$new("LinearL1")

  expect_equal(costFuncObj$intercept, TRUE)
  expect_equal(costFuncObj$tol, 1e-6)
  expect_equal(costFuncObj$maxIter, 50L)

  costFuncObj$tol = 1e-4
  costFuncObj$maxIter = 10L
  args = costFuncObj$pass()
  expect_true(setequal(names(args), c("costFunc", "intercept", "tol", "maxIter")))
  expect_equal(args$tol, 1e-4)
  expect_equal(args$maxIter, 10L)

  #Invalid `tol`/`maxIter`
  expect_error(costFuncObj$tol <- 0, "single positive value")
  expect_error(costFuncObj$tol <- -1, "single positive value")
  expect_error(costFuncObj$maxIter <- 0, "single positive integer")

  #No covariates found -> force-fit with intercept-only, with a warning, for all three algorithms
  X_constantSeg = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))

  PELTObj = PELT$new(costFunc = costFunc$new("LinearL1"))
  expect_warning(PELTObj$fit(X_constantSeg), "No `covariates` found!")
  expect_no_error(PELTObj$predict(0.1))

  binSegObj = binSeg$new(costFunc = costFunc$new("LinearL1"))
  expect_warning(binSegObj$fit(X_constantSeg), "No `covariates` found!")
  expect_no_error(binSegObj$predict(0.1))

  windowObj = Window$new(costFunc = costFunc$new("LinearL1"))
  expect_warning(windowObj$fit(X_constantSeg), "No `covariates` found!")
  expect_no_error(windowObj$predict(0.1))

})
