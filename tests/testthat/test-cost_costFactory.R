#costFactory R6 class

# ========================================================
#                  (R) Cost functions
# ========================================================

R_L1eval = function(X, start, end){
  Xe = X[(start+1):end, , drop = FALSE]
  return(sum(abs(sweep(Xe, 2, apply(Xe, 2, median), FUN = "-"))))
}

R_L2eval = function(X, start, end){
  Xe = X[(start+1):end, , drop = FALSE]
  return(sum(sweep(Xe, 2, colMeans(Xe), FUN = "-")^2))
}

R_SIGMAcov = function(X, start, end, epsilon = 1e-6){
  Xe = X[(start+1):end, , drop = FALSE]
  return(crossprod(sweep(Xe, 2, colMeans(Xe), FUN = "-"))/nrow(Xe) + diag(epsilon, ncol(Xe)))
}

# ========================================================
#                   Simulated datasets
# ========================================================

set.seed(12345)
tsMat = cbind(c(rnorm(50,0), rnorm(50,5,5)),
              c(rnorm(50,0), rnorm(50,-5)))
nr = nrow(tsMat)
X = cbind(rnorm(nr))
Y = cbind(X[,1]*rep(c(1,5), each = 50) + rnorm(nr))

nCases = 10
idx1 = sample.int(nr-10, nCases)
idx2 = integer(nCases)

for(i in 1:nCases){
  idx2[i] = sample((idx1[i]+5):nr, 1) #At least 5 observations per segment
}


test_that("Expect $eval() and $get_params() in costFactory give the correct results (L1, L2, SIGMA)", {

  L1Obj = costFactory$new(costFunc$new("L1"))
  L1Obj$fit(tsMat)
  L2Obj = costFactory$new(costFunc$new("L2"))
  L2Obj$fit(tsMat)
  SIGMAObj = costFactory$new(costFunc$new("SIGMA"))
  SIGMAObj$fit(tsMat)

  for(i in 1:nCases){
    Xe = tsMat[(idx1[i]+1):idx2[i], , drop = FALSE]
    gsCov = R_SIGMAcov(tsMat, idx1[i], idx2[i])

    expect_equal(L1Obj$eval(idx1[i],idx2[i]), R_L1eval(tsMat, idx1[i],idx2[i]))
    expect_equal(L1Obj$get_params(idx1[i],idx2[i])$median, apply(Xe, 2, median))

    expect_equal(L2Obj$eval(idx1[i],idx2[i]), R_L2eval(tsMat, idx1[i],idx2[i]))
    expect_equal(L2Obj$get_params(idx1[i],idx2[i])$mean, colMeans(Xe))

    expect_equal(SIGMAObj$eval(idx1[i],idx2[i]), nrow(Xe)*log(det(gsCov)))
    expect_equal(SIGMAObj$get_params(idx1[i],idx2[i])$cov, gsCov)
  }

})

test_that("Expect $eval() and $get_params() in costFactory give the correct results (VAR)", {

  VARObj = costFactory$new(costFunc$new("VAR", pVAR = 1))
  VARObj$fit(tsMat)

  for(i in 1:nCases){
    t = (idx1[i]+1):idx2[i] #idx1 >= 1, so every row has its lag
    fit = lm.fit(cbind(1, tsMat[t-1, , drop = FALSE]), tsMat[t, , drop = FALSE])

    expect_equal(VARObj$eval(idx1[i],idx2[i]), sum(fit$residuals^2))
    expect_equal(unname(VARObj$get_params(idx1[i],idx2[i])$coef), unname(fit$coefficients))
  }

  #Zero cost and NA coef if the segment has fewer observations than coefficients
  expect_equal(VARObj$eval(0,2), 0)
  expect_true(all(is.na(VARObj$get_params(0,2)$coef)))

})

test_that("Expect $eval() and $get_params() in costFactory give the correct results (LinearL2, LinearSIGMA)", {

  LinearL2Obj = costFactory$new(costFunc$new("LinearL2"))
  LinearL2Obj$fit(Y, X)
  LinearSIGMAObj = costFactory$new(costFunc$new("LinearSIGMA"))
  LinearSIGMAObj$fit(Y, X)

  for(i in 1:nCases){
    fit = lm(Y[(idx1[i]+1):idx2[i],] ~ X[(idx1[i]+1):idx2[i],])
    n = idx2[i] - idx1[i]
    gsCov = sum(fit$residuals^2)/n + 1e-6

    expect_equal(LinearL2Obj$eval(idx1[i],idx2[i]), sum(fit$residuals^2))
    expect_equal(c(LinearL2Obj$get_params(idx1[i],idx2[i])$coef), unname(coef(fit)))

    expect_equal(LinearSIGMAObj$eval(idx1[i],idx2[i]), n*log(gsCov))
    expect_equal(c(LinearSIGMAObj$get_params(idx1[i],idx2[i])$cov), gsCov)
  }

  #No covariates: force-fit with only an intercept
  LinearL2NoCovObj = costFactory$new(costFunc$new("LinearL2"))
  expect_warning(LinearL2NoCovObj$fit(Y),
                 regexp = "No `covariates` found! Force-fitting with only an intercept!")
  expect_equal(c(LinearL2NoCovObj$get_params(0,50)$coef), mean(Y[1:50,]))
  expect_equal(LinearL2NoCovObj$eval(0,50), R_L2eval(Y, 0, 50))

})

test_that("Expect $eval() and $get_params() in costFactory give the correct results (LinearL1, Custom)", {

  R_customL2params = function(segment, a, b) colMeans(segment)

  LinearL1Obj = costFactory$new(costFunc$new("LinearL1"))
  LinearL1Obj$fit(Y, X)
  LinearL1module = new(rupturesRcpp:::Cost_LinearL1, Y, X, TRUE, 1e-6, 1000L, FALSE)
  CustomObj = costFactory$new(costFunc$new("Custom", evalFun = function(segment, a, b) R_L2eval(segment, 0, nrow(segment)),
                                           paramFun = R_customL2params))
  CustomObj$fit(tsMat)

  for(i in 1:nCases){
    expect_equal(LinearL1Obj$eval(idx1[i],idx2[i]), LinearL1module$eval(idx1[i],idx2[i]))
    expect_equal(LinearL1Obj$get_params(idx1[i],idx2[i]), LinearL1module$get_params(idx1[i],idx2[i]))

    expect_equal(CustomObj$eval(idx1[i],idx2[i]), R_L2eval(tsMat, idx1[i],idx2[i]))
    expect_equal(CustomObj$get_params(idx1[i],idx2[i])$params,
                 colMeans(tsMat[(idx1[i]+1):idx2[i], , drop = FALSE]))
  }

  #No paramFun: empty list
  CustomObj2 = costFactory$new(costFunc$new("Custom", evalFun = function(segment, a, b) 0))
  CustomObj2$fit(tsMat)
  expect_equal(CustomObj2$get_params(0,nr), list())

})

test_that("Expect correct error messages in costFactory", {

  expect_error(costFactory$new("L2"),
               regexp = "`costFunc` must be a `R6` object of class `costFunc`")

  LinearL2MismatchObj = costFactory$new(costFunc$new("LinearL2"))
  expect_error(LinearL2MismatchObj$fit(Y, X[1:10, , drop = FALSE]),
               regexp = "Numbers of observations in `covariates` and `tsMat` do not match!")

  L2Obj = costFactory$new() #L2 by default
  L2Obj$fit(tsMat)
  expect_equal(L2Obj$eval(0,nr), R_L2eval(tsMat, 0, nr))

  expect_error(L2Obj$eval(0,nr+1),
               regexp = "out of bounds") #arma::mat indexing error

  #Give 0 if end - start = 0 or 1
  expect_equal(L2Obj$eval(0,0), 0)
  expect_equal(L2Obj$eval(0,1), 0)

  expect_error(L2Obj$get_params(0,0),
               regexp = "`start < end` must be true!")
  expect_error(L2Obj$get_params(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")

})

test_that("Expect $new()/$fit() separation and $costFunc active-binding behaviour in costFactory", {

  cf = costFactory$new(costFunc$new("L2"))

  #$eval()/$get_params() require $fit() first
  expect_error(cf$eval(0, 10), regexp = "must be run before")
  expect_error(cf$get_params(0, 10), regexp = "must be run before")

  #$fit() validates its data
  expect_error(cf$fit(), regexp = "No `tsMat` found! Please provide a `tsMat`!")
  expect_error(cf$fit("not a matrix"), regexp = "`tsMat` must be a numeric time series matrix!")

  naMat = tsMat
  naMat[1,1] = NA
  expect_error(cf$fit(naMat), regexp = "`tsMat` contains NAs!")

  cf$fit(tsMat)
  expect_equal(cf$eval(idx1[1], idx2[1]), R_L2eval(tsMat, idx1[1], idx2[1]))

  #Modifying `$costFunc` after fitting automatically re-fits with the new cost function
  expect_message(cf$costFunc <- costFunc$new("SIGMA"),
                 regexp = "`costFunc` has been updated. Re-fitting the model.")

  SIGMAObj = costFactory$new(costFunc$new("SIGMA"))
  SIGMAObj$fit(tsMat)
  expect_equal(cf$eval(idx1[1], idx2[1]), SIGMAObj$eval(idx1[1], idx2[1]))
  expect_identical(cf$costFunc$pass()$costFunc, "SIGMA")

  #Assigning a non-`costFunc` object errors and leaves the previous one in place
  expect_error(cf$costFunc <- "SIGMA",
               regexp = "`costFunc` must be a `R6` object of class `costFunc`")
  expect_identical(cf$costFunc$pass()$costFunc, "SIGMA")

})
