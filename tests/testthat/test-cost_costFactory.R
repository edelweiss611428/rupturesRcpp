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

  L1Obj = costFactory$new(costFunc$new("L1"), tsMat)
  L2Obj = costFactory$new(costFunc$new("L2"), tsMat)
  SIGMAObj = costFactory$new(costFunc$new("SIGMA"), tsMat)

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

  VARObj = costFactory$new(costFunc$new("VAR", pVAR = 1), tsMat)

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

  LinearL2Obj = costFactory$new(costFunc$new("LinearL2"), Y, X)
  LinearSIGMAObj = costFactory$new(costFunc$new("LinearSIGMA"), Y, X)

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
  expect_warning(LinearL2Obj <- costFactory$new(costFunc$new("LinearL2"), Y),
                 regexp = "No `covariates` found! Force-fitting with only an intercept!")
  expect_equal(c(LinearL2Obj$get_params(0,50)$coef), mean(Y[1:50,]))
  expect_equal(LinearL2Obj$eval(0,50), R_L2eval(Y, 0, 50))

})

test_that("Expect $eval() and $get_params() in costFactory give the correct results (LinearL1, Custom)", {

  R_customL2params = function(segment, a, b) colMeans(segment)

  LinearL1Obj = costFactory$new(costFunc$new("LinearL1"), Y, X)
  LinearL1module = new(rupturesRcpp:::Cost_LinearL1, Y, X, TRUE, 1e-6, 1000L, FALSE)
  CustomObj = costFactory$new(costFunc$new("Custom", evalFun = function(segment, a, b) R_L2eval(segment, 0, nrow(segment)),
                                           paramFun = R_customL2params), tsMat)

  for(i in 1:nCases){
    expect_equal(LinearL1Obj$eval(idx1[i],idx2[i]), LinearL1module$eval(idx1[i],idx2[i]))
    expect_equal(LinearL1Obj$get_params(idx1[i],idx2[i]), LinearL1module$get_params(idx1[i],idx2[i]))

    expect_equal(CustomObj$eval(idx1[i],idx2[i]), R_L2eval(tsMat, idx1[i],idx2[i]))
    expect_equal(CustomObj$get_params(idx1[i],idx2[i])$params,
                 colMeans(tsMat[(idx1[i]+1):idx2[i], , drop = FALSE]))
  }

  #No paramFun: empty list
  CustomObj = costFactory$new(costFunc$new("Custom", evalFun = function(segment, a, b) 0), tsMat)
  expect_equal(CustomObj$get_params(0,nr), list())

})

test_that("Expect correct error messages in costFactory", {

  expect_error(costFactory$new("L2", tsMat),
               regexp = "`costFunc` must be a `R6` object of class `costFunc`")
  expect_error(costFactory$new(costFunc$new()),
               regexp = "`tsMat` must be provided!")
  expect_error(costFactory$new(costFunc$new(), "a"),
               regexp = "`tsMat` must be a numeric time series matrix!")
  expect_error(costFactory$new(costFunc$new(), cbind(c(1, NA))),
               regexp = "`tsMat` contains NAs!")
  expect_error(costFactory$new(costFunc$new("LinearL2"), Y, X[1:10, , drop = FALSE]),
               regexp = "Numbers of observations in `covariates` and `tsMat` do not match!")

  L2Obj = costFactory$new(tsMat = tsMat) #L2 by default
  expect_equal(L2Obj$eval(0,nr), R_L2eval(tsMat, 0, nr))

  expect_error(L2Obj$eval(-1,nr),
               regexp = "`0 <= start < nSamples` must be true!")
  expect_error(L2Obj$eval(0,nr+1),
               regexp = "`0 < end <= nSamples` must be true!")
  expect_error(L2Obj$eval(0,0),
               regexp = "a must be smaller than b")
  expect_error(L2Obj$get_params(0,0),
               regexp = "a must be smaller than b")

})
