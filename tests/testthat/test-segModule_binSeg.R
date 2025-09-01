#binSeg module
set.seed(12345)


test_that("binSeg_L1 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("L1")
  binSegObj = binSeg$new(costFunc = costFuncObj)
  binSegObj$fit(tsMat)

  expect_equal( binSegObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("binSeg_L2 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("L2")
  binSegObj = binSeg$new(costFunc = costFuncObj)
  binSegObj$fit(tsMat)

  expect_equal( binSegObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("binSeg_SIGMA works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("SIGMA")
  binSegObj = binSeg$new(costFunc = costFuncObj)
  binSegObj$fit(tsMat)

  expect_equal( binSegObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("binSeg_VAR works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("VAR")
  binSegObj = binSeg$new(costFunc = costFuncObj)

  expect_warning(binSegObj$fit(tsMat), "Some systems seem singular!") #Warning once feature
  expect_equal(binSegObj$predict(pen = 0.1), seq(50,150,50))
  expect_warning(expect_false(binSegObj$eval(0,51) == 0), "seems singular")
  expect_warning(expect_true(all.equal(binSegObj$eval(0,50), 0)), "seems singular")

})


test_that("binSeg_LinearL2 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("LinearL2")
  binSegObj = binSeg$new(costFunc = costFuncObj)

  #when no covariate matrix is provided
  expect_warning(binSegObj$fit(tsMat), "an intercept")
  expect_equal(binSegObj$predict(pen = 0.1), seq(50,150,50))
  expect_false(binSegObj$eval(0,51) == 0)
  expect_true(all.equal(binSegObj$eval(0,50), 0))
})


test_that("Active binding `minSize` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a single non-negative numeric/integer
  expect_error(binSeg$new(costFunc = costFuncObj, minSize = "a"))
  expect_error(binSeg$new(costFunc = costFuncObj, minSize = NULL))
  expect_error(binSeg$new(costFunc = costFuncObj, minSize = 1:2))
  expect_error(binSeg$new(costFunc = costFuncObj, minSize = 0))

  #Getter
  binSegObj = binSeg$new(costFunc = costFuncObj)
  expect_equal(binSegObj$minSize, 1L)

  #Setter
  binSegObj$minSize = 5L
  expect_equal(binSegObj$minSize, 5L)
  expect_error(binSegObj$minSize <- "a")
  expect_error(binSegObj$minSize <- NULL)
  expect_error(binSegObj$minSize <- 1:2)
  expect_error(binSegObj$minSize <- 0)

  #Modifying `minSize` triggers refitting if fitted
  binSegObj$minSize = 1L
  binSegObj$fit(tsMat)
  ms1Seg = binSegObj$predict(0)
  expect_equal(ms1Seg, seq(1, 100, 1))

  binSegObj$minSize = 50L
  ms50Seg = binSegObj$predict(0)
  expect_equal(ms50Seg, c(50, 100))

})


test_that("Active binding `jump` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a single non-negative numeric/integer
  expect_error(binSeg$new(costFunc = costFuncObj, jump = "a"))
  expect_error(binSeg$new(costFunc = costFuncObj, jump  = NULL))
  expect_error(binSeg$new(costFunc = costFuncObj, jump  = 1:2))
  expect_error(binSeg$new(costFunc = costFuncObj, jump  = 0))
  expect_no_error(binSeg$new())

  #Getter
  binSegObj = binSeg$new(costFunc = costFuncObj)
  expect_equal(binSegObj$jump, 1L)

  #Setter
  binSegObj$jump  = 5L
  expect_equal(binSegObj$jump, 5L)
  expect_error(binSegObj$jump <- "a")
  expect_error(binSegObj$jump <- NULL)
  expect_error(binSegObj$jump <- 1:2)
  expect_error(binSegObj$jump <- 0)

  #Modifying `jump` triggers refitting if fitted
  binSegObj$jump  = 1L
  binSegObj$fit(tsMat)
  j1Seg = binSegObj$predict(0)
  expect_equal(j1Seg, seq(1, 100, 1))

  binSegObj$jump  = 5L
  j5Seg = binSegObj$predict(0)
  expect_equal(j5Seg, seq(5, 100, 5))

})


test_that("Active binding `tsMat` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  set.seed(123)
  tsMat2 = matrix(c(rnorm(50,0), rnorm(50,5)))
  tsNA = tsMat
  tsNA[1] = NA
  binSegObj = binSeg$new()

  #Wrong input: `tsMat` must be a numeric matrix without any NA
  expect_error(binSegObj$fit(tsNA))
  expect_error(binSegObj$fit(c("a", "b")))
  expect_error(binSegObj$fit(NULL))
  expect_error(binSegObj$fit(as.vector(tsMat)))
  expect_no_error(binSegObj$fit(tsMat))

  #Getter
  binSegObj$fit(tsMat)
  expect_equal(binSegObj$tsMat, tsMat)

  #Setter
  expect_error(binSegObj$tsMat <- tsNA)
  expect_error(binSegObj$tsMat <- c("a", "b"))
  expect_error(binSegObj$tsMat <- NULL)
  expect_error(binSegObj$tsMat <- as.vector(tsMat))

  binSegObj$tsMat = tsMat2
  expect_equal(binSegObj$tsMat, tsMat2)

  #Modify `tsMat` triggers refitting if fitted

  binSegObj$fit(tsMat)
  tsMat1Eval = binSegObj$eval(0,100)
  expect_equal(tsMat1Eval, sum((tsMat - mean(tsMat))^2))

  binSegObj$fit(tsMat2)
  tsMat2Eval = binSegObj$eval(0,100)
  expect_equal(tsMat2Eval, sum((tsMat2 - mean(tsMat2))^2))

})




test_that("Active binding `covariates` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  set.seed(123)
  covariateMat = as.matrix(tsMat/2 + rnorm(100))
  covariateNA = covariateMat
  covariateNA[1] = NA
  costFuncObj = costFunc$new("LinearL2")
  binSegObj = binSeg$new(costFunc = costFuncObj)

  #Wrong input: `covariate` must be a numeric matrix without any NA
  expect_error(binSegObj$fit(tsMat, as.matrix(covariateMat[-1]))) #nrows do not match
  expect_error(binSegObj$fit(tsMat, covariateNA))
  expect_error(binSegObj$fit(tsMat, c("a", "b")))
  expect_error(binSegObj$fit(tsMat, as.vector(covariateMat)))
  expect_no_error(binSegObj$fit(tsMat, covariateMat))
  expect_no_error(binSegObj$fit(tsMat, NULL))

  #Getter
  binSegObj$fit(tsMat, covariateMat)
  expect_equal(binSegObj$covariates, covariateMat)

  #Setter
  expect_error(binSegObj$covariates <- as.matrix(covariateMat[-1])) #nrows do not match
  expect_error(binSegObj$covariates <- covariateNA)
  expect_error(binSegObj$covariates <- c("a", "b"))
  expect_error(binSegObj$covariates <- NULL)
  expect_error(binSegObj$covariates <- as.vector(covariateMat))

  #Modify `tsMat` triggers refitting if fitted

  set.seed(100)
  covariateMat2 = as.matrix(tsMat/2 + rnorm(100))

  binSegObj$fit(tsMat, covariateMat)
  cM1err = binSegObj$eval(0,100)
  binSegObj$covariates = covariateMat2
  cM2err = binSegObj$eval(0,100)

  expect_equal(cM1err, sum(lm(tsMat~covariateMat)$residuals^2))
  expect_equal(cM2err, sum(lm(tsMat~covariateMat2)$residuals^2))

})

test_that("Active binding `costFunc` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a R6 object of class `costFunc`
  expect_error(binSeg$new(costFunc = 1L))
  expect_error(binSeg$new(costFunc = list(costFunc = "L2")))
  expect_error(binSeg$new(costFunc = NULL))
  expect_no_error(binSeg$new())

  #Getter
  binSegObj = binSeg$new(costFunc = costFuncObj)
  expect_equal(binSegObj$costFunc$pass()$costFunc, "L2")

  #Setter
  #Wrong input: Must be a R6 object of class `costFunc`
  expect_error(binSegObj$costFunc <- list(costFunc = "L2"))
  expect_error(binSegObj$costFunc <- NULL)
  expect_error(binSegObj$costFunc <- 1L)
  expect_no_error(binSegObj$costFunc <- costFunc$new("VAR"))

  #Modifying `costFunc` triggers refitting if fitted

  binSegObj = binSeg$new() #L2
  binSegObj$fit(tsMat)
  expect_equal(binSegObj$eval(0, 100), sum((tsMat - mean(tsMat))^2))

  binSegObj$costFunc = costFunc$new("L1")
  expect_equal(binSegObj$eval(0, 100), sum(abs(tsMat - median(tsMat))))

  binSegObj$costFunc = costFunc$new("SIGMA")
  expect_equal(binSegObj$eval(0, 100), 100*log(det(var(tsMat)*99/100+10^-6)))

})


test_that("Test that `describe()` method works properly)", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")
  binSegObj = binSeg$new(costFunc = costFuncObj)
  binSegObj$fit(tsMat)

  #Wrong input
  expect_error(binSegObj$describe(NULL))
  expect_error(binSegObj$describe(1:2))
  expect_error(binSegObj$describe("a"))

  #`describe()` returns the expected outputs
  expect_no_error(binSegObj$describe(T))
  expect_no_error(binSegObj$describe(F))
  expect_no_error(binSegObj$describe())

  expect_equal(binSegObj$describe(T), binSegObj$describe(F))
  expect_equal(binSegObj$describe(), list(minSize = 1L, jump = 1L,
                                          costFunc = costFuncObj, fitted = T,
                                          tsMat = tsMat, covariates = NULL,
                                          n = 100, p = 1))

  #Modify non-`costFunc` active bindings

  set.seed(123)
  tsMat2 = cbind(x1 = c(rnorm(50,0), rnorm(50,5)),
                 x2 = c(rnorm(50,0), rnorm(50,5)))
  binSegObj$tsMat = tsMat2
  binSegObj$minSize = 2L
  binSegObj$jump = 2L

  expect_equal(binSegObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFuncObj, fitted = T,
                                          tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2))

  #Modify `costFunc``
  costFuncObj  = costFunc$new(costFunc = "VAR", pVAR = 1L)
  binSegObj$costFunc = costFuncObj

  expect_equal(binSegObj$describe(T), binSegObj$describe(F))
  expect_equal(binSegObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFuncObj,
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2, pVAR = costFuncObj$pVAR))

  costFuncObj = costFunc$new(costFunc = "SIGMA", pVAR = 1L)
  binSegObj$costFunc = costFuncObj

  expect_equal(binSegObj$describe(T), binSegObj$describe(F))
  expect_equal(binSegObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFunc$new(costFunc = "SIGMA"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2,
                                          addSmallDiag = costFuncObj$addSmallDiag,
                                          epsilon = costFuncObj$epsilon))

  expect_true(is.null(binSegObj$describe()$pVAR))


  costFuncObj = costFunc$new(costFunc = "L1", pVAR = 1L)
  binSegObj$costFunc = costFuncObj

  expect_equal(binSegObj$describe(T), binSegObj$describe(F))
  expect_equal(binSegObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFunc$new(costFunc = "L1"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2))

  expect_true(is.null(binSegObj$describe()$pVAR))
  expect_true(is.null(binSegObj$describe()$epsilon))
  expect_true(is.null(binSegObj$describe()$addSmallDiag))


  costFuncObj = costFunc$new(costFunc = "LinearL2", pVAR = 1L)
  expect_warning(binSegObj$costFunc <- costFuncObj)

  expect_equal(binSegObj$describe(T), binSegObj$describe(F))
  expect_equal(binSegObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFunc$new(costFunc = "LinearL2"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2, intercept = TRUE))

  expect_true(is.null(binSegObj$describe()$pVAR))

})



test_that("Error handling for `eval()` works properly", {

  #Cost-specific tests are in other files

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  binSegObj = binSeg$new()
  expect_error(binSegObj$eval(0, 10)) #Not fitted

  binSegObj$fit(tsMat)
  expect_error(binSegObj$eval(NULL, 10))
  expect_error(binSegObj$eval(0, NULL))

})

test_that("Test that error handling for `predict()` works properly", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  binSegObj = binSeg$new()
  expect_error(binSegObj$predict()) #Not fitted

  binSegObj$fit(tsMat)
  expect_error(binSegObj$predict(NULL))
  expect_error(binSegObj$predict(-1))
  expect_error(binSegObj$predict(c(1:2)))
  expect_error(binSegObj$predict("a"))
  expect_error(binSegObj$predict(T))

})




test_that("Error handling for `plot()` works properly", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  binSegObj = binSeg$new()

  #Without providing endPts
  binSegObj$fit(tsMat)
  expect_error(binSegObj$plot()) #No tmpEndPts available
  binSegObj$predict(25)
  expect_no_error(binSegObj$plot())

  ## Invalid `main`
  expect_error(binSegObj$plot(main = c(1:3)))
  expect_error(binSegObj$plot(main = T))
  expect_error(binSegObj$plot(main = NA))
  expect_error(binSegObj$plot(main = NULL))
  expect_error(binSegObj$plot(main = c("a", "b")))

  ## Invalid `xlab`
  expect_error(binSegObj$plot(xlab = c(1:3)))
  expect_error(binSegObj$plot(xlab = T))
  expect_error(binSegObj$plot(xlab = NA))
  expect_error(binSegObj$plot(xlab = NULL))
  expect_error(binSegObj$plot(xlab = c("a", "b")))

  ## Invalid `d`
  expect_error(binSegObj$plot(d = 2)) #d > p
  expect_error(binSegObj$plot(d = "a"))
  expect_error(binSegObj$plot(d = T))
  expect_error(binSegObj$plot(d = NULL))
  expect_error(binSegObj$plot(d = NA))

  ## Invalid `dimNames`
  expect_error(binSegObj$plot(dimNames = c("X1", "X2")))
  expect_error(binSegObj$plot(dimNames = T))
  expect_error(binSegObj$plot(dimNames = NA))
  expect_error(binSegObj$plot(dimNames = NULL))
  expect_error(binSegObj$plot(main = 123))

  #endPts provided
  expect_error(binSegObj$plot(endPts = "a"))
  expect_error(binSegObj$plot(endPts = NA))
  expect_error(binSegObj$plot(endPts = NULL))
  expect_error(binSegObj$plot(endPts = 1:99)) #endPts not include n
  expect_error(binSegObj$plot(endPts = 0:99)) ##min endPts < 1
  expect_error(binSegObj$plot(endPts = 0:101)) ##max endPts > n
  expect_error(binSegObj$plot(endPts = c(50,50,100))) ##duplicated endpts

})


test_that("Error handling for C++ module binSeg_L2 works as intended", {
  #constructor: const arma::mat& tsMat, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::binSegCpp_L2, tsMat, 0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_L2, tsMat, 1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_L2, tsMat, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::binSegCpp_L2, tsMat2, 10, 3)) #len = minLen here

})


test_that("Error handling for C++ module binSeg_L1 works as intended", {
  #constructor: const arma::mat& tsMat, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::binSegCpp_L1_cwMed, tsMat, 0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_L1_cwMed, tsMat, 1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_L1_cwMed, tsMat, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::binSegCpp_L1_cwMed, tsMat2, 10, 3)) #len = minLen here

})


test_that("Error handling for C++ module binSeg_VAR works as intended", {
  #constructor: const arma::mat& tsMat, int pVAR, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  #pVAR = 1
  expect_error(new(rupturesRcpp:::binSegCpp_VAR, tsMat, 1, 0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_VAR, tsMat, 1, 1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_VAR, tsMat, 1, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::binSegCpp_VAR, tsMat2, 1, 10, 3)) #len = minLen here

  set.seed(123)
  tsMat3 = cbind(x1 = rnorm(5), x2 = rnorm(5), x3 = rnorm(5))
  expect_error(new(rupturesRcpp:::binSegCpp_VAR, tsMat3, 2, 1, 1)) #segment too short for fitting VAR(2)
})


test_that("Error handling for C++ module binSeg_SIGMA works as intended", {
  #constructor: const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::binSegCpp_SIGMA, tsMat, T, 10^-6, 0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_SIGMA, tsMat, T, 10^-6, 1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_SIGMA, tsMat, T, 10^-6, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::binSegCpp_SIGMA, tsMat2, T, 10^-6, 10, 3)) #len = minLen here

})


test_that("Error handling for C++ module binSeg_LinearL2 works as intended", {
  #constructor: const arma::mat& tsMat,  const arma::mat& covariates, bool intercept_, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))
  covariateMat = as.matrix(rep(1,23))
  covariateMat2 = as.matrix(rep(1,24))

  expect_error(new(rupturesRcpp:::binSegCpp_LinearL2, tsMat, covariateMat, T,  0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_LinearL2, tsMat, covariateMat, T,  1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::binSegCpp_LinearL2, tsMat, covariateMat, T, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::binSegCpp_LinearL2, tsMat2, covariateMat2, T, 10, 3)) #len = minLen here

  set.seed(1234)
  tsMat3 = as.matrix(rnorm(5))
  covariateMat3 =cbind(rnorm(5),rnorm(5),rnorm(5),rnorm(5),rnorm(5))
  expect_error(new(rupturesRcpp:::binSegCpp_LinearL2, tsMat3, covariateMat3, T, 1, 1)) #too short for fitting linear regression model
})

test_that("Compare C++ module binSeg_L2 to tdhock/binsegRcpp", {

  skip_if_not_installed("binsegRcpp")
  set.seed(12345)

  for(i in 1:10){

    tsMat = matrix(rnorm(250))
    tsVec = as.vector(tsMat)

    binSegCppModule = new(rupturesRcpp:::binSegCpp_L2, tsMat, 1L, 1L)
    binSegCppModule$fit()
    edelweiss = binSegCppModule$predict(0)
    tdhock = binsegRcpp::binseg_normal(tsVec)$splits$end[-1]

    expect_equal(edelweiss, tdhock)

  }

})



test_that("Some additional tests", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  binSegObj = binSeg$new()
  binSegObj$fit(tsMat)
  expect_equal(binSegObj$predict(999999), 100) #Too large `pen` -> only return `n`


})




