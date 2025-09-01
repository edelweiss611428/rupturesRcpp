#PELT module

X_constantSeg = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))


test_that("PELT_L2 modules in `rupturesRcpp` and `ruptures` give the same results : Test 1", {

  set.seed(59)
  X_noBkp = matrix(rnorm(250))

  skip_if_not_installed("reticulate")
  reticulate = asNamespace("reticulate")

  if (!reticulate::py_module_available("numpy") ||
      !reticulate::py_module_available("ruptures")) {
    skip("Required Python modules not available")
  }

  ruptures = reticulate::import("ruptures")
  np = reticulate::import("numpy")
  npX_noBkp = np$array(X_noBkp)

  #PyPeltL2(min_size = 1L, jump = 1L)
  PyPelt = ruptures$Pelt(model = "l2", min_size = 1L, jump = 1L)
  PyPelt$fit(npX_noBkp)

  #RPeltL2(min_size = 1L, jump = 1L)
  RPelt = PELT$new(minSize = 1L, jump = 1L, costFunc = costFunc$new("L2"))
  RPelt$fit(X_noBkp)

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 3L, jump = 3L)
  PyPelt$min_size = 3L
  PyPelt$jump = 3L

  #RPeltL2(min_size = 3L, jump = 3L)
  RPelt$minSize = 3L
  RPelt$jump = 3L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 5L, jump = 5L)
  PyPelt$min_size = 5L
  PyPelt$jump = 5L

  #RPeltL2(min_size = 5L, jump = 5L)
  RPelt$minSize = 5L
  RPelt$jump = 5L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 5L, jump = 10L)
  PyPelt$min_size = 5L
  PyPelt$jump = 10L

  #RPeltL2(min_size = 5L, jump = 10L)
  RPelt$minSize = 5L
  RPelt$jump = 10L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 10L, jump = 1L)
  PyPelt$min_size = 10L
  PyPelt$jump = 1L

  #RPeltL2(min_size = 10L, jump = 1L)
  RPelt$minSize = 10L
  RPelt$jump = 1L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 10L, jump = 2L)
  PyPelt$min_size = 10L
  PyPelt$jump = 2L

  #RPeltL2(min_size = 10L, jump = 2L)
  RPelt$minSize = 10L
  RPelt$jump = 2L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }


})



test_that("PELT_L1 works for constant segments", {

  costFuncObj = costFunc$new("L1")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("PELT_L2 modules in `rupturesRcpp` and `ruptures` give the same results : Test 2", {

  set.seed(123)
  X_noBkp = matrix(rnorm(250))

  skip_if_not_installed("reticulate")
  reticulate = asNamespace("reticulate")

  if (!reticulate::py_module_available("numpy") ||
      !reticulate::py_module_available("ruptures")) {
    skip("Required Python modules not available")
  }

  ruptures = reticulate::import("ruptures")
  np = reticulate::import("numpy")
  npX_noBkp = np$array(X_noBkp)

  #PyPeltL2(min_size = 1L, jump = 1L)
  PyPelt = ruptures$Pelt(model = "l2", min_size = 1L, jump = 1L)
  PyPelt$fit(npX_noBkp)

  #RPeltL2(min_size = 1L, jump = 1L)
  RPelt = PELT$new(minSize = 1L, jump = 1L, costFunc = costFunc$new("L2"))
  RPelt$fit(X_noBkp)

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 3L, jump = 3L)
  PyPelt$min_size = 3L
  PyPelt$jump = 3L

  #RPeltL2(min_size = 3L, jump = 3L)
  RPelt$minSize = 3L
  RPelt$jump = 3L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 5L, jump = 5L)
  PyPelt$min_size = 5L
  PyPelt$jump = 5L

  #RPeltL2(min_size = 5L, jump = 5L)
  RPelt$minSize = 5L
  RPelt$jump = 5L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 5L, jump = 10L)
  PyPelt$min_size = 5L
  PyPelt$jump = 10L

  #RPeltL2(min_size = 5L, jump = 10L)
  RPelt$minSize = 5L
  RPelt$jump = 10L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 10L, jump = 1L)
  PyPelt$min_size = 10L
  PyPelt$jump = 1L

  #RPeltL2(min_size = 10L, jump = 1L)
  RPelt$minSize = 10L
  RPelt$jump = 1L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }

  #PyPeltL2(min_size = 10L, jump = 2L)
  PyPelt$min_size = 10L
  PyPelt$jump = 2L

  #RPeltL2(min_size = 10L, jump = 2L)
  RPelt$minSize = 10L
  RPelt$jump = 2L

  for(i in 1:10){

    pen = runif(1,0,1)
    PySol = PyPelt$predict(pen)
    RSol = RPelt$predict(pen)

    expect_true(all.equal(PySol, RSol))

  }


})



test_that("PELT_L1 works for constant segments", {

  costFuncObj = costFunc$new("L1")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("PELT_L2 works for constant segments", {

  costFuncObj = costFunc$new("L2")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("PELT_SIGMA works for constant segments", {

  costFuncObj = costFunc$new("SIGMA")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

})


test_that("PELT_VAR works for constant segments", {

  costFuncObj = costFunc$new("VAR")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(X_constantSeg)

  expect_warning(expect_equal(PeltObj$predict(pen = 0.1), seq(50,150,50)),
                 "Some systems seem singular!")
  expect_warning(expect_false(PeltObj$eval(0,51) == 0), "seems singular")
  expect_warning(expect_true(all.equal(PeltObj$eval(0,50), 0)), "seems singular")

})

test_that("PELT_LinearL2 works for constant segments", {

  costFuncObj = costFunc$new("LinearL2")
  PELTObj = PELT$new(costFunc = costFuncObj)
  expect_warning(PELTObj$fit(X_constantSeg), "No `covariates` found!")

  expect_equal(PELTObj$predict(pen = 0.1), seq(50,150,50))
  expect_false(PELTObj$eval(0,51) == 0)
  expect_true(all.equal(PELTObj$eval(0,50), 0))

})


test_that("Active binding `minSize` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a single non-negative numeric/integer
  expect_error(PELT$new(costFunc = costFuncObj, minSize = "a"))
  expect_error(PELT$new(costFunc = costFuncObj, minSize = NULL))
  expect_error(PELT$new(costFunc = costFuncObj, minSize = 1:2))
  expect_error(PELT$new(costFunc = costFuncObj, minSize = 0))

  #Getter
  PELTObj = PELT$new(costFunc = costFuncObj)
  expect_equal(PELTObj$minSize, 1L)

  #Setter
  PELTObj$minSize = 5L
  expect_equal(PELTObj$minSize, 5L)
  expect_error(PELTObj$minSize <- "a")
  expect_error(PELTObj$minSize <- NULL)
  expect_error(PELTObj$minSize <- 1:2)
  expect_error(PELTObj$minSize <- 0)

  #Modifying `minSize` triggers refitting if fitted
  PELTObj$minSize = 1L
  PELTObj$fit(tsMat)
  ms1Seg = PELTObj$predict(0)
  expect_equal(ms1Seg, seq(1, 100, 1))

  PELTObj$minSize = 50L
  ms50Seg = PELTObj$predict(0)
  expect_equal(ms50Seg, c(50, 100))

})

test_that("Active binding `jump` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a single non-negative numeric/integer
  expect_error(PELT$new(costFunc = costFuncObj, jump = "a"))
  expect_error(PELT$new(costFunc = costFuncObj, jump  = NULL))
  expect_error(PELT$new(costFunc = costFuncObj, jump  = 1:2))
  expect_error(PELT$new(costFunc = costFuncObj, jump  = 0))
  expect_no_error(PELT$new())

  #Getter
  PELTObj = PELT$new(costFunc = costFuncObj)
  expect_equal(PELTObj$jump, 1L)

  #Setter
  PELTObj$jump  = 5L
  expect_equal(PELTObj$jump, 5L)
  expect_error(PELTObj$jump <- "a")
  expect_error(PELTObj$jump <- NULL)
  expect_error(PELTObj$jump <- 1:2)
  expect_error(PELTObj$jump <- 0)

  #Modifying `jump` triggers refitting if fitted
  PELTObj$jump  = 1L
  PELTObj$fit(tsMat)
  j1Seg = PELTObj$predict(0)
  expect_equal(j1Seg, seq(1, 100, 1))

  PELTObj$jump  = 5L
  j5Seg = PELTObj$predict(0)
  expect_equal(j5Seg, seq(5, 100, 5))

})

test_that("Active binding `tsMat` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  set.seed(123)
  tsMat2 = matrix(c(rnorm(50,0), rnorm(50,5)))
  tsNA = tsMat
  tsNA[1] = NA
  PELTObj = PELT$new()

  #Wrong input: `tsMat` must be a numeric matrix without any NA
  expect_error(PELTObj$fit(tsNA))
  expect_error(PELTObj$fit(c("a", "b")))
  expect_error(PELTObj$fit(NULL))
  expect_error(PELTObj$fit(as.vector(tsMat)))
  expect_no_error(PELTObj$fit(tsMat))

  #Getter
  PELTObj$fit(tsMat)
  expect_equal(PELTObj$tsMat, tsMat)

  #Setter
  expect_error(PELTObj$tsMat <- tsNA)
  expect_error(PELTObj$tsMat <- c("a", "b"))
  expect_error(PELTObj$tsMat <- NULL)
  expect_error(PELTObj$tsMat <- as.vector(tsMat))

  PELTObj$tsMat = tsMat2
  expect_equal(PELTObj$tsMat, tsMat2)

  #Modify `tsMat` triggers refitting if fitted

  PELTObj$fit(tsMat)
  tsMat1Eval = PELTObj$eval(0,100)
  expect_equal(tsMat1Eval, sum((tsMat - mean(tsMat))^2))

  PELTObj$fit(tsMat2)
  tsMat2Eval = PELTObj$eval(0,100)
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
  PELTObj = PELT$new(costFunc = costFuncObj)

  #Wrong input: `covariate` must be a numeric matrix without any NA
  expect_error(PELTObj$fit(tsMat, as.matrix(covariateMat[-1]))) #nrows do not match
  expect_error(PELTObj$fit(tsMat, covariateNA))
  expect_error(PELTObj$fit(tsMat, c("a", "b")))
  expect_error(PELTObj$fit(tsMat, as.vector(covariateMat)))
  expect_no_error(PELTObj$fit(tsMat, covariateMat))
  expect_no_error(PELTObj$fit(tsMat, NULL))

  #Getter
  PELTObj$fit(tsMat, covariateMat)
  expect_equal(PELTObj$covariates, covariateMat)

  #Setter
  expect_error(PELTObj$covariates <- as.matrix(covariateMat[-1])) #nrows do not match
  expect_error(PELTObj$covariates <- covariateNA)
  expect_error(PELTObj$covariates <- c("a", "b"))
  expect_error(PELTObj$covariates <- NULL)
  expect_error(PELTObj$covariates <- as.vector(covariateMat))

  #Modify `tsMat` triggers refitting if fitted

  set.seed(100)
  covariateMat2 = as.matrix(tsMat/2 + rnorm(100))

  PELTObj$fit(tsMat, covariateMat)
  cM1err = PELTObj$eval(0,100)
  PELTObj$covariates = covariateMat2
  cM2err = PELTObj$eval(0,100)

  expect_equal(cM1err, sum(lm(tsMat~covariateMat)$residuals^2))
  expect_equal(cM2err, sum(lm(tsMat~covariateMat2)$residuals^2))

})

test_that("Active binding `costFunc` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a R6 object of class `costFunc`
  expect_error(PELT$new(costFunc = 1L))
  expect_error(PELT$new(costFunc = list(costFunc = "L2")))
  expect_error(PELT$new(costFunc = NULL))
  expect_no_error(PELT$new())

  #Getter
  PELTObj = PELT$new(costFunc = costFuncObj)
  expect_equal(PELTObj$costFunc$pass()$costFunc, "L2")

  #Setter
  #Wrong input: Must be a R6 object of class `costFunc`
  expect_error(PELTObj$costFunc <- list(costFunc = "L2"))
  expect_error(PELTObj$costFunc <- NULL)
  expect_error(PELTObj$costFunc <- 1L)
  expect_no_error(PELTObj$costFunc <- costFunc$new("VAR"))

  #Modifying `costFunc` triggers refitting if fitted

  PELTObj = PELT$new() #L2
  PELTObj$fit(tsMat)
  expect_equal(PELTObj$eval(0, 100), sum((tsMat - mean(tsMat))^2))

  PELTObj$costFunc = costFunc$new("L1")
  expect_equal(PELTObj$eval(0, 100), sum(abs(tsMat - median(tsMat))))

  PELTObj$costFunc = costFunc$new("SIGMA")
  expect_equal(PELTObj$eval(0, 100), 100*log(det(var(tsMat)*99/100+10^-6)))

})


test_that("Test that `describe()` method works properly)", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")
  PELTObj = PELT$new(costFunc = costFuncObj)
  PELTObj$fit(tsMat)

  #Wrong input
  expect_error(PELTObj$describe(NULL))
  expect_error(PELTObj$describe(1:2))
  expect_error(PELTObj$describe("a"))

  #`describe()` returns the expected outputs
  expect_no_error(PELTObj$describe(T))
  expect_no_error(PELTObj$describe(F))
  expect_no_error(PELTObj$describe())

  expect_equal(PELTObj$describe(T), PELTObj$describe(F))
  expect_equal(PELTObj$describe(), list(minSize = 1L, jump = 1L,
                                          costFunc = costFuncObj, fitted = T,
                                          tsMat = tsMat, covariates = NULL,
                                          n = 100, p = 1))

  #Modify non-`costFunc` active bindings

  set.seed(123)
  tsMat2 = cbind(x1 = c(rnorm(50,0), rnorm(50,5)),
                 x2 = c(rnorm(50,0), rnorm(50,5)))
  PELTObj$tsMat = tsMat2
  PELTObj$minSize = 2L
  PELTObj$jump = 2L

  expect_equal(PELTObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFuncObj, fitted = T,
                                          tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2))

  #Modify `costFunc``
  costFuncObj  = costFunc$new(costFunc = "VAR", pVAR = 1L)
  PELTObj$costFunc = costFuncObj

  expect_equal(PELTObj$describe(T), PELTObj$describe(F))
  expect_equal(PELTObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFuncObj,
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2, pVAR = costFuncObj$pVAR))

  costFuncObj = costFunc$new(costFunc = "SIGMA", pVAR = 1L)
  PELTObj$costFunc = costFuncObj

  expect_equal(PELTObj$describe(T), PELTObj$describe(F))
  expect_equal(PELTObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFunc$new(costFunc = "SIGMA"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2,
                                          addSmallDiag = costFuncObj$addSmallDiag,
                                          epsilon = costFuncObj$epsilon))

  expect_true(is.null(PELTObj$describe()$pVAR))


  costFuncObj = costFunc$new(costFunc = "L1", pVAR = 1L)
  PELTObj$costFunc = costFuncObj

  expect_equal(PELTObj$describe(T), PELTObj$describe(F))
  expect_equal(PELTObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFunc$new(costFunc = "L1"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2))

  expect_true(is.null(PELTObj$describe()$pVAR))
  expect_true(is.null(PELTObj$describe()$epsilon))
  expect_true(is.null(PELTObj$describe()$addSmallDiag))


  costFuncObj = costFunc$new(costFunc = "LinearL2", pVAR = 1L)
  expect_warning(PELTObj$costFunc <- costFuncObj)

  expect_equal(PELTObj$describe(T), PELTObj$describe(F))
  expect_equal(PELTObj$describe(), list(minSize = 2L, jump = 2L,
                                          costFunc = costFunc$new(costFunc = "LinearL2"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2, intercept = TRUE))

  expect_true(is.null(PELTObj$describe()$pVAR))

})

test_that("Error handling for `eval()` works properly", {

  #Cost-specific tests are in other files

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  PELTObj = PELT$new()
  expect_error(PELTObj$eval(0, 10)) #Not fitted

  PELTObj$fit(tsMat)
  expect_error(PELTObj$eval(NULL, 10))
  expect_error(PELTObj$eval(0, NULL))

})

test_that("Test that error handling for `predict()` works properly", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  PELTObj = PELT$new()
  expect_error(PELTObj$predict()) #Not fitted

  PELTObj$fit(tsMat)
  expect_error(PELTObj$predict(NULL))
  expect_error(PELTObj$predict(-1))
  expect_error(PELTObj$predict(c(1:2)))
  expect_error(PELTObj$predict("a"))
  expect_error(PELTObj$predict(T))

})




test_that("Error handling for `plot()` works properly", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  PELTObj = PELT$new()

  #Without providing endPts
  PELTObj$fit(tsMat)
  expect_error(PELTObj$plot()) #No tmpEndPts available
  PELTObj$predict(25)
  expect_no_error(PELTObj$plot())

  ## Invalid `main`
  expect_error(PELTObj$plot(main = c(1:3)))
  expect_error(PELTObj$plot(main = T))
  expect_error(PELTObj$plot(main = NA))
  expect_error(PELTObj$plot(main = NULL))
  expect_error(PELTObj$plot(main = c("a", "b")))

  ## Invalid `xlab`
  expect_error(PELTObj$plot(xlab = c(1:3)))
  expect_error(PELTObj$plot(xlab = T))
  expect_error(PELTObj$plot(xlab = NA))
  expect_error(PELTObj$plot(xlab = NULL))
  expect_error(PELTObj$plot(xlab = c("a", "b")))

  ## Invalid `d`
  expect_error(PELTObj$plot(d = 2)) #d > p
  expect_error(PELTObj$plot(d = "a"))
  expect_error(PELTObj$plot(d = T))
  expect_error(PELTObj$plot(d = NULL))
  expect_error(PELTObj$plot(d = NA))

  ## Invalid `dimNames`
  expect_error(PELTObj$plot(dimNames = c("X1", "X2")))
  expect_error(PELTObj$plot(dimNames = T))
  expect_error(PELTObj$plot(dimNames = NA))
  expect_error(PELTObj$plot(dimNames = NULL))
  expect_error(PELTObj$plot(main = 123))

  #endPts provided
  expect_error(PELTObj$plot(endPts = "a"))
  expect_error(PELTObj$plot(endPts = NA))
  expect_error(PELTObj$plot(endPts = NULL))
  expect_error(PELTObj$plot(endPts = 1:99)) #endPts not include n
  expect_error(PELTObj$plot(endPts = 0:99)) ##min endPts < 1
  expect_error(PELTObj$plot(endPts = 0:101)) ##max endPts > n
  expect_error(PELTObj$plot(endPts = c(50,50,100))) ##duplicated endpts

})


test_that("Error handling for C++ module PELT_L2 works as intended", {
  #constructor: const arma::mat& tsMat, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::PELTCpp_L2, tsMat, 0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_L2, tsMat, 1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_L2, tsMat, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::PELTCpp_L2, tsMat2, 10, 3)) #len = minLen here

})


test_that("Error handling for C++ module PELT_L1 works as intended", {
  #constructor: const arma::mat& tsMat, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::PELTCpp_L1_cwMed, tsMat, 0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_L1_cwMed, tsMat, 1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_L1_cwMed, tsMat, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::PELTCpp_L1_cwMed, tsMat2, 10, 3)) #len = minLen here

})


test_that("Error handling for C++ module PELT_VAR works as intended", {
  #constructor: const arma::mat& tsMat, int pVAR, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  #pVAR = 1
  expect_error(new(rupturesRcpp:::PELTCpp_VAR, tsMat, 1, 0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_VAR, tsMat, 1, 1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_VAR, tsMat, 1, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::PELTCpp_VAR, tsMat2, 1, 10, 3)) #len = minLen here

  set.seed(123)
  tsMat3 = cbind(x1 = rnorm(5), x2 = rnorm(5), x3 = rnorm(5))
  expect_error(new(rupturesRcpp:::PELTCpp_VAR, tsMat3, 2, 1, 1)) #segment too short for fitting VAR(2)
})


test_that("Error handling for C++ module PELT_SIGMA works as intended", {
  #constructor: const arma::mat& tsMat, bool addSmallDiag, double epsilon, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))

  expect_error(new(rupturesRcpp:::PELTCpp_SIGMA, tsMat, T, 10^-6, 0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_SIGMA, tsMat, T, 10^-6, 1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_SIGMA, tsMat, T, 10^-6, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::PELTCpp_SIGMA, tsMat2, T, 10^-6, 10, 3)) #len = minLen here

})


test_that("Error handling for C++ module PELT_LinearL2 works as intended", {
  #constructor: const arma::mat& tsMat,  const arma::mat& covariates, bool intercept_, int minSize_, int jump_

  set.seed(123)
  tsMat = as.matrix(rnorm(23))
  #minLen = 2*jump*ceiling(minSize/jump) = 24 if jump = 3, and minSize = 10
  tsMat2 = as.matrix(rnorm(24))
  covariateMat = as.matrix(rep(1,23))
  covariateMat2 = as.matrix(rep(1,24))

  expect_error(new(rupturesRcpp:::PELTCpp_LinearL2, tsMat, covariateMat, T,  0, 1)) #minSize_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_LinearL2, tsMat, covariateMat, T,  1, 0)) #jump_ = 0
  expect_error(new(rupturesRcpp:::PELTCpp_LinearL2, tsMat, covariateMat, T, 10, 3)) #segment too short
  expect_no_error(new(rupturesRcpp:::PELTCpp_LinearL2, tsMat2, covariateMat2, T, 10, 3)) #len = minLen here

  set.seed(1234)
  tsMat3 = as.matrix(rnorm(5))
  covariateMat3 =cbind(rnorm(5),rnorm(5),rnorm(5),rnorm(5),rnorm(5))
  expect_error(new(rupturesRcpp:::PELTCpp_LinearL2, tsMat3, covariateMat3, T, 1, 1)) #too short for fitting linear regression model
})

test_that("Some additional tests", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  PELTObj = PELT$new()
  PELTObj$fit(tsMat)
  expect_equal(PELTObj$predict(999999), 100) #Too large `pen` -> only return `n`


})

