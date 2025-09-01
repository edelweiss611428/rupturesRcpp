#Window module
set.seed(12345)


test_that("Window_L1 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("L1")
  WindowObj = Window$new(costFunc = costFuncObj)
  WindowObj$fit(tsMat)

  expect_equal( WindowObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("Window_L2 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("L2")
  WindowObj = Window$new(costFunc = costFuncObj)
  WindowObj$fit(tsMat)

  expect_equal( WindowObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("Window_SIGMA works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("SIGMA")
  WindowObj = Window$new(costFunc = costFuncObj)
  WindowObj$fit(tsMat)

  expect_equal( WindowObj$predict(pen = 0.1), seq(50,150,50))

})

test_that("Window_VAR works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("VAR")
  WindowObj = Window$new(costFunc = costFuncObj)

  expect_warning(WindowObj$fit(tsMat), "Some systems seem singular!") #Warning once feature
  expect_equal(WindowObj$predict(pen = 0.1), seq(50,150,50))
  expect_warning(expect_false(WindowObj$eval(0,51) == 0), "seems singular")
  expect_warning(expect_true(all.equal(WindowObj$eval(0,50), 0)), "seems singular")

})


test_that("Window_LinearL2 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))
  costFuncObj = costFunc$new("LinearL2")
  WindowObj = Window$new(costFunc = costFuncObj)

  #when no covariate matrix is provided
  expect_warning(WindowObj$fit(tsMat), "an intercept")
  expect_equal(WindowObj$predict(pen = 0.1), seq(50,150,50))
  expect_false(WindowObj$eval(0,51) == 0)
  expect_true(all.equal(WindowObj$eval(0,50), 0))
})





test_that("Test if active bindings work properly (i.e., setter/getter/input validating abilility)", {

  X1 = matrix(rep(0,100))
  X2 = matrix(1:100)
  #L2
  costFuncObj = costFunc$new("L2")
  WindowObj = Window$new(costFunc = costFuncObj)
  WindowObj$tsMat = X2
  expect_no_error(WindowObj$fit())
  expect_true(all.equal(WindowObj$tsMat, X2))
  expect_true(all.equal(WindowObj$eval(0,2), 0.5)) #Setting tsMat via active bindings should work
  WindowObj$tsMat = X1
  expect_true(all.equal(WindowObj$eval(0,2), 0)) #Expect refitted after WindowObj$tsMat = X1

  #VAR
  WindowObj$tsMat = X2
  costFuncObj = costFunc$new("VAR")
  WindowObj$costFunc = costFuncObj #Modify cost funcc
  expect_warning(expect_false(WindowObj$eval(0,2)== 0.5), "singular") #No longer use the original cost function

  WindowObj$minSize = 5L
  expect_true(WindowObj$minSize == 5L)

  WindowObj$jump = 5L
  expect_true(WindowObj$jump == 5L)
})


test_that("Expect error when segment is too short", {

  X1 = matrix(rep(0,10))
  costFuncObj = costFunc$new("L2")
  WindowObj = Window$new(costFunc = costFuncObj, radius = 25)
  expect_error(WindowObj$fit(X1)) ## too short compared to radius

  WindowObj$radius = 1
  WindowObj$minSize = 10
  expect_error(WindowObj$fit(X1)) ## nObs must be at least 2*minSize

})


test_that("Expect error when sizes mismatch", {

  X1 = matrix(rnorm(100))
  costFuncObj = costFunc$new("LinearL2")
  WindowObj = Window$new(costFunc = costFuncObj, radius = 25)
  WindowObj$tsMat = X1
  WindowObj$covariates = matrix(rnorm(99))
  expect_error(WindowObj$fit()) #Sizes not match
  WindowObj$covariates = matrix(rnorm(100))
  expect_no_error(WindowObj$fit())
})



test_that("Active binding `minSize` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5), rnorm(50,-15)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a single non-negative numeric/integer
  expect_error(Window$new(costFunc = costFuncObj, minSize = "a"))
  expect_error(Window$new(costFunc = costFuncObj, minSize = NULL))
  expect_error(Window$new(costFunc = costFuncObj, minSize = 1:2))
  expect_error(Window$new(costFunc = costFuncObj, minSize = 0))

  #Getter
  WindowObj = Window$new(costFunc = costFuncObj)
  expect_equal(WindowObj$minSize, 1L)

  #Setter
  WindowObj$minSize = 5L
  expect_equal(WindowObj$minSize, 5L)
  expect_error(WindowObj$minSize <- "a")
  expect_error(WindowObj$minSize <- NULL)
  expect_error(WindowObj$minSize <- 1:2)
  expect_error(WindowObj$minSize <- 0)

  #Modifying `minSize` triggers refitting if fitted
  WindowObj$minSize = 1L
  WindowObj$fit(tsMat)
  ms1Seg = WindowObj$predict(0)

  WindowObj$minSize = 75L
  ms75Seg = WindowObj$predict(0)

  expect_false(identical(ms1Seg, ms75Seg))

})


test_that("Active binding `jump` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a single non-negative numeric/integer
  expect_error(Window$new(costFunc = costFuncObj, jump = "a"))
  expect_error(Window$new(costFunc = costFuncObj, jump  = NULL))
  expect_error(Window$new(costFunc = costFuncObj, jump  = 1:2))
  expect_error(Window$new(costFunc = costFuncObj, jump  = 0))
  expect_no_error(Window$new())

  #Getter
  WindowObj = Window$new(costFunc = costFuncObj)
  expect_equal(WindowObj$jump, 1L)

  #Setter
  WindowObj$jump  = 5L
  expect_equal(WindowObj$jump, 5L)
  expect_error(WindowObj$jump <- "a")
  expect_error(WindowObj$jump <- NULL)
  expect_error(WindowObj$jump <- 1:2)
  expect_error(WindowObj$jump <- 0)

  #Modifying `jump` triggers refitting if fitted
  WindowObj$jump  = 1L
  WindowObj$fit(tsMat)
  j1Seg = WindowObj$predict(0)

  WindowObj$jump  = 3L
  j5Seg = WindowObj$predict(0)

  expect_false(identical(j1Seg, j5Seg))

})

# Lack of testing for active binding `radius`

test_that("Active binding `tsMat` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  set.seed(123)
  tsMat2 = matrix(c(rnorm(50,0), rnorm(50,5)))
  tsNA = tsMat
  tsNA[1] = NA
  WindowObj = Window$new()

  #Wrong input: `tsMat` must be a numeric matrix without any NA
  expect_error(WindowObj$fit(tsNA))
  expect_error(WindowObj$fit(c("a", "b")))
  expect_error(WindowObj$fit(NULL))
  expect_error(WindowObj$fit(as.vector(tsMat)))
  expect_no_error(WindowObj$fit(tsMat))

  #Getter
  WindowObj$fit(tsMat)
  expect_equal(WindowObj$tsMat, tsMat)

  #Setter
  expect_error(WindowObj$tsMat <- tsNA)
  expect_error(WindowObj$tsMat <- c("a", "b"))
  expect_error(WindowObj$tsMat <- NULL)
  expect_error(WindowObj$tsMat <- as.vector(tsMat))

  WindowObj$tsMat = tsMat2
  expect_equal(WindowObj$tsMat, tsMat2)

  #Modify `tsMat` triggers refitting if fitted

  WindowObj$fit(tsMat)
  tsMat1Eval = WindowObj$eval(0,100)
  expect_equal(tsMat1Eval, sum((tsMat - mean(tsMat))^2))

  WindowObj$fit(tsMat2)
  tsMat2Eval = WindowObj$eval(0,100)
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
  WindowObj = Window$new(costFunc = costFuncObj)

  #Wrong input: `covariate` must be a numeric matrix without any NA
  expect_error(WindowObj$fit(tsMat, as.matrix(covariateMat[-1]))) #nrows do not match
  expect_error(WindowObj$fit(tsMat, covariateNA))
  expect_error(WindowObj$fit(tsMat, c("a", "b")))
  expect_error(WindowObj$fit(tsMat, as.vector(covariateMat)))
  expect_no_error(WindowObj$fit(tsMat, covariateMat))
  expect_no_error(WindowObj$fit(tsMat, NULL))

  #Getter
  WindowObj$fit(tsMat, covariateMat)
  expect_equal(WindowObj$covariates, covariateMat)

  #Setter
  expect_error(WindowObj$covariates <- as.matrix(covariateMat[-1])) #nrows do not match
  expect_error(WindowObj$covariates <- covariateNA)
  expect_error(WindowObj$covariates <- c("a", "b"))
  expect_error(WindowObj$covariates <- NULL)
  expect_error(WindowObj$covariates <- as.vector(covariateMat))

  #Modify `tsMat` triggers refitting if fitted

  set.seed(100)
  covariateMat2 = as.matrix(tsMat/2 + rnorm(100))

  WindowObj$fit(tsMat, covariateMat)
  cM1err = WindowObj$eval(0,100)
  WindowObj$covariates = covariateMat2
  cM2err = WindowObj$eval(0,100)

  expect_equal(cM1err, sum(lm(tsMat~covariateMat)$residuals^2))
  expect_equal(cM2err, sum(lm(tsMat~covariateMat2)$residuals^2))

})

test_that("Active binding `costFunc` works as intended", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")

  #Wrong input: Must be a R6 object of class `costFunc`
  expect_error(Window$new(costFunc = 1L))
  expect_error(Window$new(costFunc = list(costFunc = "L2")))
  expect_error(Window$new(costFunc = NULL))
  expect_no_error(Window$new())

  #Getter
  WindowObj = Window$new(costFunc = costFuncObj)
  expect_equal(WindowObj$costFunc$pass()$costFunc, "L2")

  #Setter
  #Wrong input: Must be a R6 object of class `costFunc`
  expect_error(WindowObj$costFunc <- list(costFunc = "L2"))
  expect_error(WindowObj$costFunc <- NULL)
  expect_error(WindowObj$costFunc <- 1L)
  expect_no_error(WindowObj$costFunc <- costFunc$new("VAR"))

  #Modifying `costFunc` triggers refitting if fitted

  WindowObj = Window$new() #L2
  WindowObj$fit(tsMat)
  expect_equal(WindowObj$eval(0, 100), sum((tsMat - mean(tsMat))^2))

  WindowObj$costFunc = costFunc$new("L1")
  expect_equal(WindowObj$eval(0, 100), sum(abs(tsMat - median(tsMat))))

  WindowObj$costFunc = costFunc$new("SIGMA")
  expect_equal(WindowObj$eval(0, 100), 100*log(det(var(tsMat)*99/100+10^-6)))

})


test_that("Test that `describe()` method works properly)", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  costFuncObj = costFunc$new("L2")
  WindowObj = Window$new(costFunc = costFuncObj)
  WindowObj$fit(tsMat)

  #Wrong input
  expect_error(WindowObj$describe(NULL))
  expect_error(WindowObj$describe(1:2))
  expect_error(WindowObj$describe("a"))

  #`describe()` returns the expected outputs
  expect_no_error(WindowObj$describe(T))
  expect_no_error(WindowObj$describe(F))
  expect_no_error(WindowObj$describe())

  expect_equal(WindowObj$describe(T), WindowObj$describe(F))
  expect_equal(WindowObj$describe(), list(minSize = 1L, jump = 1L, radius = 25,
                                          costFunc = costFuncObj, fitted = T,
                                          tsMat = tsMat, covariates = NULL,
                                          n = 100, p = 1))

  #Modify non-`costFunc` active bindings

  set.seed(123)
  tsMat2 = cbind(x1 = c(rnorm(50,0), rnorm(50,5)),
                 x2 = c(rnorm(50,0), rnorm(50,5)))
  WindowObj$tsMat = tsMat2
  WindowObj$minSize = 2L
  WindowObj$jump = 2L
  WindowObj$radius = 10L;
  expect_equal(WindowObj$describe(), list(minSize = 2L, jump = 2L, radius = 10L,
                                          costFunc = costFuncObj, fitted = T,
                                          tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2))

  #Modify `costFunc``
  costFuncObj  = costFunc$new(costFunc = "VAR", pVAR = 1L)
  WindowObj$costFunc = costFuncObj

  expect_equal(WindowObj$describe(T), WindowObj$describe(F))
  expect_equal(WindowObj$describe(), list(minSize = 2L, jump = 2L,  radius = 10L,
                                          costFunc = costFuncObj,
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2, pVAR = costFuncObj$pVAR))

  costFuncObj = costFunc$new(costFunc = "SIGMA", pVAR = 1L)
  WindowObj$costFunc = costFuncObj

  expect_equal(WindowObj$describe(T), WindowObj$describe(F))
  expect_equal(WindowObj$describe(), list(minSize = 2L, jump = 2L, radius = 10L,
                                          costFunc = costFunc$new(costFunc = "SIGMA"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2,
                                          addSmallDiag = costFuncObj$addSmallDiag,
                                          epsilon = costFuncObj$epsilon))

  expect_true(is.null(WindowObj$describe()$pVAR))

  costFuncObj = costFunc$new(costFunc = "L1", pVAR = 1L)
  WindowObj$costFunc = costFuncObj

  expect_equal(WindowObj$describe(T), WindowObj$describe(F))
  expect_equal(WindowObj$describe(), list(minSize = 2L, jump = 2L, radius = 10L,
                                          costFunc = costFunc$new(costFunc = "L1"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2))

  expect_true(is.null(WindowObj$describe()$pVAR))
  expect_true(is.null(WindowObj$describe()$epsilon))
  expect_true(is.null(WindowObj$describe()$addSmallDiag))


  costFuncObj = costFunc$new(costFunc = "LinearL2", pVAR = 1L)
  expect_warning(WindowObj$costFunc <- costFuncObj)

  expect_equal(WindowObj$describe(T), WindowObj$describe(F))
  expect_equal(WindowObj$describe(), list(minSize = 2L, jump = 2L, radius = 10L,
                                          costFunc = costFunc$new(costFunc = "LinearL2"),
                                          fitted = T, tsMat = tsMat2, covariates = NULL,
                                          n = 100, p = 2, intercept = TRUE))

  expect_true(is.null(WindowObj$describe()$pVAR))

})



test_that("Error handling for `eval()` works properly", {

  #Cost-specific tests are in other files

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  WindowObj = Window$new()
  expect_error(WindowObj$eval(0, 10)) #Not fitted

  WindowObj$fit(tsMat)
  expect_error(WindowObj$eval(NULL, 10))
  expect_error(WindowObj$eval(0, NULL))

})

test_that("Test that error handling for `predict()` works properly", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  WindowObj = Window$new()
  expect_error(WindowObj$predict()) #Not fitted

  WindowObj$fit(tsMat)
  expect_error(WindowObj$predict(NULL))
  expect_error(WindowObj$predict(-1))
  expect_error(WindowObj$predict(c(1:2)))
  expect_error(WindowObj$predict("a"))
  expect_error(WindowObj$predict(T))

})




test_that("Error handling for `plot()` works properly", {

  set.seed(12345)
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  WindowObj = Window$new()

  #Without providing endPts
  WindowObj$fit(tsMat)
  expect_error(WindowObj$plot()) #No tmpEndPts available
  WindowObj$predict(25)
  expect_no_error(WindowObj$plot())

  ## Invalid `main`
  expect_error(WindowObj$plot(main = c(1:3)))
  expect_error(WindowObj$plot(main = T))
  expect_error(WindowObj$plot(main = NA))
  expect_error(WindowObj$plot(main = NULL))
  expect_error(WindowObj$plot(main = c("a", "b")))

  ## Invalid `xlab`
  expect_error(WindowObj$plot(xlab = c(1:3)))
  expect_error(WindowObj$plot(xlab = T))
  expect_error(WindowObj$plot(xlab = NA))
  expect_error(WindowObj$plot(xlab = NULL))
  expect_error(WindowObj$plot(xlab = c("a", "b")))

  ## Invalid `d`
  expect_error(WindowObj$plot(d = 2)) #d > p
  expect_error(WindowObj$plot(d = "a"))
  expect_error(WindowObj$plot(d = T))
  expect_error(WindowObj$plot(d = NULL))
  expect_error(WindowObj$plot(d = NA))

  ## Invalid `dimNames`
  expect_error(WindowObj$plot(dimNames = c("X1", "X2")))
  expect_error(WindowObj$plot(dimNames = T))
  expect_error(WindowObj$plot(dimNames = NA))
  expect_error(WindowObj$plot(dimNames = NULL))
  expect_error(WindowObj$plot(main = 123))

  #endPts provided
  expect_error(WindowObj$plot(endPts = "a"))
  expect_error(WindowObj$plot(endPts = NA))
  expect_error(WindowObj$plot(endPts = NULL))
  expect_error(WindowObj$plot(endPts = 1:99)) #endPts not include n
  expect_error(WindowObj$plot(endPts = 0:99)) ##min endPts < 1
  expect_error(WindowObj$plot(endPts = 0:101)) ##max endPts > n
  expect_error(WindowObj$plot(endPts = c(50,50,100))) ##duplicated endpts

})


test_that("`Window` does not collapse to 0 change point if `pen` is too large", {

  set.seed(12345)
  #There is 1 change point. Given the default settings, Window should pick up the change point
  tsMat = matrix(c(rnorm(50,0), rnorm(50,5)))
  WindowObj = Window$new()
  WindowObj$fit(tsMat)
  #Does not collapse as `Window` does not consider the situation where there is no changepoint
  expect_false(identical(WindowObj$predict(999999), 100))


})

#Lack of testing for C++ SW modules (to-be-updated)
#Lack of comparision to gold-standard Slicing Window here (To-be-updated); however, testing has been done before

