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


test_that("Test if active bindings work properly (i.e., setter/getter/input validating abilility)", {

  X1 = matrix(rep(0,100))
  X2 = matrix(1:100)
  #L2
  costFuncObj = costFunc$new("L2")
  binSegObj = binSeg$new(costFunc = costFuncObj)
  binSegObj$tsMat = X2
  expect_no_error(binSegObj$fit())
  expect_true(all.equal(binSegObj$tsMat, X2))
  expect_true(all.equal(binSegObj$eval(0,2), 0.5)) #Setting tsMat via active bindings should work
  binSegObj$tsMat = X1
  expect_true(all.equal(binSegObj$eval(0,2), 0)) #Expect refitted after binSegObj$tsMat = X1

  #VAR
  binSegObj$tsMat = X2
  costFuncObj = costFunc$new("VAR")
  suppressWarnings(binSegObj$costFunc <- costFuncObj)
  expect_warning(expect_false(binSegObj$eval(0,2)== 0.5), "singular") #No longer use the original cost function

  binSegObj$minSize = 5L
  expect_true(binSegObj$minSize == 5L)

  binSegObj$jump = 5L
  expect_true(binSegObj$jump == 5L)
})




test_that("Expect error when segment is too short", {

  X1 = matrix(rep(0,10))
  costFuncObj = costFunc$new("L2")
  binSegObj = binSeg$new(costFunc = costFuncObj, minSize = 10)
  expect_error(binSegObj$fit(X1)) ## nObs must be at least 2*minSize

})



test_that("Expect error when sizes mismatch", {

  X1 = matrix(rnorm(100))
  costFuncObj = costFunc$new("LinearL2")
  binSegObj = binSeg$new(costFunc = costFuncObj)
  binSegObj$tsMat = X1
  binSegObj$covariates = matrix(rnorm(99))
  expect_error(binSegObj$fit()) #Sizes not match
  binSegObj$covariates = matrix(rnorm(100))
  expect_no_error(binSegObj$fit())
})

