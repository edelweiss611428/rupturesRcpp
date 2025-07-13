#test-Window.R

set.seed(12345)
test_that("Window with L2/SIGMA/VAR/LinearL2 works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))

  #L2
  costFuncObj = costFunc$new("L2")
  WindowObj = Window$new(costFunc = costFuncObj)
  WindowObj$fit(tsMat)
  expect_equal( WindowObj$predict(pen = 0.1), seq(50,150,50))

  #SIGMA
  costFuncObj = costFunc$new("SIGMA")
  WindowObj = Window$new(costFunc = costFuncObj)
  WindowObj$fit(tsMat)
  expect_equal( WindowObj$predict(pen = 0.1), seq(50,150,50))

  #VAR
  costFuncObj = costFunc$new("VAR")
  WindowObj = Window$new(costFunc = costFuncObj)
  expect_warning(WindowObj$fit(tsMat))

  expect_equal(WindowObj$predict(pen = 0.1), seq(50,150,50))
  expect_warning(expect_false(WindowObj$eval(0,51) == 0), "singular") #singular as Z^TZ where Z is stacked lag matrix is not invertible
  expect_warning(expect_true(all.equal(WindowObj$eval(0,50), 0)),  "singular") #singular as Z^TZ where Z is stacked lag matrix is not invertible

  #LinearL2
  costFuncObj = costFunc$new("LinearL2")
  WindowObj = Window$new(costFunc = costFuncObj)
  expect_warning(WindowObj$fit(tsMat), "an intercept") #without providing covariate matrix => only intercept

  expect_equal(WindowObj$predict(pen = 0.1), seq(50,150,50))
  expect_false(WindowObj$eval(0,51) == 0)
  expect_true(all.equal(WindowObj$eval(0,50), 0))

})





test_that("Test if active bindings work properl (i.e., setter/getter/input validating abilility)", {

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

  WindowObj$minSize = 5
  expect_warning(WindowObj$fit(X1), "Diameter") ##  Diameter should be at least `minSize`

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
