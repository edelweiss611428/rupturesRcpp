#test-PELT.R

set.seed(12345)

test_that("PELT with L1/L2/SIGMA/VAR works for constant segments", {

  tsMat = matrix(c(rep(0,50), rep(5, 50), rep(10, 50)))

  #L1
  costFuncObj = costFunc$new("L1")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(tsMat)
  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

  #L2
  costFuncObj = costFunc$new("L2")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(tsMat)
  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

  #SIGMA
  costFuncObj = costFunc$new("SIGMA")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(tsMat)
  expect_equal( PeltObj$predict(pen = 0.1), seq(50,150,50))

  #VAR
  costFuncObj = costFunc$new("VAR")
  PeltObj = PELT$new(costFunc = costFuncObj)
  PeltObj$fit(tsMat)

  expect_warning(expect_equal(PeltObj$predict(pen = 0.1), seq(50,150,50)),
                 "Some systems seem singular!")
  expect_warning(expect_false(PeltObj$eval(0,51) == 0), "seems singular")
  expect_warning(expect_true(all.equal(PeltObj$eval(0,50), 0)), "seems singular")

  #LinearL2
  costFuncObj = costFunc$new("LinearL2")
  PELTObj = PELT$new(costFunc = costFuncObj)
  expect_warning(PELTObj$fit(tsMat), "an intercept") #without providing covariate matrix => only intercept

  expect_equal(PELTObj$predict(pen = 0.1), seq(50,150,50))
  expect_false(PELTObj$eval(0,51) == 0)
  expect_true(all.equal(PELTObj$eval(0,50), 0))


})




test_that("Test if active bindings work properly (i.e., setter/getter/input validating abilility)", {

  X1 = matrix(rep(0,100))
  X2 = matrix(1:100)
  #L2
  costFuncObj = costFunc$new("L2")
  PELTObj = PELT$new(costFunc = costFuncObj)
  PELTObj$tsMat = X2
  expect_no_error(PELTObj$fit())
  expect_true(all.equal(PELTObj$tsMat, X2))
  expect_true(all.equal(PELTObj$eval(0,2), 0.5)) #Setting tsMat via active bindings should work
  PELTObj$tsMat = X1
  expect_true(all.equal(PELTObj$eval(0,2), 0)) #Expect refitted after PELTObj$tsMat = X1

  #VAR
  PELTObj$tsMat = X2
  costFuncObj = costFunc$new("VAR")
  PELTObj$costFunc = costFuncObj #Modify cost funcc
  expect_warning(expect_false(PELTObj$eval(0,2)== 0.5), "singular") #No longer use the original cost function

  PELTObj$minSize = 5L
  expect_true(PELTObj$minSize == 5L)

  PELTObj$jump = 5L
  expect_true(PELTObj$jump == 5L)
})


test_that("Expect error when segment is too short", {

  X1 = matrix(rep(0,10))
  costFuncObj = costFunc$new("L2")
  PELTObj = PELT$new(costFunc = costFuncObj, minSize = 10)
  expect_error(PELTObj$fit(X1)) ## nObs must be at least 2*minSize

})



test_that("Expect error when sizes mismatch", {

  X1 = matrix(rnorm(100))
  costFuncObj = costFunc$new("LinearL2")
  PELTObj = PELT$new(costFunc = costFuncObj)
  PELTObj$tsMat = X1
  PELTObj$covariates = matrix(rnorm(99))
  expect_error(PELTObj$fit()) #Sizes not match
  PELTObj$covariates = matrix(rnorm(100))
  expect_no_error(PELTObj$fit())
})

